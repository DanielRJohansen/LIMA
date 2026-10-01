#include "MDFiles.h" 
#include "Filehandling.h"

#
#include <algorithm>
#include <format>
#include <ranges>
#include "TimeIt.h"
#include "ParallelFor.h"
#include <execution>
#include <atomic>
#include <span>













using namespace FileUtils;

class TopologySectionGetter {
	int dihedralCount = 0;
	int dihedraltypesCount = 0;

	void Reset() {
		dihedralCount = 0;
		dihedraltypesCount = 0;
	}

public:
	constexpr TopologySection operator()(const std::string_view& directive) {
		if (directive == "molecules") return molecules;
		if (directive == "moleculetype") {
			Reset();
			return moleculetype;
		}
		if (directive == "atoms") return atoms;
		if (directive == "bonds") return bonds;
		if (directive == "pairs") return pairs;
		if (directive == "angles") return angles;
		if (directive == "dihedrals") {
			if (++dihedralCount <= 2) return dihedralCount == 1 ? dihedrals : impropers;
			throw std::runtime_error("Encountered 'dihedral' directive more than 2 times in .itp/.top file");
		}
		if (directive == "position_restraints") return position_restraints;
		if (directive == "system") return _system;
		if (directive == "cmap") return cmap;

		if (directive == "atomtypes") return atomtypes;
		if (directive == "pairtypes") return pairtypes;
		if (directive == "bondtypes") return bondtypes;
		if (directive == "constrainttypes") return constainttypes;
		if (directive == "angletypes") return angletypes;
		if (directive == "dihedraltypes") {
			if (++dihedraltypesCount <= 2) return dihedraltypesCount == 1 ? dihedraltypes : impropertypes;
			throw std::runtime_error("Encountered 'dihedraltypes' directive more than 2 times in .itp/ file");
		}
		if (directive == "impropertypes") return impropertypes;
		if (directive == "defaults") return defaults;
		if (directive == "cmaptypes") return cmaptypes;


		if (directive == "settles") return notimplemented;
		if (directive == "exclusions") return notimplemented;
		if (directive == "nonbond_params") return notimplemented; // TODO implement this: https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html

		throw std::runtime_error(std::format("Got unexpected topology directive: {}", directive));
	}
};

constexpr std::string extractSectionName(std::string line) {
	size_t start = line.find('[');
	size_t end = line.find(']', start);
	if (start != std::string::npos && end != std::string::npos) {
		// Extract the view of the text between '[' and ']'
		std::string sectionView(line.c_str() + start + 1, end - start - 1);

		// Remove whitespace
		size_t firstNonWhitespace = sectionView.find_first_not_of(" \t\n\r");
		size_t lastNonWhitespace = sectionView.find_last_not_of(" \t\n\r");
		if (firstNonWhitespace == std::string_view::npos) {
			return {}; // The section is all whitespace
		}
		return sectionView.substr(firstNonWhitespace, lastNonWhitespace - firstNonWhitespace + 1);
	}
	return {}; // Return an empty string_view
}


inline std::string GetCleanFilename(const fs::path& path) {
	auto filename = path.stem().string();
	const std::string prefix = "topol_";
	return filename.starts_with(prefix) ? filename.substr(prefix.length()) : filename;
}

inline std::optional<fs::path> _SearchForFile(const fs::path& dir, const std::string& filename) {
	const std::array<std::string, 2> extensions = { ".itp", ".top" };
	const std::array<std::string, 2> prefixes = { std::string(""), "topol_" };

	for (const auto& ext : extensions) {
		for (const auto& prefix : prefixes) {
			fs::path path = dir / (prefix + filename + ext);
			if (fs::exists(path)) {
				return path;
			}
		}
	}

	return std::nullopt;
}



/// <returns>True if we change section/stay on no section, in which case we should 'continue' to the next line. False otherwise</returns>
constexpr bool HandleTopologySectionStartAndStop(std::string_view line, TopologySection& currentSection, TopologySectionGetter& sectionGetter) {

	if (line.empty() && currentSection == TopologySection::title) {
		currentSection = no_section;
		return true;
	}
	else if (!line.empty() && line[0] == '[') {
		currentSection = sectionGetter(extractSectionName(std::string{ line }));
		return true;
	}

	return false; // We did not change section, and should continue parsing the current line
}
constexpr bool isOnlySpacesAndTabs(std::string_view str) {
	return std::all_of(str.begin(), str.end(), [](char c) {
		return c == ' ' || c == '\t' || c == '\r';
		});
}

// Plain loops instead of find_first_not_of, which is surprisingly slow here, and these run for every value in a topology
inline void SkipSpacesAndTabs(std::string_view& sv) noexcept {
	size_t i = 0;
	while (i < sv.size() && (sv[i] == ' ' || sv[i] == '\t'))
		i++;
	sv.remove_prefix(i);
}

inline void SkipLeadingWhitespace(std::string_view& sv) noexcept {
	// ASCII whitespace incl. CR/LF/TAB/VT/FF
	size_t i = 0;
	while (i < sv.size() && (sv[i] == ' ' || (sv[i] >= '\t' && sv[i] <= '\r')))
		i++;
	sv.remove_prefix(i);
}

// This is safe, even if there is no value present in sv
// for float, int, returns false if err
template<typename T>
inline bool ParseValue(std::string_view& sv, T& out) noexcept {
	if (sv.empty()) 
		return false;

	const char* begin = sv.data();
	const char* end = begin + sv.size();
	auto res = std::from_chars(begin, end, out);
	if (res.ec != std::errc{}) 
		return false;

	sv.remove_prefix(res.ptr - begin);
	SkipSpacesAndTabs(sv);
	return true;
}

inline bool ParseValue(std::string_view& sv, std::optional<float>& out) noexcept {
	float tmp{};	
	if (!ParseValue(sv, tmp))
		return false;
	out = tmp;
	return true;
}

inline bool ParseValue(std::string_view& sv, std::string& out) noexcept {
	if (sv.empty()) return false;

	const auto sepPos = sv.find_first_of(" \t");
	if (sepPos == 0) return false;                 // token may not start with space

	const auto tokLen = (sepPos == std::string_view::npos) ? sv.size() : sepPos;
	out = sv.substr(0, tokLen);
	sv.remove_prefix(tokLen);

	SkipSpacesAndTabs(sv);	// skip trailing spaces after the token
	return !out.empty();
}


// ------------------------------------------ Moleculetype entries ------------------------------------------ //

namespace {
	// Indexed by the 1-indexed gro id from the file, gives the 0-indexed LIMA id, or -1 if no such atom exists
	using GroIdToLimaId = std::span<const int>;

	// Returns the gro id of the atom
	int ParseAtomsEntry(std::string_view sv, TopologyFile::AtomsEntry& atom, int index /*relative to moleculetype*/) {
		SkipLeadingWhitespace(sv);

		int groId = -1;
		if (!ParseValue<int>(sv, groId) || groId < 0)
			throw std::runtime_error(std::format("Failed to read atom id in line: {}", sv));
		ParseValue(sv, atom.type);
		ParseValue<int>(sv, atom.resnr);
		ParseValue(sv, atom.residue);
		ParseValue(sv, atom.atomname);
		ParseValue(sv, atom.cgnr);
		ParseValue(sv, atom.charge);
		ParseValue(sv, atom.mass);	// This might not be present

		atom.id = index;

		if (atom.type.empty() || atom.residue.empty() || atom.atomname.empty())
			throw std::runtime_error("Atom type, residue or atomname is empty");
		return groId;
	}

	template <size_t n>
	bool LoadIds(std::string_view& sv, std::array<int, n>& ids, GroIdToLimaId groIdToLimaId) {
		for (auto& id : ids) {
			int groId;
			if (!ParseValue<int>(sv, groId) || groId < 0 || groId >= static_cast<int>(groIdToLimaId.size()) || groIdToLimaId[groId] == -1)
				return false;
			id = groIdToLimaId[groId];
		}
		return true;
	}

	// Each ParseBond returns false if the bond references atoms that don't exist, in which case the bond must be discarded
	// The parameters are optional. Some top files have them, some dont. If they exist, they take precedence over the forcefield
	bool ParseBond(std::string_view sv, TopologyFile::SingleBond& bond, GroIdToLimaId groIdToLimaId) {
		SkipLeadingWhitespace(sv);
		if (!LoadIds(sv, bond.ids, groIdToLimaId))
			return false;

		float b0, kb;
		bool err = false;
		err |= !ParseValue<int>(sv, bond.funct);
		err |= !ParseValue<float>(sv, b0);
		err |= !ParseValue<float>(sv, kb);
		if (!err)
			bond.parameters = Bondtypes::SingleBond::Parameters::CreateFromCharmm(b0, kb);
		return true;
	}

	bool ParseBond(std::string_view sv, TopologyFile::PairBond& bond, GroIdToLimaId groIdToLimaId) {
		SkipLeadingWhitespace(sv);
		if (!LoadIds(sv, bond.ids, groIdToLimaId))
			return false;

		float sigma, epsilon;
		bool err = false;
		err |= !ParseValue<int>(sv, bond.funct);
		err |= !ParseValue<float>(sv, sigma);
		err |= !ParseValue<float>(sv, epsilon);
		if (!err)
			bond.parameters = Bondtypes::PairBond::Parameters::CreateFromCharmm(sigma, epsilon);
		return true;
	}

	bool ParseBond(std::string_view sv, TopologyFile::AngleBond& bond, GroIdToLimaId groIdToLimaId) {
		SkipLeadingWhitespace(sv);
		if (!LoadIds(sv, bond.ids, groIdToLimaId))
			return false;

		float theta0, ktheta, ub0, kUb;
		bool err = false;
		err |= !ParseValue<int>(sv, bond.funct);
		err |= !ParseValue<float>(sv, theta0);
		err |= !ParseValue<float>(sv, ktheta);
		err |= !ParseValue<float>(sv, ub0);
		err |= !ParseValue<float>(sv, kUb);
		if (!err)
			bond.parameters = Bondtypes::AngleUreyBradleyBond::Parameters::CreateFromCharmm(theta0, ktheta, ub0, kUb, bond.funct);
		return true;
	}

	bool ParseBond(std::string_view sv, TopologyFile::DihedralBond& bond, GroIdToLimaId groIdToLimaId) {
		SkipLeadingWhitespace(sv);
		if (!LoadIds(sv, bond.ids, groIdToLimaId))
			return false;

		bool err = false;
		float phi0, kphi; int n;
		err |= !ParseValue<int>(sv, bond.funct);
		err |= !ParseValue<float>(sv, phi0);
		err |= !ParseValue<float>(sv, kphi);
		err |= !ParseValue<int>(sv, n);
		if (!err)
			bond.parameters = Bondtypes::DihedralBond::Parameters::CreateFromCharmm(phi0, kphi, n);
		return true;
	}

	bool ParseBond(std::string_view sv, TopologyFile::ImproperDihedralBond& bond, GroIdToLimaId groIdToLimaId) {
		SkipLeadingWhitespace(sv);
		if (!LoadIds(sv, bond.ids, groIdToLimaId))
			return false;

		float phi0, kphi;
		bool err = false;
		err |= !ParseValue<int>(sv, bond.funct);
		err |= !ParseValue<float>(sv, phi0);
		err |= !ParseValue<float>(sv, kphi);
		if (!err)
			bond.parameters = Bondtypes::ImproperDihedralBond::Parameters::CreateFromCharmm(phi0, kphi);
		return true;
	}

	bool ParseBond(std::string_view sv, TopologyFile::CmapBond& bond, GroIdToLimaId groIdToLimaId) {
		SkipLeadingWhitespace(sv);
		if (!LoadIds(sv, bond.ids, groIdToLimaId))
			return false;
		ParseValue<int>(sv, bond.funct);
		return true;
	}


	// Calls f(line) for each line in text, without the trailing '\n' or "\r\n"
	template <typename F>
	void ForEachLine(std::string_view text, F&& f) {
		while (!text.empty()) {
			const size_t newline = text.find('\n');
			std::string_view line = text.substr(0, newline);
			text.remove_prefix(newline == std::string_view::npos ? text.size() : newline + 1);
			if (!line.empty() && line.back() == '\r')
				line.remove_suffix(1);
			f(line);
		}
	}

	// A line with something other than whitespace or a comment
	bool IsDataLine(std::string_view line) {
		SkipLeadingWhitespace(line);
		return !line.empty() && line[0] != ';';
	}

	// Parsing is parallel per moleculetype, so only large moleculetypes are worth parallelizing internally
	template <typename F>
	void ForEachEntry(size_t n, F&& f) {
		ParallelUtils::ParallelForBlocked(n, f);
	}

	template <typename Bond>
	void ParseBonds(const std::vector<std::string_view>& lines, std::vector<Bond>& bonds, GroIdToLimaId groIdToLimaId) {
		bonds.resize(lines.size());
		std::vector<uint8_t> valid(lines.size());
		ForEachEntry(lines.size(), [&](size_t i) { valid[i] = ParseBond(lines[i], bonds[i], groIdToLimaId); });

		// Discard bonds that reference atoms which don't exist
		size_t nValid = 0;
		for (size_t i = 0; i < bonds.size(); i++) {
			if (!valid[i]) continue;
			if (nValid != i) bonds[nValid] = bonds[i];
			nValid++;
		}
		bonds.resize(nValid);
	}
}


// ------------------------------------------ Topology file parsing ------------------------------------------ //

// Parses a topology and all its includes in 3 phases:
//	1. Scan (parallel per file): read every reachable file, and locate its preprocessor directives and [ section ] headers.
//	   This does not depend on the #defines, so we scan the includes of both branches of an #ifdef
//	2. Resolve (sequential): walk the directives in include order with a single set of #defines, like the C preprocessor does.
//	   This decides which parts of the files are active, and is cheap since it only visits the directives
//	3. Parse: walk the active text in order to build the moleculetypes, system and forcefield (sequential, only visits headers
//	   and the small sections), then parse the atoms and bonds of all moleculetypes (parallel per moleculetype)
// As before, each included file starts with its own section state, and a moleculetype can not span multiple files
class TopologyParser {
	struct Directive {
		enum class Kind { Define, Undef, Ifdef, Ifndef, If, Elif, Else, Endif, Include, Other };
		Kind kind = Kind::Other;
		size_t lineBegin = 0, lineEnd = 0;	// lineEnd is past the newline
		std::string argument;				// The symbol of a define/undef/ifdef/ifndef, the filename of an include

		// Only for includes
		bool ignoredInclude = false;
		std::optional<fs::path> includePath;	// nullopt if the file could not be found
		int includeFile = -1;
	};
	struct Header {
		size_t lineBegin = 0, lineEnd = 0;	// lineEnd is past the newline
		std::string name;
	};
	struct SourceFile {
		fs::path path;
		std::string contents;
		std::vector<Directive> directives;
		std::vector<Header> headers;
		std::atomic<int> unparsedMoleculetypes = 0; // Contents can be released when this reaches 0
	};

	// An active text range of the file, or an included file
	struct Piece {
		size_t begin = 0, end = 0;
		int childInstance = -1;
	};
	// A file may be included multiple times, and different parts may be active each time
	struct Instance {
		int file = -1;
		std::optional<fs::path> includeName;	// nullopt for the main file
		std::vector<Piece> pieces;
	};

	// The atoms/bonds/... sections of a moleculetype, which we parse in parallel in the last phase
	static constexpr std::array entrySections{ atoms, bonds, pairs, angles, dihedrals, impropers, cmap };
	struct MoleculetypeWork {
		std::shared_ptr<TopologyFile::Moleculetype> moleculetype;
		int file = -1;
		std::array<std::vector<std::string_view>, entrySections.size()> chunks;
	};

	struct SectionState {
		TopologySection section = TopologySection::title;
		TopologySectionGetter sectionGetter{};
		int work = -1;	// Index of the most recent moleculetype in this file
	};

	TopologyFile& topology;
	std::vector<std::unique_ptr<SourceFile>> files;
	std::unordered_map<std::string, int> fileIndices;
	std::vector<Instance> instances;
	std::vector<MoleculetypeWork> works;

public:
	TopologyParser(TopologyFile& topology) : topology(topology) {}

	void Parse(const fs::path& path) {
		topology.defines.insert("FLEXIBLE");// Cant handle gromacs definition of rigid water right now

		ScanAll(path);
		const int mainInstance = Resolve(0, std::nullopt, 0);
		Walk(mainInstance);

		for (auto& file : files)
			if (file->unparsedMoleculetypes == 0)
				ReleaseContents(*file);

		ParallelUtils::ParallelFor(works.size(), [&](size_t i) {
			ParseMoleculetype(works[i]);
			if (--files[works[i].file]->unparsedMoleculetypes == 0)
				ReleaseContents(*files[works[i].file]);
			});
	}

private:
	// ------------------------------------------ Phase 1: Scan ------------------------------------------ //

	int AddFile(const fs::path& path) {
		const std::string key = path.lexically_normal().string();
		if (auto it = fileIndices.find(key); it != fileIndices.end())
			return it->second;

		files.emplace_back(std::make_unique<SourceFile>());
		files.back()->path = path;
		fileIndices.insert({ key, static_cast<int>(files.size() - 1) });
		return static_cast<int>(files.size() - 1);
	}

	// Scans files breadth first, as we only learn of the includes of a file after scanning it
	void ScanAll(const fs::path& mainFile) {
		AddFile(mainFile);
		for (size_t scanned = 0; scanned < files.size();) {
			const size_t end = files.size();
			ParallelUtils::ParallelFor(end - scanned, [&](size_t i) { Scan(*files[scanned + i]); });

			for (size_t i = scanned; i < end; i++)
				for (auto& directive : files[i]->directives)
					if (directive.includePath)
						directive.includeFile = AddFile(*directive.includePath);
			scanned = end;
		}
	}

	static void ReleaseContents(SourceFile& file) {
		std::string{}.swap(file.contents);
	}

	static std::string_view FirstToken(std::string_view sv) {
		SkipLeadingWhitespace(sv);
		return sv.substr(0, sv.find_first_of(" \t;"));
	}

	static Directive ParseDirective(std::string_view line, const fs::path& dir) {
		// line starts at the '#'
		line.remove_prefix(1);
		SkipLeadingWhitespace(line);
		const size_t wordEnd = line.find_first_of(" \t;\"<");
		const std::string_view word = line.substr(0, wordEnd);
		const std::string_view rest = wordEnd == std::string_view::npos ? std::string_view{} : line.substr(wordEnd);

		using Kind = Directive::Kind;
		Directive directive;
		if (word == "define") { directive.kind = Kind::Define; directive.argument = FirstToken(rest); }
		else if (word == "undef") { directive.kind = Kind::Undef; directive.argument = FirstToken(rest); }
		else if (word == "ifdef") { directive.kind = Kind::Ifdef; directive.argument = FirstToken(rest); }
		else if (word == "ifndef") { directive.kind = Kind::Ifndef; directive.argument = FirstToken(rest); }
		else if (word == "if") { directive.kind = Kind::If; }
		else if (word == "elif") { directive.kind = Kind::Elif; }
		else if (word == "else") { directive.kind = Kind::Else; }
		else if (word == "endif") { directive.kind = Kind::Endif; }
		else if (word == "include") {
			directive.kind = Kind::Include;
			const size_t open = rest.find_first_of("\"<");
			const size_t close = open == std::string_view::npos ? open : rest.find_first_of("\">", open + 1);
			if (close == std::string_view::npos || close == open + 1)
				throw std::runtime_error(std::format("Include is not formatted as expected: {}", line));
			directive.argument = rest.substr(open + 1, close - open - 1);

			const std::string& filename = directive.argument;
			if (filename.find("posre") != std::string::npos || filename.find(".itp") == std::string::npos) {
				directive.ignoredInclude = true;	// Position restraints are not yet supported
			}
			else if (fs::exists(dir / filename))
				directive.includePath = dir / filename;
			else if (fs::exists(FileUtils::GetLimaDir() / "resources/forcefields" / filename))
				directive.includePath = FileUtils::GetLimaDir() / "resources/forcefields" / filename;
			// Else the file is missing, which is only an error if the include is active
		}
		return directive;
	}

	static void Scan(SourceFile& file) {
		file.contents = FileUtils::ReadFileToString(file.path);
		const std::string_view text = file.contents;
		const fs::path dir = file.path.parent_path();

		for (size_t lineBegin = 0; lineBegin < text.size();) {
			const size_t newline = text.find('\n', lineBegin);
			const size_t lineEnd = newline == std::string_view::npos ? text.size() : newline + 1;
			size_t first = lineBegin;
			while (first < lineEnd && (text[first] == ' ' || text[first] == '\t'))
				first++;

			if (first < lineEnd && (text[first] == '#' || text[first] == '[')) {
				std::string_view line = text.substr(first, (newline == std::string_view::npos ? text.size() : newline) - first);
				if (!line.empty() && line.back() == '\r')
					line.remove_suffix(1);

				if (line[0] == '#') {
					Directive directive = ParseDirective(line, dir);
					directive.lineBegin = lineBegin;
					directive.lineEnd = lineEnd;
					file.directives.emplace_back(std::move(directive));
				}
				else {
					file.headers.emplace_back(Header{ lineBegin, lineEnd, extractSectionName(std::string{ line }) });
				}
			}
			lineBegin = lineEnd;
		}
	}

	// ------------------------------------------ Phase 2: Resolve ------------------------------------------ //

	// Returns the index of the new instance
	int Resolve(int fileIndex, std::optional<fs::path> includeName, int depth) {
		if (depth > 64)
			throw std::runtime_error(std::format("Include depth exceeded 64 in file {}, is there a recursive #include?", files[fileIndex]->path.string()));

		const int instanceIndex = static_cast<int>(instances.size());
		instances.emplace_back(Instance{ fileIndex, includeName, {} });

		const SourceFile& file = *files[fileIndex];
		std::vector<Piece> pieces;	// Not directly in the instance, as recursion may reallocate instances

		struct Conditional { bool parentActive; bool condition; bool inElse; };
		std::vector<Conditional> conditionals;
		bool active = true;
		size_t cursor = 0;	// Start of the text since the previous directive

		auto Error = [&](std::string_view message) {
			return std::runtime_error(std::format("{} in file {}", message, file.path.string()));
			};

		using Kind = Directive::Kind;
		for (const Directive& directive : file.directives) {
			if (active && directive.lineBegin > cursor)
				pieces.emplace_back(Piece{ cursor, directive.lineBegin });
			cursor = directive.lineEnd;

			switch (directive.kind) {
			case Kind::Ifdef:
			case Kind::Ifndef: {
				const bool defined = topology.defines.contains(directive.argument);
				conditionals.emplace_back(Conditional{ active, directive.kind == Kind::Ifdef ? defined : !defined, false });
				active = active && conditionals.back().condition;
				break;
			}
			case Kind::If:
				if (active)
					throw Error("#if is not supported, only #ifdef and #ifndef");
				conditionals.emplace_back(Conditional{ false, false, false });
				break;
			case Kind::Elif:
				if (conditionals.empty())
					throw Error("#elif without #if");
				if (conditionals.back().parentActive)
					throw Error("#elif is not supported");
				break;
			case Kind::Else:
				if (conditionals.empty() || conditionals.back().inElse)
					throw Error("Unexpected #else");
				conditionals.back().inElse = true;
				active = conditionals.back().parentActive && !conditionals.back().condition;
				break;
			case Kind::Endif:
				if (conditionals.empty())
					throw Error("#endif without #ifdef");
				active = conditionals.back().parentActive;
				conditionals.pop_back();
				break;
			case Kind::Define:
				if (active)
					topology.defines.insert(directive.argument);
				break;
			case Kind::Undef:
				if (active)
					topology.defines.erase(directive.argument);
				break;
			case Kind::Include:
				if (!active || directive.ignoredInclude)
					break;
				if (directive.includeFile == -1)
					throw std::runtime_error(std::format("Could not find file \"{}\" in directory \"{}\"", directive.argument, file.path.parent_path().string()));
				pieces.emplace_back(Piece{ 0, 0, Resolve(directive.includeFile, directive.argument, depth + 1) });
				break;
			case Kind::Other:
				break;
			}
		}
		if (!conditionals.empty())
			throw Error("Missing #endif");
		if (active && cursor < file.contents.size())
			pieces.emplace_back(Piece{ cursor, file.contents.size() });

		instances[instanceIndex].pieces = std::move(pieces);
		return instanceIndex;
	}

	// ------------------------------------------ Phase 3: Parse ------------------------------------------ //

	void Walk(int instanceIndex) {
		const Instance& instance = instances[instanceIndex];
		const SourceFile& file = *files[instance.file];
		SectionState state;

		for (const Piece& piece : instance.pieces) {
			if (piece.childInstance != -1) {
				Walk(piece.childInstance);
				continue;
			}

			size_t pos = piece.begin;
			auto header = std::lower_bound(file.headers.begin(), file.headers.end(), piece.begin,
				[](const Header& h, size_t offset) { return h.lineBegin < offset; });
			for (; header != file.headers.end() && header->lineBegin < piece.end; ++header) {
				ProcessSectionText(state, instance, std::string_view{ file.contents }.substr(pos, header->lineBegin - pos));
				EnterSection(state, instance, header->name);
				pos = header->lineEnd;
			}
			ProcessSectionText(state, instance, std::string_view{ file.contents }.substr(pos, piece.end - pos));
		}
	}

	void EnterSection(SectionState& state, const Instance& instance, const std::string& name) {
		state.section = state.sectionGetter(name);

		if (state.section == defaults) {
			// This file is a forcefield
			if (topology.forcefieldInclude != std::nullopt)
				throw std::runtime_error("Trying to include a forcefield, but topology already has 1!");
			topology.forcefieldInclude.emplace(TopologyFile::ForcefieldInclude(fs::path{ instance.includeName.value_or("forcefield.itp") }));
		}
	}

	void ProcessSectionText(SectionState& state, const Instance& instance, std::string_view text) {
		if (text.empty())
			return;

		if (const auto it = std::ranges::find(entrySections, state.section); it != entrySections.end()) {
			// The expensive sections, we parse these later in parallel
			if (state.work != -1)
				works[state.work].chunks[it - entrySections.begin()].push_back(text);
			return;
		}

		switch (state.section) {
		case TopologySection::title:
			ForEachLine(text, [&](std::string_view line) {
				if (state.section != TopologySection::title)
					return;
				if (line.empty()) {
					state.section = no_section;
					return;
				}
				if (!instance.includeName.has_value() && !isOnlySpacesAndTabs(line))	// Only use main top title
					topology.title.append(std::string{ line } + "\n");
				});
			return;
		case TopologySection::moleculetype:
		case TopologySection::_system:
		case TopologySection::molecules:
		case TopologySection::defaults:
		case TopologySection::atomtypes:
		case TopologySection::pairtypes:
		case TopologySection::bondtypes:
		case TopologySection::constainttypes:
		case TopologySection::angletypes:
		case TopologySection::dihedraltypes:
		case TopologySection::impropertypes:
			ForEachLine(text, [&](std::string_view line) {
				if (IsDataLine(line))
					ProcessLine(state, instance, line);
				});
			return;
		default:
			return;	// Sections we dont use
		}
	}

	void ProcessLine(SectionState& state, const Instance& instance, std::string_view line) {
		const fs::path& path = files[instance.file]->path;

		switch (state.section) {
		case TopologySection::moleculetype: {
			std::string_view sv = line;
			SkipLeadingWhitespace(sv);
			std::string moleculetypename;
			int nrexcl = 0;
			ParseValue(sv, moleculetypename);
			ParseValue<int>(sv, nrexcl);

			if (moleculetypename.empty())
				throw std::runtime_error("Moleculetype name is empty in file: " + path.string());
			if (topology.moleculetypes.contains(moleculetypename)) {
				// The first definition wins, and the entries of this one are discarded.
				// GROMACS would reject this, but our own ToGmx outputs eg. SOL both from SOL.itp and tip3p.itp
				state.work = -1;
				break;
			}

			auto moleculetype = std::make_shared<TopologyFile::Moleculetype>(moleculetypename, nrexcl, instance.includeName);
			topology.moleculetypes.insert({ moleculetypename, moleculetype });
			works.emplace_back(MoleculetypeWork{ moleculetype, instance.file });
			files[instance.file]->unparsedMoleculetypes++;
			state.work = static_cast<int>(works.size() - 1);
			break;
		}
		case TopologySection::_system:
			topology.SetSystem(std::string{ line });
			break;
		case TopologySection::molecules: {
			std::string_view sv = line;
			SkipLeadingWhitespace(sv);
			std::string molname;
			int count = 0;
			ParseValue(sv, molname);
			ParseValue<int>(sv, count);

			if (!topology.HasSystem())
				throw std::runtime_error("Molecule section encountered before system section in file: " + path.string());
			const auto moleculetype = topology.moleculetypes.find(molname);
			if (moleculetype == topology.moleculetypes.end())
				throw std::runtime_error(std::format("Moleculetype {} not defined before being used in file: {}", molname, path.string()));
			if (count > 0)
				topology.m_system.molecules.emplace_back(TopologyFile::MoleculeEntry{ molname, moleculetype->second, count });
			break;
		}
		default:	// Forcefield sections
			if (!topology.forcefieldInclude)
				throw std::runtime_error(std::format("Forcefield section encountered before [ defaults ] in file: {}", path.string()));
			topology.forcefieldInclude->AddEntry(state.section, std::string{ line });
			break;
		}
	}

	static std::vector<std::string_view> DataLines(const std::vector<std::string_view>& chunks) {
		std::vector<std::string_view> lines;
		for (const auto chunk : chunks)
			ForEachLine(chunk, [&](std::string_view line) {
				if (IsDataLine(line))
					lines.push_back(line);
				});
		return lines;
	}

	static void ParseMoleculetype(MoleculetypeWork& work) {
		TopologyFile::Moleculetype& moleculetype = *work.moleculetype;
		auto Lines = [&](TopologySection section) {
			return DataLines(work.chunks[std::ranges::find(entrySections, section) - entrySections.begin()]);
			};

		const auto atomLines = Lines(atoms);
		if (atomLines.size() >= 999'999)
			throw std::runtime_error("file contained more that 999'999 atoms. This makes their id non-unique, due to limitations of the format. Please split your file into multiple topologies.");

		moleculetype.atoms.resize(atomLines.size());
		std::vector<int> limaIdToGroId(atomLines.size());
		ForEachEntry(atomLines.size(), [&](size_t i) {
			limaIdToGroId[i] = ParseAtomsEntry(atomLines[i], moleculetype.atoms[i], static_cast<int>(i));
			});

		// Invert the mapping. If a groId appears multiple times, the last atom wins
		const int maxGroId = limaIdToGroId.empty() ? 0 : std::ranges::max(limaIdToGroId);
		std::vector<int> groIdToLimaId(maxGroId + 1, -1);
		for (int i = 0; i < static_cast<int>(limaIdToGroId.size()); i++)
			groIdToLimaId[limaIdToGroId[i]] = i;

		ParseBonds(Lines(bonds), moleculetype.singlebonds, groIdToLimaId);
		ParseBonds(Lines(pairs), moleculetype.pairbonds, groIdToLimaId);
		ParseBonds(Lines(angles), moleculetype.anglebonds, groIdToLimaId);
		ParseBonds(Lines(dihedrals), moleculetype.dihedralbonds, groIdToLimaId);
		ParseBonds(Lines(impropers), moleculetype.improperdihedralbonds, groIdToLimaId);
		ParseBonds(Lines(cmap), moleculetype.cmapbonds, groIdToLimaId);
	}
};


TopologyFile::TopologyFile() = default;
TopologyFile::TopologyFile(const fs::path& path) : path(path)
{
	if (!(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp"))
		throw std::runtime_error("Expected .top or .itp extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	TopologyParser{ *this }.Parse(path);
}



GenericItpFile::GenericItpFile(const fs::path& path) {
	if (!(path.extension().string() == ".itp" || path.extension().string() == ".top")) { throw std::runtime_error(std::format("Expected .itp extension with file {}", path.string())); }
	if (!fs::exists(path)) { throw std::runtime_error(std::format("File \"{}\" was not found", path.string())); }

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	TopologySection current_section{ TopologySection::title };
	TopologySectionGetter getTopolSection{};
	bool newSection = true;

	std::string line{}, word{};
	std::string sectionname = "";
	std::vector<int> groIdToLimaId;


	while (getline(file, line)) {

		if (line.empty())
			continue;

		if (line.find("#include") != std::string::npos) {
			if (!sections.contains(includes))
				sections.insert({ includes, {} });

			sections.at(includes).emplace_back(FileUtils::ExtractBetweenQuotemarks(line));
			continue;
		}

		if (HandleTopologySectionStartAndStop(std::string_view( line ), current_section, getTopolSection)) {
			newSection = true;
			continue;
		}

		// Check if current line is commented
		if (firstNonspaceCharIs(line, TopologyFile::commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) { continue; }	// Only title-sections + atoms reads the comments

		// Check if this line contains another illegal keyword
		if (firstNonspaceCharIs(line, '#')) {// Skip ifdef, include etc TODO: this needs to be implemented at some point
			continue;
		}


		if (newSection) {
			if (!GetSection(current_section).empty()) {
				throw std::runtime_error("Found the same section muliple times in the same file");
			}
			sections.insert({ current_section, {} });
			newSection = false;
		}

		GetSection(current_section).emplace_back(line);	// OPTIM: i prolly should cache this address instead
	}
}

void GenericItpFile::printToFile(const fs::path& path) const {
	if (path.extension() != ".itp") {
		throw std::runtime_error(std::format("Expected .itp extension with file {}", path.string()));
	}
	std::ofstream file(path);
	if (!file) throw std::runtime_error(std::format("Failed to create {}", path.string()));

	for (const auto& include : GetSection(includes)) file << "#include \"" << include << "\"\n";
	constexpr std::array sectionNames{
		std::pair{ defaults, "defaults" }, std::pair{ atomtypes, "atomtypes" },
		std::pair{ pairtypes, "pairtypes" }, std::pair{ bondtypes, "bondtypes" },
		std::pair{ constainttypes, "constrainttypes" }, std::pair{ angletypes, "angletypes" },
		std::pair{ dihedraltypes, "dihedraltypes" }, std::pair{ impropertypes, "dihedraltypes" },
		std::pair{ cmaptypes, "cmaptypes" }, std::pair{ moleculetype, "moleculetype" },
		std::pair{ atoms, "atoms" }, std::pair{ bonds, "bonds" }, std::pair{ pairs, "pairs" },
		std::pair{ angles, "angles" }, std::pair{ dihedrals, "dihedrals" },
		std::pair{ impropers, "dihedrals" }, std::pair{ cmap, "cmap" },
		std::pair{ position_restraints, "position_restraints" },
		std::pair{ _system, "system" }, std::pair{ molecules, "molecules" }
	};
	for (const auto& [section, name] : sectionNames) {
		const auto& entries = GetSection(section);
		if (entries.empty()) continue;
		file << "\n[ " << name << " ]\n";
		for (const auto& entry : entries) file << entry << '\n';
	}
}



void TopologyFile::ForcefieldInclude::AddEntry(TopologySection section, const std::string& entry) {
	contents.GetSection(section).emplace_back(entry);
}


void TopologyFile::ForcefieldInclude::SaveToDir(const fs::path& directory) const {
	if (!fs::is_directory(directory)) {
		throw std::runtime_error(std::format("Directory \"{}\" does not exist", directory.string()));
	}
	if (filename.extension() != ".itp")
		throw std::runtime_error("Forcefield include name must have .itp extension");	
	std::ofstream file(directory / "forcefield.itp");
	if (!file.is_open()) {
		throw std::runtime_error(std::format("Failed to open file {}", (directory / "forcefield.itp").string()));
	}

	file << "[ defaults ]\n";
	for (const auto& entry : contents.GetSection(defaults)) {
		file << entry << '\n';
	}

	file << "[ atomtypes ]\n";
	for (const auto& entry : contents.GetSection(atomtypes)) {
		file << entry << '\n';
	}

	file << "[ pairtypes ]\n";
	for (const auto& entry : contents.GetSection(pairtypes)) {
		file << entry << '\n';
	}

	file << "[ bondtypes ]\n";
	for (const auto& entry : contents.GetSection(bondtypes)) {
		file << entry << '\n';
	}

	file << "[ angletypes ]\n";
	for (const auto& entry : contents.GetSection(angletypes)) {
		file << entry << '\n';
	}

	file << "[ dihedraltypes ]\n";
	for (const auto& entry : contents.GetSection(dihedraltypes)) {
		file << entry << '\n';
	}

	file << "[ dihedraltypes ]\n";
	for (const auto& entry : contents.GetSection(impropertypes)) {
		file << entry << '\n';
	}
}



void TopologyFile::AppendMolecule(const std::string& moleculename) {
	if (!m_system.IsInit()) {
		throw std::runtime_error("System is not initialized");
	}
	if (!moleculetypes.contains(moleculename)) {
		throw std::runtime_error(std::format("Moleculetype {} not found in topology", moleculename));
	}

	if (!m_system.molecules.empty() && m_system.molecules.back().name == moleculename)
		m_system.molecules.back().count++;
	else
		m_system.molecules.emplace_back(MoleculeEntry{ moleculename, moleculetypes.at(moleculename) });
}
void TopologyFile::AppendMoleculetype(const std::shared_ptr<const Moleculetype> moleculetype, std::optional<ForcefieldInclude> inputForcefieldInclude) {
	if (!moleculetypes.contains(moleculetype->name)) {
		moleculetypes.insert({ moleculetype->name, std::make_shared<Moleculetype>(*moleculetype) });	//COPY

		if (inputForcefieldInclude.has_value()) {
			if (!forcefieldInclude.has_value())
				forcefieldInclude.emplace(inputForcefieldInclude.value());
			/*else 
				assert(forcefieldInclude->filename == inputForcefieldInclude->filename);			*/ // TODO: Handle this somehow?? For now its not an issue since it occurs when i've given appropriate forcefields.. 
		}
	}
	AppendMolecule(moleculetype->name);
}
//
//void TopologyFile::AppendSolvents(int count, const fs::path& solventFF) {
//
//	ParseFileIntoTopology(*this, solventFF);
//
//
//	//TopologyFile solventTopology(solventFF);
//
//	for (int i = 0; i < count; i++)
//		AppendMolecule("SOL");
//		//AppendMoleculetype(solventTopology.GetMoleculeTypePtr(), std::nullopt);
//
//	/*for (const auto& molecule : m_system.molecules) {
//		
//	}*/
//}

void TopologyFile::printToFile(const std::filesystem::path& path) const {
	const auto ext = path.extension().string();
	if (ext != ".top" && ext != ".itp") { throw std::runtime_error(std::format("Got {} extension, expected [.top/.itp]", ext)); }
	if (!path.parent_path().empty())
		fs::create_directories(path.parent_path());
	{
		std::ofstream file(path);
		if (!file.is_open()) {
			throw std::runtime_error(std::format("Failed to open file {}", path.string()));
		}
		
		file << "; " << title << "\n\n";

		// TODO: Have multiple forcefields, just only 1 with the [ defaults ] directive
		{				
			if (forcefieldInclude) {
				bool usesInternalForcefield = fs::exists(GetLimaDir() / "resources/forcefields" / forcefieldInclude->filename);
				if (usesInternalForcefield) {
					file << ("#include \"" + forcefieldInclude->filename.string() + "\"\n");
				}
				else {
					forcefieldInclude.value().SaveToDir(path.parent_path());
					file << ("#include \"forcefield.itp\"\n");
				}
				file << "\n";
			}			
		}
		for (const auto& [_, moleculetype] : moleculetypes) {
			moleculetype->ToFile(path.parent_path());
			file << "#include \"" << moleculetype->includePath.value_or(fs::path(moleculetype->name + ".itp")).string() << "\"\n";
		}
		for (const auto& include : otherIncludes) file << "#include \"" << include << "\"\n";
		file << "\n";

		if (m_system.IsInit()) {
			file << "[ system ]\n";
			file << m_system.title << "\n\n";

			file << "[ molecules ]\n";
			// If we have the same molecule multiple times in a row, we only print it once together with the total count
			for (size_t i = 0; i < m_system.molecules.size(); i++) {
				int count = m_system.molecules[i].count;
				while (i + 1 < m_system.molecules.size() && m_system.molecules[i].name == m_system.molecules[i + 1].name) {
					count += m_system.molecules[i + 1].count;
					i++;
				}
				file << m_system.molecules[i].name << " " << count << "\n";
			}
		}
		file << "\n";

		//if (!molecules.entries.empty()) { file << molecules.composeString(); }
		// If we have the same submolecule multiple times in a row, we only print it once together with a count of how many there are

		//if (!molecules.empty())
		//	file << molecules.title << "\n" << molecules.legend << "\n";
		//for (int i = 0; i < molecules.entries.size(); i++) {
		//	std::ostringstream oss;
		//	int count = 1;
		//	while (i + 1 < molecules.entries.size() && molecules.entries[i].name == molecules.entries[i + 1].name) {
		//		count++;
		//		i++;
		//	}
		//	file << molecules.entries[i].includeTopologyFile->name << " " << count << "\n";
		//}
	}

	// Also cache the file
	//WriteFileToBinaryCache(*this, path);

	/*for (const auto& [name, include] : includeTopologies) {
		include->printToFile(path.parent_path() / ("topol_" + name + ".itp"), printForcefieldinclude);
	}*/
}






std::string generateLegend(const std::vector<std::string>& elements)
{
	std::ostringstream legend;
	legend << ';'; // Start with a semicolon

	for (const auto& element : elements) {
		legend << std::setw(10) << std::right << element;
	}
	return legend.str();
}



template <typename T>
std::string composeString(const std::vector<T>&elements) {
	std::ostringstream oss;
	for (const auto& entry : elements) {
		entry.composeString(oss);
	}
	oss << '\n';
	return oss.str();
}

void TopologyFile::Moleculetype::ToFile(const fs::path& dir) const {

	const fs::path path = dir / includePath.value_or(fs::path{ name + ".itp" });
	
	if (atoms.empty())
		throw(std::runtime_error("Trying to print moleculetype to file, but it has no atoms"));

	{
		std::ofstream file(path);
		if (!file.is_open()) {
			throw std::runtime_error(std::format("Failed to open file {}", path.string()));
		}

		file << "; " << name << "\n\n";

		file << "[ moleculetype ]\n";
		file << generateLegend({ "name", "nrexcl" }) + "\n";
		file << std::right << std::setw(10) << name << std::setw(10) << nrexcl << "\n\n";
		
		
		
		
		file << "[ atoms ]\n" << generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) + "\n";
		file << composeString(atoms);
		file << "[ bonds ]\n" << generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) + "\n";
		file << composeString(singlebonds);
		file << "[ pairs ]\n" << generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) + "\n";
		file << composeString(pairbonds);
		file << "[ angles ]\n" << generateLegend({ "ai", "aj", "ak", "funct", "c0", "c1", "c2", "c3" }) + "\n";
		file << composeString(anglebonds);
		file << "[ dihedrals ]\n" << generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3", "c4", "c5" }) + "\n";
		file << composeString(dihedralbonds);
		file << "[ dihedrals ]\n" << generateLegend({ "ai", "aj", "ak", "al", "funct", "c0", "c1", "c2", "c3" }) + "\n";
		file << composeString(improperdihedralbonds);
		file << "[ cmap ]\n" << generateLegend({ "ai", "aj", "ak", "al", "am", "funct" }) + "\n";
		file << composeString(cmapbonds);
		if (positionRestraintsInclude) {
			file << "#ifdef POSRES\n#include \"" << positionRestraintsInclude->string() << "\"\n#endif\n";
		}
	}
}
void TopologyFile::AtomsEntry::composeString(std::ostringstream& oss) const {
	if (section_name) {
		oss << section_name.value() << "\n";
	}
	oss << std::right
		<< std::setw(10) << id + 1 // convert back to 1-indexed
		<< std::setw(10) << type
		<< std::setw(10) << resnr
		<< std::setw(10) << residue
		<< std::setw(10) << atomname
		<< std::setw(10) << cgnr;
	if (charge.has_value())
		oss << std::setw(10) << std::fixed << std::setprecision(2) << charge.value();
	if (mass.has_value())
		oss << std::setw(10) << std::fixed << std::setprecision(3) << mass.value();
	oss << '\n';
}

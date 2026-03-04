#include "MDFiles.h" 
#include "Filehandling.h"

#
#include <algorithm>
#include <format>
#include <ranges>
#include "TimeIt.h"
#include <execution>













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

constexpr std::string_view extractSectionName(const std::string& line) {
	size_t start = line.find('[');
	size_t end = line.find(']', start);
	if (start != std::string::npos && end != std::string::npos) {
		// Extract the view of the text between '[' and ']'
		std::string_view sectionView(line.c_str() + start + 1, end - start - 1);

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


constexpr bool VerifyAllParticlesInBondExists(const std::unordered_map<int, int>& groIdToLimaId, std::span<const int> ids) {
	for (const auto id : ids) {
		if (!groIdToLimaId.contains(id) || groIdToLimaId.at(id) == -1)
			return false;
	}
	return true;
}

/// <returns>True if we change section/stay on no section, in which case we should 'continue' to the next line. False otherwise</returns>
constexpr bool HandleTopologySectionStartAndStop(const std::string& line, TopologySection& currentSection, TopologySectionGetter& sectionGetter) {

	if (line.empty() && currentSection == TopologySection::title) {
		currentSection = no_section;
		return true;
	}
	else if (!line.empty() && line[0] == '[') {
		currentSection = sectionGetter(extractSectionName(line));
		return true;
	}

	return false; // We did not change section, and should continue parsing the current line
}
constexpr bool isOnlySpacesAndTabs(const std::string& str) {
	return std::all_of(str.begin(), str.end(), [](char c) {
		return c == ' ' || c == '\t' || c == '\r';
		});
}

inline void SkipLeadingWhitespace(std::string_view& sv) noexcept {
	// ASCII whitespace incl. CR/LF/TAB/VT/FF
	const auto p = sv.find_first_not_of(" \t\r\n\v\f");
	sv.remove_prefix(p == std::string_view::npos ? sv.size() : p);
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
	// skip spaces
	auto pos = sv.find_first_not_of(" \t");
	sv.remove_prefix(pos == std::string_view::npos ? sv.size() : pos);
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

	// skip trailing spaces after the token
	const auto skip = sv.find_first_not_of(" \t");
	sv.remove_prefix(skip == std::string_view::npos ? sv.size() : skip);
	return !out.empty();
}

void TopologyFile::ParseAtomsEntry(std::string_view sv, TopologyFile::AtomsEntry& atom, std::vector<int>& limaIdToGroId, int index /*relative to moleculetype*/) {
	SkipLeadingWhitespace(sv);

	if (firstNonspaceCharIs(sv, ';')) {
		// Skip the very first line which is the legend
		/*if (sv.find("cgnr") != std::string_view::npos) {
			return;
		}*/
		/*if (sv.find("residue") != std::string::npos || sv.find("lipid_section") != std::string::npos)
			moleculetype.mostRecentAtomsSectionName = sv;*/
	}
	else {
		int groId;
		ParseValue<int>(sv, groId);
		ParseValue(sv, atom.type);
		ParseValue<int>(sv, atom.resnr);
		ParseValue(sv, atom.residue);
		ParseValue(sv, atom.atomname);
		ParseValue(sv, atom.cgnr);
		ParseValue(sv, atom.charge);
		ParseValue(sv, atom.mass);	// This might not be present

		/*if (groIdToLimaId.size() < groId + 1)
			throw std::runtime_error(std::format("Atom with groId {} found, but only {} atoms have been defined so far. This is most likely due to the numbering not starting at 1, or not being sequential", groId, groIdToLimaId.size()));*/

		limaIdToGroId[index] = groId;
			//groIdToLimaId[groId] = index;
		atom.id = index; //groIdToLimaId[groId];	

		if (atom.type.empty() || atom.residue.empty() || atom.atomname.empty())
			throw std::runtime_error("Atom type, residue or atomname is empty");
	}
}

template <int n>
bool LoadIds(std::string_view& sv, std::array<int, n>& ids, const std::unordered_map<int, int>& groIdToLimaId, bool&err) {
	for (int i = 0; i < n; i++) {
		if (!ParseValue<int>(sv, ids[i])) {
			err = true;
			return false;
		}
	}

	if (!VerifyAllParticlesInBondExists(groIdToLimaId, ids)) {
		err = true;
		return false;
	}

	for (int i = 0; i < n; i++)
		ids[i] = groIdToLimaId.at(ids[i]);

	return true;
}
void TopologyFile::ParseSingleBond(std::string_view sv, TopologyFile::SingleBond& bond, const std::unordered_map<int, int>& groIdToLimaId, bool& error) {
	SkipLeadingWhitespace(sv);
    if (!LoadIds<2>(sv, bond.ids, groIdToLimaId, error))
		return;

	float b0, kb;	
	bool err = false;
	err |= !ParseValue<int>(sv, bond.funct);
	err |= !ParseValue<float>(sv, b0);
	err |= !ParseValue<float>(sv, kb);
	if (!err)	// Some top files have this data, some dont. If it exists, it takes precedence over forcefield
		bond.parameters = Bondtypes::SingleBond::Parameters::CreateFromCharmm(b0, kb);
	//singlebond.sourceLine = line;
}

void TopologyFile::ParsePairBond(std::string_view sv, TopologyFile::PairBond& bond, const std::unordered_map<int, int>& groIdToLimaId, bool& error) {
	SkipLeadingWhitespace(sv);
    if (!LoadIds<2>(sv, bond.ids, groIdToLimaId, error))
		return;

	float sigma, epsilon;
	bool err = false;
	err |= !ParseValue<int>(sv, bond.funct);
	err |= !ParseValue<float>(sv, sigma);
	err |= !ParseValue<float>(sv, epsilon);
	if (!err)
		bond.parameters = Bondtypes::PairBond::Parameters::CreateFromCharmm(sigma, epsilon);
	//pairbond.sourceLine = line;
}

void TopologyFile::ParseAngleBond(std::string_view sv, TopologyFile::AngleBond& bond, const std::unordered_map<int, int>& groIdToLimaId, bool& error) {
	SkipLeadingWhitespace(sv);
    if (!LoadIds<3>(sv, bond.ids, groIdToLimaId, error))
		return;
	

	float theta0, ktheta, ub0, kUb;
	bool err = false;
	err |= !ParseValue<int>(sv, bond.funct);
	err |= !ParseValue<float>(sv, theta0);
	err |= !ParseValue<float>(sv, ktheta);
	err |= !ParseValue<float>(sv, ub0);
	err |= !ParseValue<float>(sv, kUb);
	if (!err)
		bond.parameters = Bondtypes::AngleUreyBradleyBond::Parameters::CreateFromCharmm(theta0, ktheta, ub0, kUb, bond.funct);
}

void TopologyFile::ParseDihedralBond(std::string_view sv, TopologyFile::DihedralBond& bond, const std::unordered_map<int, int>& groIdToLimaId, bool& error) {
	SkipLeadingWhitespace(sv);
    if (!LoadIds<4>(sv, bond.ids, groIdToLimaId, error))
		return;

	bool err = false;
	float phi0, kphi; int n;
	err |= !ParseValue<int>(sv, bond.funct);
	err |= !ParseValue<float>(sv, phi0);
	err |= !ParseValue<float>(sv, kphi);
	err |= !ParseValue<int>(sv, n);
	if (!err)
		bond.parameters = Bondtypes::DihedralBond::Parameters::CreateFromCharmm(phi0, kphi, n);
}

void TopologyFile::ParseImproperDihedralBond(std::string_view sv, TopologyFile::ImproperDihedralBond& bond, const std::unordered_map<int, int>& groIdToLimaId, bool& error) {
	SkipLeadingWhitespace(sv);
    if (!LoadIds<4>(sv, bond.ids, groIdToLimaId, error))
		return;

	float phi0, kphi;
	bool err = false;
	err |= !ParseValue<int>(sv, bond.funct);
	err |= !ParseValue<float>(sv, phi0);
	err |= !ParseValue<float>(sv, kphi);
	if (!err)
		bond.parameters = Bondtypes::ImproperDihedralBond::Parameters::CreateFromCharmm(phi0, kphi);
}

void TopologyFile::ParseMoleculetypeEntry(TopologySection section, const std::string& line, std::shared_ptr<Moleculetype> moleculetype) {
	std::istringstream iss(line);

	switch (section)
	{
	case TopologySection::atoms:
	{
		//if (firstNonspaceCharIs(line, ';')) {
		//	// Skip the very first line which is the legend
		//	if (line.find("cgnr") != std::string::npos) {
		//		break;
		//	}
		//	if (line.find("residue") != std::string::npos || line.find("lipid_section") != std::string::npos)
		//		moleculetype->mostRecentAtomsSectionName = line;
		//}
		//else {
		//	TopologyFile::AtomsEntry atom;
		//	int groId;
		//	iss >> groId >> atom.type >> atom.resnr >> atom.residue >> atom.atomname >> atom.cgnr >> atom.charge >> atom.mass;

		//	if (moleculetype->groIdToLimaId.size() < groId + 1)
		//		moleculetype->groIdToLimaId.resize(groId + 1, -1);
		//	moleculetype->groIdToLimaId[groId] = moleculetype->atoms.size();
		//	atom.id = moleculetype->groIdToLimaId[groId];
		//	moleculetype->atoms.emplace_back(atom);

		//	if (moleculetype->mostRecentAtomsSectionName != "") {
		//		moleculetype->atoms.back().section_name = moleculetype->mostRecentAtomsSectionName;
		//		moleculetype->mostRecentAtomsSectionName = "";
		//	}
		//	if (atom.type.empty() || atom.residue.empty() || atom.atomname.empty())
		//		throw std::runtime_error("Atom type, residue or atomname is empty");

		//}
		//break;
	}
	case TopologySection::bonds: {
		TopologyFile::SingleBond singlebond{};
		int groIds[2];
		iss >> groIds[0] >> groIds[1] >> singlebond.funct;
		if (!VerifyAllParticlesInBondExists(moleculetype->groIdToLimaId, groIds))
			break;
		for (int i = 0; i < 2; i++)
			singlebond.ids[i] = moleculetype->groIdToLimaId[groIds[i]];
		//singlebond.sourceLine = line;

		float b0, kb;
		if (iss >> b0 >> kb) {																						// TODO LONG: it is a very bad idea that we interprete ff params both here and in forcefield.cpp, we should ONLY do that 1 plac
			singlebond.parameters = Bondtypes::SingleBond::Parameters::CreateFromCharmm(b0, kb);
		}

		moleculetype->singlebonds.emplace_back(singlebond);
		break;
	}
	case TopologySection::pairs: {
		TopologyFile::PairBond pairbond{};
		int groIds[2];
		iss >> groIds[0] >> groIds[1] >> pairbond.funct;
		if (!VerifyAllParticlesInBondExists(moleculetype->groIdToLimaId, groIds))
			break;
		for (int i = 0; i < 2; i++)
			pairbond.ids[i] = moleculetype->groIdToLimaId.at(groIds[i]);

		float sigma, epsilon;
		if (iss >> sigma >> epsilon) {
			pairbond.parameters = Bondtypes::PairBond::Parameters::CreateFromCharmm(sigma, epsilon);
		}

		moleculetype->pairbonds.emplace_back(pairbond);
		break;
	}
	case TopologySection::angles: {
		TopologyFile::AngleBond angle{};
		int groIds[3];
		iss >> groIds[0] >> groIds[1] >> groIds[2] >> angle.funct;
		if (!VerifyAllParticlesInBondExists(moleculetype->groIdToLimaId, groIds))
			break;
		for (int i = 0; i < 3; i++)
			angle.ids[i] = moleculetype->groIdToLimaId.at(groIds[i]);

		float theta0, ktheta, ub0, kUb;
		if (iss >> theta0 >> ktheta >> ub0 >> kUb) {
			angle.parameters = Bondtypes::AngleUreyBradleyBond::Parameters::CreateFromCharmm(theta0, ktheta, ub0, kUb, angle.funct);
		}

		moleculetype->anglebonds.emplace_back(angle);
		break;
	}
	case TopologySection::dihedrals: {
		TopologyFile::DihedralBond dihedral{};
		int groIds[4];
		iss >> groIds[0] >> groIds[1] >> groIds[2] >> groIds[3] >> dihedral.funct;
		if (!VerifyAllParticlesInBondExists(moleculetype->groIdToLimaId, groIds))
			break;
		for (int i = 0; i < 4; i++)
			dihedral.ids[i] = moleculetype->groIdToLimaId.at(groIds[i]);

		float phi0, kphi;
		int n;
		if (iss >> phi0 >> kphi >> n) {
			dihedral.parameters = Bondtypes::DihedralBond::Parameters::CreateFromCharmm(phi0, kphi, n);
		}

		moleculetype->dihedralbonds.emplace_back(dihedral);
		break;
	}
	case TopologySection::impropers: {
		TopologyFile::ImproperDihedralBond improper{};
		int groIds[4];
		iss >> groIds[0] >> groIds[1] >> groIds[2] >> groIds[3] >> improper.funct;
		if (!VerifyAllParticlesInBondExists(moleculetype->groIdToLimaId, groIds))
			break;
		for (int i = 0; i < 4; i++)
			improper.ids[i] = moleculetype->groIdToLimaId.at(groIds[i]);

		float phi0, kphi;
		if (iss >> phi0 >> kphi) {
			improper.parameters = Bondtypes::ImproperDihedralBond::Parameters::CreateFromCharmm(phi0, kphi);
		}

		moleculetype->improperdihedralbonds.emplace_back(improper);
		//moleculetype->improperdihedralbonds.back().sourceLine = line;
		break;
	}
	default: {
		// We shouldnt get here
	}
	}
}

void TopologyFile::ParseFileIntoTopology(TopologyFile& topology, const fs::path& path, std::optional<fs::path> includefileName) {
	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	topology.defines.insert("FLEXIBLE");// Cant handle gromacs definition of rigid water right now

	TopologySection current_section{ TopologySection::title };
	TopologySectionGetter getTopolSection{};
	std::shared_ptr<Moleculetype> mostRecentMoleculetype = nullptr;

	std::string line{};


	// Data for processing
	/*std::string superbuffer;
	std::vector<std::string_view> singlebondStrings;*/
	std::vector<::std::string> atomStrings;
	std::vector<std::string> singlebondStrings;
	std::vector<std::string> pairbondStrings;
	std::vector<std::string> anglebondStrings;
	std::vector<std::string> dihedralbondStrings;
	std::vector<std::string> improperbondStrings;

	while (getline(file, line)) {
		if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {

			// Directives where the directive itself is enough
			if (current_section == defaults) {
				// This file is a forcefield. We add it to the includes, and return to parent topol
				if (topology.forcefieldInclude != std::nullopt)
					throw std::runtime_error("Trying to include a forcefield, but topology already has 1!");

				topology.forcefieldInclude.emplace(ForcefieldInclude(fs::path{includefileName.value_or("forcefield.itp")}));
			}
			continue;
		}

		if (line.empty() || isOnlySpacesAndTabs(line))
			continue;

		// Check if current line is commented
		//if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) {	
		if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title) {// Currently skipping these lines from topologiues: ; residue   1 MET rtp MET  q +1.0 
			continue;
		}	// Only title-sections + atoms reads the comments

		if (line[0] == '#') {
			if (FileUtils::ChechlineForDefine(line)) {
				topology.defines.insert(FileUtils::ChechlineForDefine(line).value());
				continue;
			}

			if (FileUtils::ChecklineForIfdefAndSkipIfFound(file, line, topology.defines))
				continue;

			if (line.size() > 8 && line.substr(0, 8) == "#include") {
				// take second word, remove "
				std::istringstream iss(line);
				std::string _, pathWithQuotes;
				iss >> _ >> pathWithQuotes;
				if (pathWithQuotes.size() < 3)
					throw std::runtime_error("Include is not formatted as expected: " + line);

				std::string filename = pathWithQuotes.substr(1, pathWithQuotes.size() - 2);

				// TODO: Check that we havent' already parsed this file

				if (filename.find("posre") != std::string::npos) {
					// Do nothing, not yet supported
				}
				else if (filename.find(".itp") != std::string::npos) {
					const fs::path filepath(path.parent_path() / filename);
					if (fs::exists(filepath))
						ParseFileIntoTopology(topology, path.parent_path() / filename, filename);
					else if (fs::exists(FileUtils::GetLimaDir() / "resources/forcefields" / filename))
						ParseFileIntoTopology(topology, FileUtils::GetLimaDir() / "resources/forcefields" / filename, filename);
					else
						throw std::runtime_error(std::format("Could not find file \"{}\" in directory \"{}\"", filename, path.parent_path().string()));
				}
			}
			continue;
		}

		//TimeIt time("entry");
		// Directives where w eread the contents
		switch (current_section)
		{
		case TopologySection::title:
			if (!includefileName.has_value())	// Only use main top title
				topology.title.append(line + "\n");	// +\n because getline implicitly strips it away.
			break;
		case TopologySection::moleculetype:
		{
			std::istringstream iss(line);
			std::string moleculetypename;
			int nrexcl;
			iss >> moleculetypename >> nrexcl;

			if (moleculetypename.empty())
				throw std::runtime_error("Moleculetype name is empty in file: " + path.string());

			mostRecentMoleculetype = std::make_shared<Moleculetype>(moleculetypename, nrexcl, includefileName);
			//mostRecentMoleculetype->name = moleculetypename;
			//mostRecentMoleculetype->nrexcl = nrexcl;

			//auto nextSection = ParseMoleculetype(file, moleculetype);
			assert(!topology.moleculetypes.contains(moleculetypename));
			topology.moleculetypes.insert({ moleculetypename, mostRecentMoleculetype });

			//current_section = nextSection;
			break;
		}		
		case TopologySection::_system: {
			topology.SetSystem(line);
			break;
		}
		case TopologySection::molecules: {
			std::istringstream iss(line);

			std::string molname;
			int cnt = 0;
			iss >> molname >> cnt;

			if (topology.m_system.title == "noSystem")
				throw std::runtime_error("Molecule section encountered before system section in file: " + path.string());
			if (!topology.moleculetypes.contains(molname))
				throw std::runtime_error(std::format("Moleculetype {} not defined before being used in file: {}", molname, path.string()));
			for (int i = 0; i < cnt; i++)
				topology.m_system.molecules.emplace_back(MoleculeEntry{ molname, topology.moleculetypes.at(molname) });
			break;
		}
		case TopologySection::atoms:
			atomStrings.push_back(std::move(line));
			break;
		/*{
			TimeIt time("entry");
			if (mostRecentMoleculetype == nullptr)
				throw std::invalid_argument("Moleculetype not set before parsing atoms/bonds/pairs/angles/dihedrals/impropers");
			ParseMoleculetypeEntry(current_section, line, mostRecentMoleculetype);
			time.stop();		
			break;
		}*/
		case TopologySection::bonds:
			singlebondStrings.push_back(std::move(line));
			break;
		case TopologySection::pairs:
			pairbondStrings.push_back(std::move(line));
			break;
		case TopologySection::angles:
			anglebondStrings.push_back(std::move(line));
			break;
		case TopologySection::dihedrals:
			dihedralbondStrings.push_back(std::move(line));
			break;
		case TopologySection::impropers:
			improperbondStrings.push_back(std::move(line));
			break;
		case TopologySection::defaults:
		case TopologySection::atomtypes:
		case TopologySection::pairtypes:
		case TopologySection::bondtypes:
		case TopologySection::constainttypes:
		case TopologySection::angletypes:
		case TopologySection::dihedraltypes:
		case TopologySection::impropertypes:
			topology.forcefieldInclude->AddEntry(current_section, line);
			break;
		default:
			// Do nothing
			//throw std::runtime_error("Illegal state");
			break;
		}
		// This switch steals the line, DO NOT ADD CODE AFTER HERE
		//time.stop();
	}
	

	if (!mostRecentMoleculetype)
		return;





	//TimeIt time("Parsing");

	bool error = false;


	if (atomStrings.size() >= 999'999)
		throw std::runtime_error("file contained more that 999'999 atoms. This makes their id non-unique, due to limitations of the format. Please split your file into multiple topologies.");



	mostRecentMoleculetype->atoms.resize(atomStrings.size());
	mostRecentMoleculetype->groIdToLimaId.reserve(atomStrings.size());
	//mostRecentMoleculetype->groIdToLimaId.resize(atomStrings.size() + 1, -1); // as groid is 1-indexed
	std::vector<int>limaIdToGroId(atomStrings.size());
	{
		auto indices = std::views::iota(size_t{ 0 }, atomStrings.size());
		std::for_each(
			std::execution::par_unseq,
			indices.begin(),
			indices.end(),
			[&](int i) {
				ParseAtomsEntry(atomStrings[i], mostRecentMoleculetype->atoms[i], limaIdToGroId, i);
			});
	}

	// Now invert the mapping
	for (int i = 0; i < limaIdToGroId.size(); i++) {
		int groId = limaIdToGroId[i];
		mostRecentMoleculetype->groIdToLimaId[groId] = i;
	}


	//for (int i = 0; i < atomStrings.size(); i++) {	// This loop sadly must be sequential for now, due to the way we read residue names..
	//	ParseAtomsEntry(atomStrings[i], mostRecentMoleculetype->atoms[i], mostRecentMoleculetype->groIdToLimaId, i);
	//}
	// groIdToLimaId valid after this loop


	auto ParseBonds = [mostRecentMoleculetype](const auto& strings, auto& bonds, auto Parser, bool& err) {
		auto indices = std::views::iota(size_t{ 0 }, strings.size());
		std::for_each(
			std::execution::par_unseq,
			indices.begin(),
			indices.end(),
			[&](int i) {
				Parser(
					strings[i],
					bonds[i],
					mostRecentMoleculetype->groIdToLimaId,
					err
				);
			});
		};



	mostRecentMoleculetype->singlebonds.resize(singlebondStrings.size());
	mostRecentMoleculetype->pairbonds.resize(pairbondStrings.size());
	mostRecentMoleculetype->anglebonds.resize(anglebondStrings.size());
	mostRecentMoleculetype->dihedralbonds.resize(dihedralbondStrings.size());
	mostRecentMoleculetype->improperdihedralbonds.resize(improperbondStrings.size());

	ParseBonds(singlebondStrings, mostRecentMoleculetype->singlebonds, ParseSingleBond, error);
	ParseBonds(pairbondStrings, mostRecentMoleculetype->pairbonds, ParsePairBond, error);
	ParseBonds(anglebondStrings, mostRecentMoleculetype->anglebonds, ParseAngleBond, error);
	ParseBonds(dihedralbondStrings, mostRecentMoleculetype->dihedralbonds, ParseDihedralBond, error);
	ParseBonds(improperbondStrings, mostRecentMoleculetype->improperdihedralbonds, ParseImproperDihedralBond, error);
}
//
//void TopologyFile::ParsePreprocessedFileIntoTopology(const std::string& preprocessedFile) {
//	//topology.defines.insert("FLEXIBLE");// Cant handle gromacs definition of rigid water right now
//
//	std::istringstream file(preprocessedFile);
//
//	TopologySection current_section{ TopologySection::title };
//	TopologySectionGetter getTopolSection{};
//	std::shared_ptr<Moleculetype> mostRecentMoleculetype = nullptr;
//
//	std::string line{};
//
//
//	while (getline(file, line)) {
//		if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {
//
//			// Directives where the directive itself is enough
//			if (current_section == defaults) {
//				// This file is a forcefield. We add it to the includes, and return to parent topol
//				if (forcefieldInclude != std::nullopt)
//					throw std::runtime_error("Trying to include a forcefield, but topology already has 1!");
//
//				// TODO: I dunno wtf this is, maybe not have this at all anymore?
//				//forcefieldInclude.emplace(ForcefieldInclude(fs::path{ includefileName.value_or("forcefield.itp") }));
//			}
//			continue;
//		}
//
//		if (line.empty() || isOnlySpacesAndTabs(line))
//			continue;
//
//		// Check if current line is commented
//		if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) {
//			continue;
//		}	// Only title-sections + atoms reads the comments
//
//		/*if (FileUtils::ChechlineForDefine(line)) {
//			topology.defines.insert(FileUtils::ChechlineForDefine(line).value());
//			continue;
//		}*/
//
//		/*if (FileUtils::ChecklineForIfdefAndSkipIfFound(file, line, topology.defines))
//			continue;*/
//
//
//
//		//if (line[0] == '#') {
//
//		//	if (line.size() > 8 && line.substr(0, 8) == "#include") {
//		//		// take second word, remove "
//		//		std::istringstream iss(line);
//		//		std::string _, pathWithQuotes;
//		//		iss >> _ >> pathWithQuotes;
//		//		if (pathWithQuotes.size() < 3)
//		//			throw std::runtime_error("Include is not formatted as expected: " + line);
//
//		//		std::string filename = pathWithQuotes.substr(1, pathWithQuotes.size() - 2);
//
//		//		// TODO: Check that we havent' already parsed this file
//
//		//		if (filename.find("posre") != std::string::npos) {
//		//			// Do nothing, not yet supported
//		//		}
//		//		else if (filename.find(".itp") != std::string::npos) {
//		//			const fs::path filepath(path.parent_path() / filename);
//		//			if (fs::exists(filepath))
//		//				ParseFileIntoTopology(topology, path.parent_path() / filename, filename);
//		//			else if (fs::exists(FileUtils::GetLimaDir() / "resources/forcefields" / filename))
//		//				ParseFileIntoTopology(topology, FileUtils::GetLimaDir() / "resources/forcefields" / filename, filename);
//		//			else
//		//				throw std::runtime_error(std::format("Could not find file \"{}\" in directory \"{}\"", filename, path.parent_path().string()));
//		//		}
//		//	}
//		//	continue;
//		//}
//
//		// Directives where w eread the contents
//		switch (current_section)
//		{
//		case TopologySection::title:
//			title.append(line + "\n");	// +\n because getline implicitly strips it away.
//			break;
//		case TopologySection::moleculetype:
//		{
//			std::istringstream iss(line);
//			std::string moleculetypename;
//			int nrexcl;
//			iss >> moleculetypename >> nrexcl;
//
//			if (moleculetypename.empty())
//				throw std::runtime_error("Moleculetype name is empty");
//
//			mostRecentMoleculetype = std::make_shared<Moleculetype>();
//			mostRecentMoleculetype->name = moleculetypename;
//			mostRecentMoleculetype->nrexcl = nrexcl;
//
//			//auto nextSection = ParseMoleculetype(file, moleculetype);
//			assert(!moleculetypes.contains(moleculetypename));
//			moleculetypes.insert({ moleculetypename, mostRecentMoleculetype });
//
//			//current_section = nextSection;
//			break;
//		}
//		case TopologySection::_system: {
//			SetSystem(line);
//			break;
//		}
//		case TopologySection::molecules: {
//			std::istringstream iss(line);
//
//			std::string molname;
//			int cnt = 0;
//			iss >> molname >> cnt;
//
//			if (m_system.title == "noSystem")
//				throw std::runtime_error("Molecule section encountered before system section in file: " + path.string());
//			if (!moleculetypes.contains(molname))
//				throw std::runtime_error(std::format("Moleculetype {} not defined before being used in file: {}", molname, path.string()));
//			for (int i = 0; i < cnt; i++)
//				m_system.molecules.emplace_back(MoleculeEntry{ molname, moleculetypes.at(molname) });
//			break;
//		}
//		case TopologySection::atoms:
//		case TopologySection::bonds:
//		case TopologySection::pairs:
//		case TopologySection::angles:
//		case TopologySection::dihedrals:
//		case TopologySection::impropers:
//			if (mostRecentMoleculetype == nullptr)
//				throw std::invalid_argument("Moleculetype not set before parsing atoms/bonds/pairs/angles/dihedrals/impropers");
//			ParseMoleculetypeEntry(current_section, line, mostRecentMoleculetype);
//			break;
//		case TopologySection::atomtypes:
//		case TopologySection::pairtypes:
//		case TopologySection::bondtypes:
//		case TopologySection::constainttypes:
//		case TopologySection::angletypes:
//		case TopologySection::dihedraltypes:
//		case TopologySection::impropertypes:
//			forcefieldInclude->AddEntry(current_section, line);
//			break;
//		default:
//			// Do nothing
//			//throw std::runtime_error("Illegal state");
//			break;
//		}
//	}
//}

TopologyFile::TopologyFile() {}
TopologyFile::TopologyFile(const fs::path& path) : path(path)
{
	if (!(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp"))
		throw std::runtime_error("Expected .top or .itp extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	ParseFileIntoTopology(*this, path);


	/*TimeIt::PrintTaskStats("entry");
	TimeIt::PrintTaskStats("Parsing");*/
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

		if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {
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

	m_system.molecules.emplace_back(MoleculeEntry{ moleculename, moleculetypes.at(moleculename) });
}
void TopologyFile::AppendMoleculetype(const std::shared_ptr<const Moleculetype> moleculetype, std::optional<ForcefieldInclude> inputForcefieldInclude) {
	if (!moleculetypes.contains(moleculetype->name)) {
		moleculetypes.insert({ moleculetype->name, std::make_shared<Moleculetype>(*moleculetype) });	//COPY

		if (inputForcefieldInclude.has_value()) {
			if (!forcefieldInclude.has_value())
				forcefieldInclude.emplace(inputForcefieldInclude.value());
			else
				assert(forcefieldInclude->filename == inputForcefieldInclude->filename);
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
	{
		std::ofstream file(path);
		if (!file.is_open()) {
			throw std::runtime_error(std::format("Failed to open file {}", path.string()));
		}
		
		file << "; " << title << "\n\n";

		// TODO: Have multiple forcefields, just only 1 with the [ defaults ] directive
		if (forcefieldInclude) {
			forcefieldInclude.value().SaveToDir(path.parent_path());
			file << ("#include \"forcefield.itp\"\n");
		}
		file << "\n";

		for (const auto& [_, moleculetype] : moleculetypes) {
			moleculetype->ToFile(path.parent_path());
			file << "#include \"" << moleculetype->includePath.value_or(fs::path(moleculetype->name + ".itp")).string() << "\"\n";
		}
		file << "\n";

		if (m_system.IsInit()) {
			file << "[ system ]\n";
			file << m_system.title << "\n\n";

			file << "[ molecules ]\n";
			for (int i = 0; i < m_system.molecules.size(); i++) {
				std::ostringstream oss;
				int count = 1;
				while (i + 1 < m_system.molecules.size() && m_system.molecules[i].name == m_system.molecules[i + 1].name) {
					count++;
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

	const fs::path path = dir / (name + ".itp");
	
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
		<< std::setw(10) << cgnr
		<< std::setw(10) << std::fixed << std::setprecision(2) << charge;
	if (mass.has_value())
		oss << std::setw(10) << std::fixed << std::setprecision(3) << mass.value();
	oss << '\n';
}

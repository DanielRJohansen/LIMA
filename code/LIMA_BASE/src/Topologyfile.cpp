#include "MDFiles.h" 
#include "Filehandling.h"

#include <format>
#include <algorithm>

using namespace FileUtils;

class TopologySectionGetter {
	int dihedralCount = 0;
	int dihedraltypesCount = 0;
public:
	TopologySection operator()(const std::string& directive) {
		if (directive == "molecules") return molecules;
		if (directive == "moleculetype") return moleculetype;
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

		throw std::runtime_error(std::format("Got unexpected topology directive: {}", directive));
	}
};

inline std::string extractSectionName(const std::string& line) {
	size_t start = line.find('[');
	size_t end = line.find(']', start);
	if (start != std::string::npos && end != std::string::npos) {
		// Extract the text between '[' and ']'
		std::string sectionName = line.substr(start + 1, end - start - 1);
		removeWhitespace(sectionName);
		return sectionName;
	}
	return "";
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

inline fs::path SearchForFile(const fs::path& workDir, const std::string& includeName) {
	// First look relative to the current topology file
	std::optional<fs::path> includePath = _SearchForFile(workDir, includeName);

	// If no value, look in the default includes dir
	if (!includePath.has_value())
		includePath = _SearchForFile(FileUtils::GetLimaDir() / "resources/Slipids", includeName);

	if (!includePath.has_value())
		throw std::runtime_error(std::format("Could not find file \"{}\" in directory \"{}\"", includeName, workDir.string()));

	return includePath.value();
}



template <int n>
inline bool VerifyAllParticlesInBondExists(const std::vector<int>& groIdToLimaId, int ids[n]) {
	for (int i = 0; i < n; i++) {
		if (ids[i] >= groIdToLimaId.size() || groIdToLimaId[ids[i]] == -1)
			return false;
	}
	return true;
}

/// <returns>True if we change section/stay on no section, in which case we should 'continue' to the next line. False otherwise</returns>
inline bool HandleTopologySectionStartAndStop(const std::string& line, TopologySection& currentSection, TopologySectionGetter& sectionGetter) {

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
inline bool isOnlySpacesAndTabs(const std::string& str) {
	return std::all_of(str.begin(), str.end(), [](char c) {
		return c == ' ' || c == '\t';
		});
}






TopologySection TopologyFile::ParseMoleculetype(std::ifstream& file, std::shared_ptr<Moleculetype> moleculetype) {

	TopologySection current_section{ TopologySection::moleculetype };
	TopologySectionGetter getTopolSection{};

	std::vector<int> groIdToLimaId;

	std::string line, atomsSectionName;
	while (getline(file, line)) {
		bool exit = false;
		if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {
			if (current_section != TopologySection::atoms && current_section != TopologySection::bonds && current_section != TopologySection::pairs
				&& current_section != TopologySection::angles && current_section != TopologySection::dihedrals && current_section != TopologySection::impropers)
				return current_section;

			continue;
		}

		if (line.empty() || isOnlySpacesAndTabs(line))
			continue;

		// Check if current line is commented
		if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) {
			continue;
		}	// Only title-sections + atoms reads the comments

		std::istringstream iss(line);


		switch (current_section)
		{
		case TopologySection::atoms:
		{
			if (firstNonspaceCharIs(line, ';')) {
				// TODO: Test for residue or lipid_section in the [1] position of the comment instead

				// Skip the very first line which is the legend
				if (line.find("cgnr") != std::string::npos) {
					break;

				}
				if (line.find("residue") != std::string::npos || line.find("lipid_section") != std::string::npos)
					atomsSectionName= line;
			}
			else {
				TopologyFile::AtomsEntry atom;
				int groId;
				iss >> groId >> atom.type >> atom.resnr >> atom.residue >> atom.atomname >> atom.cgnr >> atom.charge >> atom.mass;

				if (groIdToLimaId.size() < groId + 1)
					groIdToLimaId.resize(groId + 1, -1);
				groIdToLimaId[groId] = moleculetype->atoms.size();
				atom.id = groIdToLimaId[groId];
				moleculetype->atoms.emplace_back(atom);

				if (atomsSectionName!= "") {
					moleculetype->atoms.back().section_name = atomsSectionName;
					atomsSectionName= "";
				}
				if (atom.type.empty() || atom.residue.empty() || atom.atomname.empty())
					throw std::runtime_error("Atom type, residue or atomname is empty");

			}
			break;
		}
		case TopologySection::bonds: {
			TopologyFile::SingleBond singlebond{};
			int groIds[2];
			iss >> groIds[0] >> groIds[1] >> singlebond.funct;
			if (!VerifyAllParticlesInBondExists<2>(groIdToLimaId, groIds))
				break;
			for (int i = 0; i < 2; i++)
				singlebond.ids[i] = groIdToLimaId[groIds[i]];
			singlebond.sourceLine = line;
			moleculetype->singlebonds.emplace_back(singlebond);
			break;
		}
		case TopologySection::pairs: {
			TopologyFile::Pair pair{};
			int groIds[2];
			iss >> groIds[0] >> groIds[1] >> pair.funct;
			if (!VerifyAllParticlesInBondExists<2>(groIdToLimaId, groIds))
				break;
			for (int i = 0; i < 2; i++)
				pair.ids[i] = groIdToLimaId.at(groIds[i]);
			moleculetype->pairs.emplace_back(pair);
			break;
		}
		case TopologySection::angles: {
			TopologyFile::AngleBond angle{};
			int groIds[3];
			iss >> groIds[0] >> groIds[1] >> groIds[2] >> angle.funct;
			if (!VerifyAllParticlesInBondExists<3>(groIdToLimaId, groIds))
				break;
			for (int i = 0; i < 3; i++)
				angle.ids[i] = groIdToLimaId.at(groIds[i]);
			moleculetype->anglebonds.emplace_back(angle);
			break;
		}
		case TopologySection::dihedrals: {
			TopologyFile::DihedralBond dihedral{};
			int groIds[4];
			iss >> groIds[0] >> groIds[1] >> groIds[2] >> groIds[3] >> dihedral.funct;
			if (!VerifyAllParticlesInBondExists<4>(groIdToLimaId, groIds))
				break;
			for (int i = 0; i < 4; i++)
				dihedral.ids[i] = groIdToLimaId.at(groIds[i]);
			moleculetype->dihedralbonds.emplace_back(dihedral);
			break;
		}
		case TopologySection::impropers: {
			TopologyFile::ImproperDihedralBond improper{};
			int groIds[4];
			iss >> groIds[0] >> groIds[1] >> groIds[2] >> groIds[3] >> improper.funct;
			if (!VerifyAllParticlesInBondExists<4>(groIdToLimaId, groIds))
				break;
			for (int i = 0; i < 4; i++)
				improper.ids[i] = groIdToLimaId.at(groIds[i]);
			moleculetype->improperdihedralbonds.emplace_back(improper);
			moleculetype->improperdihedralbonds.back().sourceLine = line;
			break;
		}
		default: {
			// We shouldnt get here
		}
		}
	}
	return current_section;
}







void TopologyFile::ParseFileIntoTopology(TopologyFile& topology, const fs::path& path) {
	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	TopologySection current_section{ TopologySection::title };
	TopologySectionGetter getTopolSection{};

	std::string line{};


	while (getline(file, line)) {
		if (HandleTopologySectionStartAndStop(line, current_section, getTopolSection)) {
			continue;
		}

		if (line.empty() || isOnlySpacesAndTabs(line))
			continue;

		// Check if current line is commented
		if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) {
			continue;
		}	// Only title-sections + atoms reads the comments

		std::istringstream iss(line);


		if (line[0] == '#') {
			if (line.size() > 8 && line.substr(0, 8) == "#include") {
				// take second word, remove "

				std::string _, pathWithQuotes;
				iss >> _ >> pathWithQuotes;
				if (pathWithQuotes.size() < 3)
					throw std::runtime_error("Include is not formatted as expected: " + line);

				std::string filename = pathWithQuotes.substr(1, pathWithQuotes.size() - 2);

				// TODO: Check that we havent' already parsed this file

				if (filename.find("forcefield.itp") != std::string::npos) {
					assert(!topology.forcefieldInclude);
					topology.forcefieldInclude = ForcefieldInclude(filename, path.parent_path());
					//ParseFileIntoTopology(topology, path.parent_path() / filename);
				}
				else if (filename.find("charmm27.ff") != std::string::npos) {
					// This is when including ions.itp or top3p.itp or other structures which we dont yet support!
					// TODO Warn user here!
				}
				else if (filename.find("posre") != std::string::npos) {
					// Do nothing, not yet supported
				}
				else if (filename.find("topol_") != std::string::npos || filename.find(".itp") != std::string::npos) {
					ParseFileIntoTopology(topology, path.parent_path() / filename);
				}
			}
			continue;
		}


		switch (current_section)
		{
		case TopologySection::title:
			topology.title.append(line + "\n");	// +\n because getline implicitly strips it away.
			break;
		case TopologySection::defaults: {

		}
		case TopologySection::moleculetype:
		{
			std::string moleculetypename;
			int nrexcl;
			iss >> moleculetypename >> nrexcl;

			auto moleculetype = std::make_shared<Moleculetype>();
			moleculetype->name = moleculetypename;
			moleculetype->nrexcl = nrexcl;

			auto nextSection = ParseMoleculetype(file, moleculetype);
			assert(!topology.moleculetypes.contains(moleculetypename));
			topology.moleculetypes.insert({ moleculetypename, moleculetype });

			current_section = nextSection;
			break;
		}		
		case TopologySection::_system: {
			topology.SetSystem(line);
			break;
		}
		case TopologySection::molecules: {
			std::string molname;
			int cnt = 0;
			iss >> molname >> cnt;

			if (molname == "SOL")
				continue;

			if (topology.m_system.title == "noSystem")
				throw std::runtime_error("Molecule section encountered before system section in file: " + path.string());
			if (!topology.moleculetypes.contains(molname))
				throw std::runtime_error(std::format("Moleculetype {} not defined before being used in file: {}", molname, path.string()));
			for (int i = 0; i < cnt; i++)
				topology.m_system.molecules.emplace_back(MoleculeEntry{ molname, topology.moleculetypes.at(molname) });
		}
		default:
			// Do nothing
			//throw std::runtime_error("Illegal state");
			break;
		}
	}
}


TopologyFile::TopologyFile() {}
TopologyFile::TopologyFile(const fs::path& path, TopologyFile* parentTop) : path(path)
{
	if (!(path.extension().string() == std::string{ ".top" } || path.extension().string() == ".itp"))
		throw std::runtime_error("Expected .top or .itp extension");
	if (!fs::exists(path))
		throw std::runtime_error(std::format("File \"{}\" was not found", path.string()));

	ParseFileIntoTopology(*this, path);
}


GenericItpFile::GenericItpFile(const fs::path& path) {
	if (path.extension().string() != ".itp") { throw std::runtime_error(std::format("Expected .itp extension with file {}", path.string())); }
	if (!fs::exists(path)) { throw std::runtime_error(std::format("File \"{}\" was not found", path.string())); }

	std::ifstream file;
	file.open(path);
	if (!file.is_open() || file.fail()) {
		throw std::runtime_error(std::format("Failed to open file {}\n", path.string()));
	}

	//lastModificationTimestamp = fs::last_write_time(path);

	//if (UseCachedBinaryFile(path, lastModificationTimestamp) && ENABLE_FILE_CACHING) {
	//	readTopFileFromBinaryCache(path, *this);

	//	// Any include topologies will not be in this particular cached binary, so iterate through them, and load them
	//	for (auto& molecule : molecules.entries) {
	//		assert(molecule.includeTopologyFile == nullptr);
	//		assert(molecule.name != "");

	//		if (includedFiles.count(molecule.name) == 0)
	//			includedFiles.emplace(molecule.name, SearchForFile(path.parent_path(), molecule.name));
	//		molecule.includeTopologyFile = includedFiles.at(molecule.name).Get();
	//	}
	//}
	//else {


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

			sections.at(includes).emplace_back(line);
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

std::optional<std::string> extractStringBetweenQuotationMarks(const std::string& in) {
	size_t first_quote = in.find('"');
	if (first_quote == std::string::npos) {
		return std::nullopt;
	}

	size_t second_quote = in.find('"', first_quote + 1);
	if (second_quote == std::string::npos) {
		return std::nullopt;
	}

	return in.substr(first_quote + 1, second_quote - first_quote - 1);
}

// If a forcefield with the specified name is available relative to the topology's path, then we take that user supplied path and use
// Otherwise we assume the name is simply pointing to a LIMA forcefield.
TopologyFile::ForcefieldInclude::ForcefieldInclude(const std::string& name, const fs::path& ownerDir) :
	name(name),
	path(fs::exists(ownerDir / name)
		? ownerDir / name
		: FileUtils::GetLimaDir() / "resources" / "forcefields" / name)
{}
void TopologyFile::ForcefieldInclude::CopyToDirectory(const fs::path& directory, const fs::path& ownerDir) const {
	if (!fs::is_directory(directory)) {
		throw std::runtime_error(std::format("Directory \"{}\" does not exist", directory.string()));
	}

	// Create the target path for the main file
	const fs::path toplevelForcefieldTargetPath = directory / name;

	if (fs::canonical(toplevelForcefieldTargetPath) == fs::canonical(Path()))
		return; // This forcefield has already been copied to the target location

	// Copy the main file
	if (!fs::exists(toplevelForcefieldTargetPath.parent_path())) {
		fs::create_directories(toplevelForcefieldTargetPath.parent_path()); // Create parent directories if they don't exist
	}
	fs::copy_file(Path(), toplevelForcefieldTargetPath, fs::copy_options::overwrite_existing);


	std::vector<fs::path> subIncludes;
	GenericItpFile ffInclude(Path());

	for (const std::string& subInclude : ffInclude.GetSection(includes)) {
		auto subIncludeName = extractStringBetweenQuotationMarks(subInclude);
		if (!subIncludeName.has_value()) {
			throw std::runtime_error("Could not extract include name from include directive: " + subInclude);
		}
		subIncludes.emplace_back(Path().parent_path() / subIncludeName.value());
	}

	// Copy sub-includes
	for (const auto& include : subIncludes) {
		fs::path includeTargetPath = toplevelForcefieldTargetPath.parent_path() / include.filename();
		if (!fs::exists(includeTargetPath.parent_path())) {
			fs::create_directories(includeTargetPath.parent_path()); // Ensure directory structure exists
		}

		fs::copy_file(include, includeTargetPath, fs::copy_options::overwrite_existing);
	}
}
void TopologyFile::ForcefieldInclude::LoadFullPath(const fs::path& ownerDir) {
	path = fs::exists(ownerDir / name)
		? ownerDir / name
		: FileUtils::GetLimaDir() / "resources" / "forcefields" / name;
}
const fs::path& TopologyFile::ForcefieldInclude::Path() const {
	if (!path.has_value())
		throw std::runtime_error("Path has not been set");
	return path.value();
}






std::optional<fs::path> TopologyFile::GetForcefieldPath() const {
	if (forcefieldInclude)
		return forcefieldInclude->Path();
	/*else if (parentTopology != nullptr)
		return parentTopology->GetForcefieldPath();*/

	return std::nullopt;
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
				forcefieldInclude = inputForcefieldInclude.value();
			else
				assert(forcefieldInclude->name == inputForcefieldInclude->name);
		}
	}
	AppendMolecule(moleculetype->name);
}
//void TopologyFile::AppendMolecule(const MoleculeEntry& molecule) {
//	if (!m_system.IsInit()) {
//		throw std::runtime_error("System is not initialized");
//	}
//
//	if (!moleculetypes.contains(molecule.name)) {
//		moleculetypes.insert({ molecule.name, molecule.moleculetype });
//	}
//
//	if (m_system.molecules.back().name == molecule.name) {
//		m_system.molecules.back().count += molecule.count;
//	}
//	else {
//		m_system.molecules.emplace_back(molecule);
//	}
//}
//void TopologyFile::AppendMolecules(const std::vector<MoleculeEntry>& molecules) {
//	for (const auto& molecule : molecules)
//		AppendMolecule(molecule);
//}

void TopologyFile::printToFile(const std::filesystem::path& path) const {
	const auto ext = path.extension().string();
	if (ext != ".top" && ext != ".itp") { throw std::runtime_error(std::format("Got {} extension, expected [.top/.itp]", ext)); }
	{
		std::ofstream file(path);
		if (!file.is_open()) {
			throw std::runtime_error(std::format("Failed to open file {}", path.string()));
		}
		
		file << title << "\n\n";

		// TODO: Have multiple forcefields, just only 1 with the [ defaults ] directive
		if (forcefieldInclude) {
			forcefieldInclude.value().CopyToDirectory(path.parent_path(), path.parent_path());
			file << ("#include \"" + forcefieldInclude->name.string() + "\"\n");
		}
		file << "\n";

		for (const auto& moleculetype : moleculetypes) {
			moleculetype.second->ToFile(path.parent_path());
			file << "#include \"" << moleculetype.second->name << ".itp\"\n";
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

		file << name << "\n\n";

		file << "[ moleculetype ]\n";
		file << generateLegend({ "name", "nrexcl" }) + "\n";
		file << std::right << std::setw(10) << name << std::setw(10) << nrexcl << "\n\n";
		
		
		
		
		file << "[ atoms ]\n" << generateLegend({ "nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass" }) + "\n";
		file << composeString(atoms);
		file << "[ bonds ]\n" << generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) + "\n";
		file << composeString(singlebonds);
		file << "[ pairs ]\n" << generateLegend({ "ai", "aj", "funct", "c0", "c1", "c2", "c3" }) + "\n";
		file << composeString(pairs);
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
		<< std::setw(10) << std::fixed << std::setprecision(2) << charge
		<< std::setw(10) << std::fixed << std::setprecision(3) << mass
		<< '\n';
}

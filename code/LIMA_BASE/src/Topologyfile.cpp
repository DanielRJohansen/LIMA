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


		if (directive == "settles") return notimplemented;
		if (directive == "exclusions") return notimplemented;
		if (directive == "nonbond_params") return notimplemented; // TODO implement this: https://manual.gromacs.org/current/reference-manual/topologies/parameter-files.html

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

		if (FileUtils::ChecklineForIfdefAndSkipIfFound(file, line))
			continue;

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







void TopologyFile::ParseFileIntoTopology(TopologyFile& topology, const fs::path& path, std::optional<std::string> includefileName) {
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

			// Directives where the directive itself is enough
			if (current_section == defaults) {
				// This file is a forcefield. We add it to the includes, and return to parent topol
				if (topology.forcefieldInclude != std::nullopt)
					throw std::runtime_error("Trying to include a forcefield, but topology already has 1!");

				topology.forcefieldInclude.emplace(ForcefieldInclude(fs::path{includefileName.value_or("forcefield.itp")}));
				//topology.forcefieldInclude = ForcefieldInclude(includefileName.value_or(path.filename().string()), path);				
			}


			continue;
		}

		if (line.empty() || isOnlySpacesAndTabs(line))
			continue;

		// Check if current line is commented
		if (firstNonspaceCharIs(line, commentChar) && current_section != TopologySection::title && current_section != TopologySection::atoms) {
			continue;
		}	// Only title-sections + atoms reads the comments
		
		if (FileUtils::ChecklineForIfdefAndSkipIfFound(file, line))
			continue;

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

		// Directives where w eread the contents
		switch (current_section)
		{
		case TopologySection::title:
			topology.title.append(line + "\n");	// +\n because getline implicitly strips it away.
			break;
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
			break;
		}
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

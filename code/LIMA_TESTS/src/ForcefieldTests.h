#include "Programs.h"
#include "TestUtils.h"





#include <regex>

using namespace TestUtils;








#include <variant>


void ParseForcefieldFromTpr(const std::string& filePath, std::vector<SingleBond::Parameters>& bonds, std::vector<UreyBradley::Parameters>& angles,
    std::vector<DihedralBond::Parameters>& dihedrals, std::vector<ImproperDihedralBond::Parameters>& improperDihedrals, std::vector<std::array<int, 4>>& dihIds) 
{
    std::vector<std::variant<SingleBond::Parameters, AngleBond::Parameters, UreyBradley::Parameters, DihedralBond::Parameters, ImproperDihedralBond::Parameters>> functypes;

    std::ifstream file(filePath);
    std::string line;


    std::regex ntypesRegex(R"(ntypes=(\d+))");
    std::smatch match;
    while (std::getline(file, line)) {
        if (std::regex_search(line, match, ntypesRegex)) {
            int ntypes = std::stoi(match[1].str());
            functypes.resize(ntypes);
            break;
        }
    }
    assert(functypes.size() > 0);




    //std::regex bondRegex(R"(b0A\s*=\s*([\d\.\+e\-]+),\s*cbA\s*=\s*([\d\.\+e\-]+))");
    //std::regex angleRegex(R"(t0A\s*=\s*([\d\.\+e\-]+),\s*ctA\s*=\s*([\d\.\+e\-]+))");
    //std::regex dihedralRegex(R"(phiA\s*=\s*([\d\.\+e\-]+),\s*cpA\s*=\s*([\d\.\+e\-]+),\s*phiB\s*=\s*([\d\.\+e\-]+),\s*cpB\s*=\s*([\d\.\+e\-]+),\s*mult\s*=\s*(\d+))");
    //std::regex improperDihedralRegex(R"(psi0\s*=\s*([\d\.\+e\-]+),\s*kPsi\s*=\s*([\d\.\+e\-]+))");


    std::regex bondRegex(R"(functype\[(\d+)\]=BONDS,\s*b0A\s*=\s*([\d\.\+e\-]+),\s*cbA\s*=\s*([\d\.\+e\-]+),\s*b0B\s*=\s*([\d\.\+e\-]+),\s*cbB\s*=\s*([\d\.\+e\-]+))");
    std::regex angleRegex(R"(functype\[(\d+)\]=ANGLES,\s*t0A\s*=\s*([\d\.\+e\-]+),\s*ctA\s*=\s*([\d\.\+e\-]+))");
    std::regex ureyBradleyRegex(R"(functype\[(\d+)\]=UREY_BRADLEY,\s*thetaA\s*=\s*([\d\.\+e\-]+),\s*kthetaA\s*=\s*([\d\.\+e\-]+),\s*r13A\s*=\s*([\d\.\+e\-]+),\s*kUBA\s*=\s*([\d\.\+e\-]+))");    
    std::regex dihedralRegex(R"(functype\[(\d+)\]=PDIHS,\s*phiA\s*=\s*([\d\.\+e\-]+),\s*cpA\s*=\s*([\d\.\+e\-]+),\s*phiB\s*=\s*([\d\.\+e\-]+),\s*cpB\s*=\s*([\d\.\+e\-]+),\s*mult\s*=\s*(\d+))");
    std::regex improperDihedralRegex(R"(functype\[(\d+)\]=IDIHS,\s*xiA\s*=\s*([\d\.\+e\-]+),\s*cxA\s*=\s*([\d\.\+e\-]+))");



    std::regex bondInstanceRegex(R"(type=(\d+)\s*\(BONDS\))");
    std::regex ureybradleyInstanceRegex(R"(type=(\d+)\s*\(UREY_BRADLEY\))");
    //std::regex dihedralsInstanceRegex(R"(type=(\d+)\s*\(PDIHS\))");
    std::regex dihedralsInstanceRegex(R"(type=(\d+)\s*\(PDIHS\)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+))");
    std::regex improperDihedralsInstanceRegex(R"(type=(\d+)\s*\(IDIHS\))");

    // First find all the types
    while (std::getline(file, line)) {
        std::smatch match;
        if (std::regex_search(line, match, bondRegex)) {
            int index = std::stoi(match[1].str());
            functypes[index] = SingleBond::Parameters{ std::stof(match[2].str()), std::stof(match[3].str())*KILO};
        }
       /* else if (std::regex_search(line, match, angleRegex)) {
            int index = std::stoi(match[1].str());
            functypes[index] = AngleBond::Parameters{ std::stof(match[2].str()), std::stof(match[3].str()) };
        }*/
        else if (std::regex_search(line, match, ureyBradleyRegex)) {
            int index = std::stoi(match[1].str());
            UreyBradley::Parameters params{
                std::stof(match[2].str()), std::stof(match[3].str()), std::stof(match[4].str()), std::stof(match[5].str())
            };
            params.theta0 *= DEG_TO_RAD;
            params.kTheta *= KILO;
            functypes[index] = params;
        }
        else if (std::regex_search(line, match, dihedralRegex)) {
            int index = std::stoi(match[1].str());
            functypes[index] = DihedralBond::Parameters{std::stof(match[2].str()) * DEG_TO_RAD, std::stof(match[3].str()) * KILO, std::stoi(match[6].str())};
        }
        else if (std::regex_search(line, match, improperDihedralRegex)) {
            int index = std::stoi(match[1].str());
            functypes[index] = ImproperDihedralBond::Parameters{ std::stof(match[2].str()) * DEG_TO_RAD, std::stof(match[3].str()) * KILO };
        }

        // Now find the instances and which type they point to
        // We can safely do so in the same loop, as these will come after
        else if (std::regex_search(line, match, bondInstanceRegex)) {
            bonds.emplace_back(std::get<SingleBond::Parameters>(functypes[std::stoi(match[1].str())]));
        }
        else if (std::regex_search(line, match, ureybradleyInstanceRegex)) {
            angles.emplace_back(std::get<UreyBradley::Parameters>(functypes[std::stoi(match[1].str())]));
        }
        else if (std::regex_search(line, match, dihedralsInstanceRegex)) {
			dihedrals.emplace_back(std::get<DihedralBond::Parameters>(functypes[std::stoi(match[1].str())]));
            dihIds.emplace_back(std::array<int,4>{ std::stoi(match[2].str()), std::stoi(match[3].str()) , std::stoi(match[4].str()) , std::stoi(match[5].str()) });
		}
		else if (std::regex_search(line, match, improperDihedralsInstanceRegex)) {
			improperDihedrals.emplace_back(std::get<ImproperDihedralBond::Parameters>(functypes[std::stoi(match[1].str())]));
		}
    }
}

void ParseForcefieldFromItp(
    const std::string& filePath, 
    std::vector<SingleBond::Parameters>& bonds, 
    std::vector<AngleBond::Parameters>& angles,
    std::vector<DihedralBond::Parameters>& dihedrals, 
    std::vector<ImproperDihedralBond::Parameters>& improperDihedrals,
    std::vector<std::array<std::string, 2>>& atomnamesSinglebonds,
    std::vector<std::array<std::string, 3>>& atomnamesAnglebonds,
    std::vector<std::array<std::string, 4>>& atomnamesDihedralbonds,
    std::vector<std::array<std::string, 4>>& atomnamesImproperDihedralbonds)
{
    const GenericItpFile file{ filePath };
    
    int ignoredCount = 0;

    for (const std::string& line : file.GetSection(TopologySection::bondtypes)) {
        std::istringstream iss(line);

        std::array<std::string, 2> atomNames;
        SingleBond::Parameters bondparams{};
        iss >> atomNames[0] >> atomNames[1]
            >> bondparams.b0	
            >> bondparams.kb;	
        bonds.emplace_back(bondparams);
        atomnamesSinglebonds.emplace_back(atomNames);
	}
    for (const std::string& line : file.GetSection(TopologySection::angletypes)) {
		std::istringstream iss(line);

		std::array<std::string, 3> atomNames;
		AngleBond::Parameters angleparams{};
		iss >> atomNames[0] >> atomNames[1] >> atomNames[2]
			>> angleparams.theta_0
			>> angleparams.k_theta;
		angles.emplace_back(angleparams);
		atomnamesAnglebonds.emplace_back(atomNames);
	}
    for (const std::string& line : file.GetSection(TopologySection::dihedraltypes)) {
        std::istringstream iss(line);

        std::array<std::string, 4> atomNames;
        DihedralBond::Parameters dihedralparams{};
        float phi0{};
        float kPhi{};
        float n{};
        iss >> atomNames[0] >> atomNames[1] >> atomNames[2] >> atomNames[3]
            >> phi0
            >> kPhi
            >> n;
        dihedralparams.phi_0 = phi0;
        dihedralparams.k_phi = kPhi;
        dihedralparams.n = n;

        if (kPhi == 0) {
			ignoredCount++;
			continue;
		}

        dihedrals.emplace_back(dihedralparams);
        atomnamesDihedralbonds.emplace_back(atomNames);
    }
    for (const std::string& line : file.GetSection(TopologySection::impropertypes)) {
        std::istringstream iss(line);

        std::array<std::string, 4> atomNames;
        ImproperDihedralBond::Parameters improperDihedralparams{};
        iss >> atomNames[0] >> atomNames[1] >> atomNames[2] >> atomNames[3]
            >> improperDihedralparams.psi_0
            >> improperDihedralparams.k_psi;
        improperDihedrals.emplace_back(improperDihedralparams);
    }
}


LimaUnittestResult TestLimaChosesSameBondparametersAsGromacs(EnvMode envmode) 
{
	//Programs::GetForcefieldParams(GroFile{ TestUtils::simulations_dir / "T4Lysozyme/molecule/conf.gro" }, 
 //    TopologyFile{ TestUtils::simulations_dir / "T4Lysozyme/molecule/topol.top" },
	//	TestUtils::simulations_dir / "Forcefieldtests");
	

    std::vector<SingleBond::Parameters> bondparamsGromacs, bondparamsLima;
    std::vector<UreyBradley::Parameters> angleparamsGromacs;
    std::vector<AngleBond::Parameters> angleparamsLima;
    std::vector<DihedralBond::Parameters> dihedralparamsGromacs, dihedralparamsLima;
    std::vector<ImproperDihedralBond::Parameters> improperDihedralparamsGromacs, improperDihedralparamsLima;
    std::vector<std::array<std::string, 2>> atomnamesSinglebonds;
    std::vector<std::array<std::string, 3>> atomnamesAnglebonds;
    std::vector<std::array<std::string, 4>> atomnamesDihedralbonds;
    std::vector<std::array<std::string, 4>> atomnamesImproperDihedralbonds;

    std::vector<std::array<int, 4>> dihIds;


    ParseForcefieldFromTpr((TestUtils::simulations_dir / "Forcefieldtests/tpr_content.txt").string(), bondparamsGromacs, 
        angleparamsGromacs, dihedralparamsGromacs, improperDihedralparamsGromacs, dihIds);

    ParseForcefieldFromItp((TestUtils::simulations_dir / "Forcefieldtests/appliedForcefield.itp").string(), bondparamsLima,
        angleparamsLima, dihedralparamsLima, improperDihedralparamsLima, atomnamesSinglebonds, atomnamesAnglebonds, atomnamesDihedralbonds, atomnamesImproperDihedralbonds);





    TopologyFile top{ TestUtils::simulations_dir / "T4Lysozyme/molecule/topol.top" };
    int index= 0;
    for (auto dih : top.GetAllDihedralbonds()) {
        for (int i = 0; i < 4; i++) {
            float k = dihedralparamsLima[index].k_phi;
            float phi = dihedralparamsLima[index].phi_0;
            float n = dihedralparamsLima[index].n;
            if (dih.ids[i] != dihIds[index][i])
                int c = 0;
        }
        index++;
    }








    const float maxError = 0.01;

    ASSERT(bondparamsGromacs.size() == bondparamsLima.size(), "Singlebond parameters size mismatch");
    ASSERT(bondparamsGromacs.size() == atomnamesSinglebonds.size(), "Singlebond names and params size mismatch");
    for (int i = 0; i < bondparamsGromacs.size(); i++) {
        const auto g = bondparamsGromacs[i];
        const auto l = bondparamsLima[i];

        const std::string errMsg = std::format("Singlebond mismatch at index {}\n\t {} {}\n\t GMX params: {:.4e} {:.4e}\n\tLIMA params: {:.4e} {:.4e}", 
            i, atomnamesSinglebonds[i][0], atomnamesSinglebonds[i][1], g.b0, g.kb, l.b0, l.kb);

        ASSERT(std::abs(l.b0 - g.b0) < 0.001f && std::abs((l.kb - g.kb)/std::max(g.kb, 1.f)) < maxError, errMsg);
    }
    if (envmode == Full) {
        printf("%d Singlebond parameters verified\n", bondparamsGromacs.size());
    }


    ASSERT(angleparamsGromacs.size() == angleparamsLima.size(), "Anglebond parameters size mismatch");
    ASSERT(angleparamsGromacs.size() == atomnamesAnglebonds.size(), "Anglebond names and params size mismatch");
    for (int i = 0; i < angleparamsGromacs.size(); i++) {
		const auto g = angleparamsGromacs[i];
		const auto l = angleparamsLima[i];

		const std::string errMsg = std::format("Anglebond mismatch at index {}\n\t {} {} {}\n\t GMX params: {:.4e} {:.4e}\n\tLIMA params: {:.4e} {:.4e}",
			i, atomnamesAnglebonds[i][0], atomnamesAnglebonds[i][1], atomnamesAnglebonds[i][2], g.theta0, g.kTheta, l.theta_0, l.k_theta);

		ASSERT(std::abs(l.theta_0 - g.theta0) < 0.001f && std::abs((l.k_theta - g.kTheta) / std::max(g.kTheta, 1.f)) < maxError, errMsg);
	}
    if (envmode == Full) {
		printf("%d Anglebond parameters verified\n", angleparamsGromacs.size());
	}


 //   ASSERT(dihedralparamsGromacs.size() == dihedralparamsLima.size(), "Dihedralbond parameters size mismatch");
 //   ASSERT(dihedralparamsGromacs.size() == atomnamesDihedralbonds.size(), "Dihedralbond names and params size mismatch");
 //   for (int i = 0; i < dihedralparamsGromacs.size(); i++) {
	//	const auto g = dihedralparamsGromacs[i];
	//	const auto l = dihedralparamsLima[i];

	//	const std::string errMsg = std::format("Dihedralbond mismatch at index {}\n\t {} {} {} {}\n\t GMX params: {:.4e} {:.4e} {}\n\tLIMA params: {:.4e} {:.4e} {}",
	//		i, atomnamesDihedralbonds[i][0], atomnamesDihedralbonds[i][1], atomnamesDihedralbonds[i][2], atomnamesDihedralbonds[i][3], 
 //           (float)g.phi_0, (float)g.k_phi, (float)g.n, (float)l.phi_0, (float)l.k_phi, (float)l.n);

	//	ASSERT(std::abs((float)l.phi_0 - (float)g.phi_0) < 0.001f 
 //           && std::abs(((float)l.k_phi - (float)g.k_phi) / std::max((float)g.k_phi, 1.f)) < maxError 
 //           && (float)l.n == (float)g.n, errMsg);
	//}
 //   if (envmode == Full) {
 //       printf("%d Dihedralbond parameters verified\n", dihedralparamsGromacs.size());
 //   }


    return LimaUnittestResult{ true, "", envmode == Full };
}
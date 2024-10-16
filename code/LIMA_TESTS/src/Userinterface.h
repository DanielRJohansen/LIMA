#pragma once

#include "TestUtils.h"
#include "Statistics.h"




namespace UserinterfaceTests {
    static const EnvMode envmode = Headless;


    fs::path executable() {
        // Since this is being called by LIMA_TESTS.exe in build/debug|release/code/LIMA_TESTS
        return fs::current_path().parent_path() / "LIMA" / "lima.exe";
    }

    static const fs::path test_data_dir = FileUtils::GetLimaDir() / "tests" / "userinterfacetests";



    // Utility function to run a command and capture its return code
    int run_command(const std::string& command) {
        return std::system(command.c_str());
    }

    // Test 1: Verify that buildmembrane with valid lipids generates output files
    LimaUnittestResult TestBuildMembraneWithValidLipids() {
        // Clean up any existing output files
        fs::remove("membrane.gro");
        fs::remove("membrane.top");

        // Command to run
        std::string command = executable().string() + " buildmembrane -lipids DPPC 60 DOPC 40 -centerz 3.0 -boxsize 10.0";

        // Run the command
        int ret = run_command(command);
        ASSERT(ret == 0, "buildmembrane command failed with return code " + std::to_string(ret));

        // Check that the output files exist
        ASSERT(fs::exists("membrane.gro"), "Output file membrane.gro not found");
        ASSERT(fs::exists("membrane.top"), "Output file membrane.top not found");

        // Optionally, inspect the content of the files (e.g., check the number of lipids)
        // This requires domain-specific knowledge and parsing of the .gro and .top files

        return LimaUnittestResult(true, "", false);
    }

    // Test 2: Verify that buildmembrane fails with invalid lipid percentages (e.g., percentages do not sum to 100)
    LimaUnittestResult TestBuildMembraneWithInvalidPercentages() {
        // Command to run with invalid percentages
        std::string command = executable().string() + " buildmembrane -lipids DPPC 70 DOPC 40 -centerz 3.0 -boxsize 10.0";

        // Run the command
        int ret = run_command(command);

        // Expecting a non-zero return code
        ASSERT(ret != 0, "buildmembrane command should have failed due to invalid lipid percentages");

        // Check that output files are not generated
        ASSERT(!fs::exists("membrane.gro"), "membrane.gro should not have been created");
        ASSERT(!fs::exists("membrane.top"), "membrane.top should not have been created");

        return LimaUnittestResult(true, "", false);
    }

    // Test 3: Verify that buildmembrane fails with an odd number of arguments after -lipids
    LimaUnittestResult TestBuildMembraneWithIncompleteLipidsArgs() {
        // Command with incomplete lipid arguments
        std::string command = executable().string() + " buildmembrane -lipids DPPC 70 -centerz 3.0 -boxsize 10.0";

        // Run the command
        int ret = run_command(command);

        // Expecting a non-zero return code
        ASSERT(ret != 0, "buildmembrane command should have failed due to incomplete lipid arguments");

        // Check that output files are not generated
        ASSERT(!fs::exists("membrane.gro"), "membrane.gro should not have been created");
        ASSERT(!fs::exists("membrane.top"), "membrane.top should not have been created");

        return LimaUnittestResult(true, "", false);
    }

    // Test 4: Verify that buildmembrane handles invalid lipid names gracefully
    LimaUnittestResult TestBuildMembraneWithInvalidLipidNames() {
        // Command with invalid lipid names
        std::string command = executable().string() + " buildmembrane -lipids INVALID_LIPID 100 -centerz 3.0 -boxsize 10.0";

        // Run the command
        int ret = run_command(command);

        // Expecting a non-zero return code
        ASSERT(ret != 0, "buildmembrane command should have failed due to invalid lipid names");

        // Check that output files are not generated
        ASSERT(!fs::exists("membrane.gro"), "membrane.gro should not have been created");
        ASSERT(!fs::exists("membrane.top"), "membrane.top should not have been created");

        return LimaUnittestResult(true, "", false);
    }

    // Test 5: Verify that the -help option displays help text and exits
    LimaUnittestResult TestBuildMembraneHelpOption() {
        // Command to display help
        std::string command = executable().string() + " buildmembrane -help";

        // Run the command
        int ret = run_command(command);

        // Expecting a zero return code
        ASSERT(ret == 0, "buildmembrane -help command failed with return code " + std::to_string(ret));

        // Optionally, capture and check the output if possible

        return LimaUnittestResult(true, "", false);
    }

    // Test 6: Verify that buildmembrane correctly sets the box size in the output files
    LimaUnittestResult TestBuildMembraneBoxSizeAndCenterz() {
        // Clean up any existing output files
        fs::remove("membrane.gro");
        fs::remove("membrane.top");

        // Command with specific box size
        std::string command = executable().string() + " buildmembrane -lipids DPPC 100 -boxsize 6.0 -centerz 3.0";

        // Run the command
        int ret = run_command(command);
        ASSERT(ret == 0, "buildmembrane command failed with return code " + std::to_string(ret));

        // Check that the output files exist
        ASSERT(fs::exists("membrane.gro"), "Output file membrane.gro not found");

        // Open membrane.gro and check the box size
        GroFile grofile{ "membrane.gro" };

        ASSERT(grofile.box_size == Float3(6.f), "Box size does not match expected value");

        // I guess go through .gro file and calc the avg z value and compare to centerz    
        std::vector<Float3> positions(grofile.atoms.size());
        for (const auto& a : grofile.atoms) {
            positions.push_back(a.position);
        }
        const float avg_z = Statistics::Mean(positions).z;

        ASSERT(std::abs(avg_z - 3.0) < 1e-2, "Average z value does not match expected value");

        return LimaUnittestResult(true, "", false);
    }


    // Main function to run all tests
    LimaUnittestResult TestBuildmembranesInterface(EnvMode _) {
        //LimaUnittestResult result;

        //result = TestBuildMembraneWithValidLipids();
        //if (!result.success) return result;

        auto result = TestBuildMembraneWithInvalidPercentages();
        if (!result.success) return result;

        result = TestBuildMembraneWithIncompleteLipidsArgs();
        if (!result.success) return result;

        result = TestBuildMembraneWithInvalidLipidNames();
        if (!result.success) return result;

        result = TestBuildMembraneHelpOption();
        if (!result.success) return result;

        result = TestBuildMembraneBoxSizeAndCenterz();
        if (!result.success) return result;

        return result;
    }
}
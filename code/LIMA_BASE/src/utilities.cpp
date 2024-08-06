#include "Utilities.h"
#include <filesystem>
#include <iostream>
#include "Printer.h"

using namespace LIMA_Print;
using std::string;

LimaLogger::LimaLogger(LogMode lm, EnvMode em, const string& name, const fs::path& workfolder)
    : logmode(lm)
    , envmode{em}
    , enable_logging(workfolder !="")
    , log_dir((workfolder / "logs/").string())
{
    if (enable_logging) {
        std::filesystem::create_directories(log_dir);
        const auto logFilePath = log_dir + name + "_log.txt";
        logFile.open(logFilePath, std::ios::out);
    }    
}

LimaLogger::~LimaLogger() {
    if (enable_logging && logFile.is_open()) {
        logFile.close();
    }
}

void LimaLogger::startSection(const std::string& input)
{
    if (envmode != Headless) {
        printH2(input, true, false);
    }

    logToFile(input);
}

void LimaLogger::print(const std::string& input, const bool log) 
{
    if (envmode != Headless) {
        string input_copy = input;
        if (logmode == compact) {
            if (clear_next) { clearLine(); }
            if (!input_copy.empty() && input_copy.back() == '\n') {
                input_copy.back() = ' ';
                clear_next = true;  // Only clear of this was end of line
            }
        }
        std::cout << input_copy;
    }


    if (log) {
        logToFile(input);
    }
}

void LimaLogger::finishSection(const string& str) {
    if (envmode != Headless) {
        if (logmode == compact) {
            std::cout << "\n";
        }

        printH2(str, false, true);
    }
    
    logToFile(str);
}

void LimaLogger::logToFile(const string& str)
{
	if (!enable_logging) { return; }
	if (!logFile.is_open()) {
		std::cerr << "Error: Log file is not open" << std::endl;
		return;
	}

	logFile << str << std::endl;
}

void LimaLogger::clearLine() {
    //int a = static_cast<int>(std::cerr.tellp()) - 1;
    //std::cout << "\r" << std::flush;
    //std::cout << std::string(std::max(0, static_cast<int>(std::cerr.tellp()) - 1), ' ') << "\r" << std::flush;
    std::cout << "\033[2K\r" << std::flush;
    clear_next = false;
}



/// <summary>
/// Implementation of SutherlandHodgman algorithm to find the intersection of two convex hulls
/// https://github.com/Alamot/code-snippets/blob/master/graphics/SutherlandHodgman/Linux/SutherlandHodgman.cpp
/// </summary>
/// <param name="ch1"></param>
/// <param name="ch2"></param>
/// <returns></returns>
ConvexHull LAL::FindIntersectionConvexhullFrom2Convexhulls(const ConvexHull& ch1, const ConvexHull& ch2) {
    float D1, D2 = 0;
    ConvexHull polygon = ConvexHull(ch1);
    const float EPSILON = 0.00001f;

    const auto& clippingPlanes = ch2.GetFacets();

    for (unsigned int c = 0; c < clippingPlanes.size(); c++) 
    {
        ConvexHull clippedPolygon{};

        std::vector<Float3> points = polygon.GetVertices();
        for (unsigned int i = 0; i < points.size() - 1; i++) {
            D1 = clippingPlanes[c].distance(points[i]);
            D2 = clippingPlanes[c].distance(points[i + 1]);
            if ((D1 <= 0) && (D2 <= 0))
            {
                //clippedPolygon.add(points[i + 1]);
                clippedPolygon.Add(points[i + 1]);
            }
            else if ((D1 > 0) && ((D2 > -EPSILON) && (D2 < EPSILON)))
            {
                //clippedPolygon.add(points[i + 1]);
                clippedPolygon.Add(points[i + 1]);

            }
            else if (((D1 > -EPSILON) && (D1 < EPSILON)) && (D2 > 0))
            {
                continue;
            }
            else if ((D1 <= 0) && (D2 > 0))
            {
                clippedPolygon.Add(clippingPlanes[c].intersectionPoint(points[i], points[i + 1]));
            }
            else if ((D1 > 0) && (D2 <= 0))
            {                
                clippedPolygon.Add(clippingPlanes[c].intersectionPoint(points[i], points[i + 1]));
                clippedPolygon.Add(points[i + 1]);
            }
        }
        polygon = clippedPolygon; // keep on working with the new polygon
    }
    return polygon;
}
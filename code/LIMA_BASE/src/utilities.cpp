#include "Utilities.h"
#include <filesystem>
#include <iostream>
#include "Printer.h"

using namespace LIMA_Print;
namespace fs = std::filesystem;

LimaLogger::LimaLogger(LogMode lm, EnvMode em, const std::string& name, const fs::path& workfolder)
    : logmode(lm)
    , envmode{em}
    //, enable_logging(workfolder !="")
    , enable_logging(false)
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
        std::string input_copy = input;
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

void LimaLogger::finishSection(const std::string& str) {
    if (envmode != Headless) {
        if (logmode == compact) {
            std::cout << "\n";
        }

        printH2(str, false, true);
    }
    
    logToFile(str);
}

void LimaLogger::logToFile(const std::string& str)
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

//float LAL::LargestDiff(const Float3 queryPoint, const std::span<Float3>& points) {
//	float maxDiff = 0.f;
//    for (const Float3& p : points) {
//        maxDiff = std::max(maxDiff, (p - queryPoint).len()); // OPTIM: compute lenSq, and only len for the index with smallest lenSq
//    }
//
//	return maxDiff;
//}



std::string StringUtils::FormatTime(
    std::chrono::duration<double> duration,
    int decimalPlacesBeforePoint,
    int decimalPlacesAfterPoint
) {
    using Seconds = std::chrono::duration<double>;

    struct Unit {
        Seconds scale;
        const char* suffix;
    };

    static constexpr std::array<Unit, 9> units{
        Unit{ std::chrono::duration<double, std::micro>{1.0}, "us" },
        Unit{ std::chrono::duration<double, std::milli>{1.0}, "ms" },
        Unit{ std::chrono::seconds{1},                         "s" },
        Unit{ std::chrono::minutes{1},                         "min" },
        Unit{ std::chrono::hours{1},                           "hr" },
        Unit{ std::chrono::hours{24},                          "days" },
        Unit{ std::chrono::hours{24 * 7},                      "weeks" },
        Unit{ std::chrono::duration<double, std::ratio<2629746>>{1.0}, "months" }, // 30.44 days
        Unit{ std::chrono::duration<double, std::ratio<31557600>>{1.0}, "years" }   // 365.25 days
    };

    decimalPlacesBeforePoint = std::max(0, decimalPlacesBeforePoint);
    decimalPlacesAfterPoint = std::max(0, decimalPlacesAfterPoint);

    const auto absDuration = duration >= Seconds::zero() ? duration : -duration;

    const Unit* selectedUnit = &units.front();
    for (const auto& unit : units) {
        if (absDuration >= unit.scale) {
            selectedUnit = &unit;
        }
        else {
            break;
        }
    }

    const double value = (duration / selectedUnit->scale);

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(decimalPlacesAfterPoint) << value;

    std::string formatted = oss.str();

    if (decimalPlacesBeforePoint > 0) {
        const std::size_t signOffset = !formatted.empty() && formatted[0] == '-' ? 1u : 0u;
        const std::size_t dotPos = formatted.find('.');
        const std::size_t integerEnd = dotPos == std::string::npos ? formatted.size() : dotPos;
        const std::size_t integerDigits = integerEnd - signOffset;

        if (integerDigits < static_cast<std::size_t>(decimalPlacesBeforePoint)) {
            formatted.insert(
                signOffset,
                static_cast<std::size_t>(decimalPlacesBeforePoint) - integerDigits,
                '0'
            );
        }
    }

    return std::format("{} [{}]", formatted, selectedUnit->suffix);
}

std::vector<std::string> StringUtils::SplitWords(std::string_view line) {
    	std::istringstream input{ std::string(line) };
    	std::vector<std::string> result;
    	for (std::string word; input >> word;) result.push_back(std::move(word));
    	return result;
    }
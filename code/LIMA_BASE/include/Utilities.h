#pragma once

#include <assert.h>
#include <string>
#include <fstream>

namespace LIMA_UTILS {

	/*template <typename T>
	assertEqual(T a, T b) */
	static int roundUp(int numToRound, int multiple)
	{
		assert(multiple);
		return ((numToRound + multiple - 1) / multiple) * multiple;
	}
}




class LimaLogger {
public:
    enum LogMode {
        normal,
        compact
    };
    LimaLogger(const LimaLogger& log) = delete;
    LimaLogger(const LogMode mode, const std::string& name, const std::string& folder="");
    ~LimaLogger();

    void print(const std::string& input, bool log=true);
    void finishSection();

private:
    LogMode mode;
    //std::string logFilePath;
    std::ofstream logFile;
    const bool enable_logging;

    void clearLine();
    bool clear_next = false;
};
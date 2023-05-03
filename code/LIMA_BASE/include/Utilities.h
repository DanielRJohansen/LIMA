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
    LimaLogger(const LimaLogger&) {};
    LimaLogger(const LogMode mode, const std::string& name, const std::string& workfolder="");
    ~LimaLogger();

    void print(const std::string& input, bool log=true);
    void finishSection();
    
    template <typename T>
    static void printToFile(std::string file_path, std::vector<T> data) {
        FILE* file;

        if (!fopen_s(&file, file_path, "wb")) {
            fwrite(data.data(), sizeof(T), data.size(), file);
            fclose(file);
        }
    }


private:
    LogMode mode;
    //std::string logFilePath;
    std::ofstream logFile;
    const bool enable_logging{ false };

    void clearLine();
    bool clear_next = false;
};
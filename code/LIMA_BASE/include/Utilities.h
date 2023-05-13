#pragma once

#include <assert.h>
#include <string>
#include <fstream>
#include <vector>

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

    LimaLogger(const LogMode mode, const std::string& name, const std::string& workfolder="");
    ~LimaLogger();

    void print(const std::string& input, bool log=true);
    void finishSection();
    
    template <typename T>
    void printToFile(const std::string& filename, const std::vector<T>& data) const {
        FILE* file;

        const std::string filepath = (log_dir + filename);
        if (!fopen_s(&file, filepath.c_str(), "wb")) {
            fwrite(data.data(), sizeof(T), data.size(), file);
            fclose(file);
        }
    }


private:
    LogMode mode;
    //std::string logFilePath;
    const std::string log_dir;
    std::ofstream logFile;
    const bool enable_logging{ false };

    void clearLine();
    bool clear_next = false;
};


namespace LAL {


    template<typename T>
    class optional {
    private:
        bool hasValue;
        T _value;

    public:
        optional() : hasValue(false), _value() {}

        optional(const T& value) : hasValue(true), _value(value) {}

        optional(T&& value) : hasValue(true), _value(std::move(value)) {}

        optional& operator=(const T& value) {
            this->_value = value;
            hasValue = true;
            return *this;
        }

        optional& operator=(T&& value) {
            this->_value = std::move(value);
            hasValue = true;
            return *this;
        }

        T value() const {
            if (!hasValue) {
                throw std::exception("No value present in optional.");
            }
            return _value;
        }

        explicit operator bool() const {
            return hasValue;
        }
    };    
}
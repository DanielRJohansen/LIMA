#pragma once
//
//#include <chrono>
//#include <string>
//
//
//
//class TimeIt {
//public:
//	TimeIt(const std::string& taskName = "Task", bool printUponDestruction=false)
//		: printUponDestruction(printUponDestruction), taskName(taskName), start(std::chrono::high_resolution_clock::now()), manuallyStopped(false) {
//		end = start;
//	}
//
//	std::chrono::nanoseconds GetTiming() {
//		return end - start;
//	}
//
//	std::chrono::nanoseconds stop() {
//		if (!manuallyStopped) {
//			end = std::chrono::high_resolution_clock::now();
//			manuallyStopped = true;
//		}
//		return GetTiming();
//	}
//
//	std::chrono::milliseconds elapsed() const {
//		auto currentTime = std::chrono::high_resolution_clock::now();
//		return std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - start);
//	}
//
//	~TimeIt() {
//		if (!manuallyStopped) {
//			end = std::chrono::high_resolution_clock::now();
//		}		
//
//		if (printUponDestruction) {
//			auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
//
//			if (elapsed < std::chrono::milliseconds(2)) {
//				std::cout << taskName << " took " << elapsed.count() << " microseconds.\n";
//			}
//			else {
//				std::cout << taskName << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds.\n";
//			}
//		}
//	}
//
//private:
//	const std::string taskName;
//	const std::chrono::time_point<std::chrono::high_resolution_clock> start;
//	std::chrono::time_point<std::chrono::high_resolution_clock> end;
//	bool manuallyStopped;
//	const bool printUponDestruction;
//};

#include <iostream>
#include <chrono>
#include <string>
#include <unordered_map>
#include <mutex>

class TimeIt {
public:
    TimeIt(const std::string& taskName = "Task", bool printUponDestruction = false);

    std::chrono::nanoseconds GetTiming() const;

    std::chrono::nanoseconds stop();

    std::chrono::milliseconds elapsed() const;

    ~TimeIt();

    static void PrintTaskStats(const std::string& taskName);

private:
    const std::string taskName;
    const std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    bool manuallyStopped;
    const bool printUponDestruction;

    struct TaskRecord {
        std::chrono::nanoseconds totalTime = std::chrono::nanoseconds(0);
        size_t count = 0;
    };

    static std::unordered_map<std::string, TaskRecord> taskRecords;
    static std::mutex mutex_;

    void updateRecord();
};

//std::unordered_map<std::string, TimeIt::TaskRecord> TimeIt::taskRecords;
//std::mutex TimeIt::mutex_;

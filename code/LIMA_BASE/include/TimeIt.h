#pragma once

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
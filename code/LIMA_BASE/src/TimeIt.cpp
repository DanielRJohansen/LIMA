#include "TimeIt.h"

#include <iostream>
#include <chrono>
#include <string>
#include <unordered_map>
#include <mutex>


TimeIt::TimeIt(const std::string& taskName, bool printUponDestruction)
	: taskName(taskName),
	printUponDestruction(printUponDestruction),
	start(std::chrono::high_resolution_clock::now()),
	manuallyStopped(false) 
{
	end = start;
}

std::chrono::nanoseconds TimeIt::GetTiming() const {
	return end - start;
}

std::chrono::nanoseconds TimeIt::stop() {
	if (!manuallyStopped) {
		end = std::chrono::high_resolution_clock::now();
		manuallyStopped = true;
		updateRecord();
	}
	return GetTiming();
}

std::chrono::milliseconds TimeIt::elapsed() const {
	auto currentTime = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - start);
}

TimeIt::~TimeIt() {
	if (!manuallyStopped) {
		end = std::chrono::high_resolution_clock::now();
		updateRecord();
	}

	if (printUponDestruction) {
		auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		if (elapsed < std::chrono::milliseconds(2)) {
			std::cout << taskName << " took " << elapsed.count() << " microseconds.\n";
		}
		else {
			std::cout << taskName << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " milliseconds.\n";
		}
	}
}

void TimeIt::PrintTaskStats(const std::string& taskName) {
	std::lock_guard<std::mutex> lock(mutex_);
	if (taskRecords.find(taskName) != taskRecords.end()) {
		auto& record = taskRecords[taskName];
		double avgTime = static_cast<double>(record.totalTime.count()) / record.count;

		// Convert to milliseconds or microseconds for readability
		auto totalTimeInMicroseconds = std::chrono::duration_cast<std::chrono::microseconds>(record.totalTime);
		auto totalTimeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(record.totalTime);

		if (totalTimeInMilliseconds.count() > 0) {
			// Use milliseconds if the total time is greater than 1 millisecond
			std::cout << "Task \"" << taskName << "\" - Total Time: " << totalTimeInMilliseconds.count()
				<< " milliseconds, Average Time: " << avgTime / 1e6 << " milliseconds.\n";
		}
		else {
			// Otherwise, use microseconds
			std::cout << "Task \"" << taskName << "\" - Total Time: " << totalTimeInMicroseconds.count()
				<< " microseconds, Average Time: " << avgTime / 1e3 << " microseconds.\n";
		}
	}
	else {
		std::cout << "No records found for task \"" << taskName << "\".\n";
	}
}


void TimeIt::updateRecord() {
	auto elapsedTime = end - start;
	std::lock_guard<std::mutex> lock(mutex_);
	taskRecords[taskName].totalTime += elapsedTime;
	taskRecords[taskName].count++;
}


std::unordered_map<std::string, TimeIt::TaskRecord> TimeIt::taskRecords;
std::mutex TimeIt::mutex_;

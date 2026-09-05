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

std::chrono::duration<double> TimeIt::Elapsed() const {
	auto currentTime = std::chrono::high_resolution_clock::now();
	return currentTime - start;
}

// Returns the time with 2 decimals in either ms or s
std::string TimeIt::ElapsedPretty() const {
	const auto t = elapsed();
	const auto ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(t).count();

	if (ms < 1000.0)
		return std::format("{:.2f} [ms]", ms);

	const double s = ms / 1000.0;
	return std::format("{:.2f} [s]", s);
}



void TimeIt::PrintTaskStats(const std::string& taskName, const TaskRecord& record) {
	const double totalUs = std::chrono::duration<double, std::micro>(record.totalTime).count();
	const bool useMs = totalUs >= 1000.0;
	const double scale = useMs ? 1e-3 : 1.0;
	const char* unit = useMs ? "ms" : "us";

	std::cout << std::format("Task \"{}\" - Calls: {}, Total: {:.3f} {}, Average: {:.3f} {}\n",
		taskName, record.count, totalUs * scale, unit, totalUs * scale / record.count, unit);
}

void TimeIt::PrintTaskStats(const std::string& taskName) {
	std::lock_guard<std::mutex> lock(mutex_);
	if (taskRecords.find(taskName) != taskRecords.end()) {
		PrintTaskStats(taskName, taskRecords[taskName]);
	}
	else {
		std::cout << "No records found for task \"" << taskName << "\".\n";
	}
}

void TimeIt::PrintAllTaskStats() {
	std::lock_guard<std::mutex> lock(mutex_);
	for (auto & [taskName, record] : taskRecords) {
		PrintTaskStats(taskName, record);
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

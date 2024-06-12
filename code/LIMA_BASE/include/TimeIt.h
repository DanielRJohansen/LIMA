#pragma once

#include <chrono>
#include <string>



class TimeIt {
public:
	TimeIt(const std::string& taskName = "Task", bool printUponDestruction=false)
		: printUponDestruction(printUponDestruction), taskName(taskName), start(std::chrono::high_resolution_clock::now()), manuallyStopped(false) {
		end = start;
	}

	std::chrono::nanoseconds GetTiming() {
		return end - start;
	}

	std::chrono::nanoseconds stop() {
		if (!manuallyStopped) {
			end = std::chrono::high_resolution_clock::now();
			manuallyStopped = true;
		}
		return GetTiming();
	}

	~TimeIt() {
		if (!manuallyStopped) {
			end = std::chrono::high_resolution_clock::now();
		}
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		if (printUponDestruction)
			std::cout << taskName << " took " << elapsed << " milliseconds.\n";
	}

private:
	const std::string taskName;
	const std::chrono::time_point<std::chrono::high_resolution_clock> start;
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	bool manuallyStopped;
	const bool printUponDestruction;
};
#pragma once

#include "LimaTypes.cuh"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <vector>

struct XDRFILE;

// A GROMACS .xtc or .trr trajectory, read one frame at a time. Opening indexes where each frame starts, so any frame
// can be read directly, for scrubbing back and forth without holding the trajectory in memory
class TrajectoryFile {
public:
	struct Frame {
		std::vector<Float3> positions;	// [nm] In the atom order of the structure file
		Float3 boxSize{};				// [nm]
		int step = 0;
		float time = 0.f;				// [ps]
	};

	explicit TrajectoryFile(const std::filesystem::path& path);
	~TrajectoryFile();
	TrajectoryFile(const TrajectoryFile&) = delete;
	TrajectoryFile& operator=(const TrajectoryFile&) = delete;

	int AtomCount() const { return atomCount; }
	int FrameCount() const { return static_cast<int>(frameOffsets.size()); }
	Frame Read(int frame);

private:
	XDRFILE* file = nullptr;
	bool isXtc = false;
	int atomCount = 0;
	std::vector<int64_t> frameOffsets;
	std::vector<float> buffer;	// x, y, z per atom

	bool ReadNext(int& step, float& time, float box[3][3]);
};

// Playing a trajectory in the display, in place of a live simulation
struct TrajectoryPlayback {
	std::unique_ptr<TrajectoryFile> file;
	int frame = 0;
	int shownFrame = -1;	// The frame whose positions the display has, -1 for none yet
	bool playing = true;
	float framesPerSecond = 20.f;
	float frameTime = 0.f;	// [ps] Of the shown frame
	std::chrono::steady_clock::time_point lastAdvance{};
};

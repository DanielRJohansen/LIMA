#include "TrajectoryFile.h"

#include "Format.h"
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

// From xdr_seek.h, which lacks extern "C" for use from C++
extern "C" {
	int64_t xdr_tell(XDRFILE* xd);
	int xdr_seek(XDRFILE* xd, int64_t pos, int whence);
}

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

TrajectoryFile::TrajectoryFile(const std::filesystem::path& path)
{
	std::string extension = path.extension().string();
	std::ranges::transform(extension, extension.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	if (extension != ".xtc" && extension != ".trr")
		throw std::runtime_error(Lima::Format("Trajectories must be .xtc or .trr, got {}", path.string()));
	isXtc = extension == ".xtc";

	std::string pathString = path.string();
	const int result = isXtc ? read_xtc_natoms(pathString.data(), &atomCount) : read_trr_natoms(pathString.data(), &atomCount);
	if (result != exdrOK || atomCount <= 0)
		throw std::runtime_error(Lima::Format("Could not read the trajectory {}", path.string()));

	file = xdrfile_open(pathString.c_str(), "r");
	if (!file)
		throw std::runtime_error(Lima::Format("Could not open the trajectory {}", path.string()));
	buffer.resize(3 * static_cast<std::size_t>(atomCount));

	int step = 0;
	float time = 0.f;
	float box[3][3]{};
	for (int64_t offset = xdr_tell(file); ReadNext(step, time, box); offset = xdr_tell(file))
		frameOffsets.push_back(offset);
	if (frameOffsets.empty())
		throw std::runtime_error(Lima::Format("The trajectory {} has no frames", path.string()));
}

TrajectoryFile::~TrajectoryFile()
{
	if (file)
		xdrfile_close(file);
}

bool TrajectoryFile::ReadNext(int& step, float& time, float box[3][3])
{
	auto* positions = reinterpret_cast<rvec*>(buffer.data());
	if (isXtc) {
		float precision = 0.f;
		return read_xtc(file, atomCount, &step, &time, box, positions, &precision) == exdrOK;
	}
	float lambda = 0.f;
	return read_trr(file, atomCount, &step, &time, &lambda, box, positions, nullptr, nullptr) == exdrOK;
}

TrajectoryFile::Frame TrajectoryFile::Read(int frame)
{
	if (frame < 0 || frame >= FrameCount())
		throw std::runtime_error(Lima::Format("Trajectory frame {} is out of range, it has {} frames", frame, FrameCount()));

	Frame result;
	float box[3][3]{};
	if (xdr_seek(file, frameOffsets[frame], SEEK_SET) != exdrOK || !ReadNext(result.step, result.time, box))
		throw std::runtime_error(Lima::Format("Could not read trajectory frame {}", frame));
	result.boxSize = Float3{ box[0][0], box[1][1], box[2][2] };
	result.positions.resize(atomCount);
	for (int i = 0; i < atomCount; ++i)
		result.positions[i] = Float3{ buffer[3 * i], buffer[3 * i + 1], buffer[3 * i + 2] };
	return result;
}

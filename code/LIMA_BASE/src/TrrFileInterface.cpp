#include "MDFiles.h"
#include "Format.h"

#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"

#include <algorithm>
#include <cctype>
#include <limits>


MDFiles::TrajectoryWriter::TrajectoryWriter(const fs::path& path, int nAtoms) : path(path), nAtoms(nAtoms) {
	std::string extension = path.extension().string();
	std::ranges::transform(extension, extension.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	if (extension != ".trr" && extension != ".xtc")
		throw std::runtime_error(Lima::Format("Trajectories are written as .trr or .xtc, got {}", path.string()));
	xtc = extension == ".xtc";
	if (!path.parent_path().empty())
		fs::create_directories(path.parent_path());
	file = xdrfile_open(path.string().c_str(), "w");
	if (!file)
		throw std::runtime_error(Lima::Format("Failed to open {} for writing", path.string()));
}

MDFiles::TrajectoryWriter::~TrajectoryWriter() {
	if (file)
		xdrfile_close(file);
}

void MDFiles::TrajectoryWriter::WriteFrame(int64_t step, float timePs, Float3 box, const Float3* positions) {
	if (step > std::numeric_limits<int>::max())
		throw std::runtime_error(Lima::Format("Step {} does not fit the 32-bit step of a {} frame", step, xtc ? "xtc" : "trr"));
	matrix boxMatrix = {
		{ box.x, 0.f, 0.f },
		{ 0.f, box.y, 0.f },
		{ 0.f, 0.f, box.z }
	};
	// The xdrfile API takes non-const pointers, but only reads the positions
	rvec* x = reinterpret_cast<rvec*>(const_cast<Float3*>(positions));
	constexpr float xtcPrecision = 1000.f;	// [1/nm]
	const int result = xtc
		? write_xtc(file, nAtoms, static_cast<int>(step), timePs, boxMatrix, x, xtcPrecision)
		: write_trr(file, nAtoms, static_cast<int>(step), timePs, 0.f, boxMatrix, x, nullptr, nullptr);
	if (result != exdrOK)
		throw std::runtime_error(Lima::Format("Failed to write a frame to {}", path.string()));
}

void MDFiles::Dump(Trajectory& trajectory, const fs::path& path) {
	TrajectoryWriter writer{ path, trajectory.nAtoms };
	if (trajectory.nAtoms == 0)
		return;
	for (int frame = 0; frame < trajectory.nFrames; frame++)
		writer.WriteFrame(static_cast<int64_t>(trajectory.stepInterval) * frame, trajectory.frameSpacingPs * frame,
			trajectory.boxSize, trajectory.GetFrame(frame));
}

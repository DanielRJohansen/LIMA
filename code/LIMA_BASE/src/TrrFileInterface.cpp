//#include "ChemfilesInterface.h"
//#include "EngineUtils.cuh"
#include "MDFiles.h"


#include "xdrfile.h"
#include "xdrfile_trr.h"


void MDFiles::Dump(Trajectory& trajectory, const fs::path& path) {
	XDRFILE* file = xdrfile_open(path.string().c_str(), "w");

	if (trajectory.nAtoms == 0 || trajectory.nFrames == 0)
		return;

	
	
	

	matrix box = {
		{trajectory.boxSize.x, 0.f, 0.f},
		{0.f, trajectory.boxSize.y, 0.f},
		{ 0.f, 0.f, trajectory.boxSize.z }
	};


	//rvec* xData = new rvec[trajectory.nAtoms];	
	for (int step = 0; step < trajectory.nFrames; step++) {
		rvec* x = (rvec*)trajectory.GetFrame(step);
		const float t = trajectory.lambda * step;

		write_trr(file, trajectory.nAtoms, step, t, trajectory.lambda, box, x, nullptr, nullptr);
	}

	xdrfile_close(file);


	// Check that the file was written correctly
	{
		int natoms;
		unsigned long nframes;
		read_trr_nframes(&path.string()[0], &nframes);
		read_trr_natoms(&path.string()[0], &natoms);

		if (nframes != trajectory.nFrames) {
			throw std::runtime_error("Number of frames in TRR file does not match the number of frames in the simulation");
		}
		if (natoms != trajectory.nAtoms) {
			throw std::runtime_error("Number of atoms in TRR file does not match the number of atoms in the simulation");
		}

	}
}
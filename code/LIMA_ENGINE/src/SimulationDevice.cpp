#include "SimulationDevice.cuh"


// Assumes "this" is a ptr to a device object
void BoxDevice::CopyDataToHost(Box& boxHost) {
	BoxDevice boxtemp;
	cudaMemcpy(&boxtemp, this, sizeof(BoxDevice), cudaMemcpyDeviceToHost);

	assert(boxHost.compounds.size() == boxtemp.boxparams.n_compounds);
	cudaMemcpy(boxHost.compounds.data(), boxtemp.compounds, sizeof(Compound) * boxHost.compounds.size(), cudaMemcpyDeviceToHost);
	cudaMemcpy(boxHost.solvents.data(), boxtemp.solvents, sizeof(Solvent) * boxHost.solvents.size(), cudaMemcpyDeviceToHost);

	boxHost.solventblockgrid_circularqueue->CopyDataFromDevice(boxtemp.solventblockgrid_circularqueue);
	boxHost.compoundcoordsCircularQueue->CopyDataFromDevice(boxtemp.compoundcoordsCircularQueue);
}

// Free this immediately after calling this function
void BoxDevice::DeleteBox() {
	BoxDevice boxtemp;
	cudaMemcpy(&boxtemp, this, sizeof(BoxDevice), cudaMemcpyDeviceToHost);

	cudaFree(boxtemp.compounds);
	cudaFree(boxtemp.compoundcoordsCircularQueue);
	cudaFree(boxtemp.solvents);
	cudaFree(boxtemp.solventblockgrid_circularqueue);
	cudaFree(boxtemp.bridge_bundle);
	cudaFree(boxtemp.bpLUTs);
}

BoxDevice* MakeBox(const Box& box) {
	BoxDevice boxTemp;
	cudaMalloc(&boxTemp.compounds, sizeof(Compound) * box.boxparams.n_compounds);
	cudaMalloc(&boxTemp.solvents, sizeof(Solvent) * box.boxparams.n_solvents);
	cudaMalloc(&boxTemp.bridge_bundle, sizeof(CompoundBridgeBundleCompact));

	boxTemp.boxparams = box.boxparams;
	boxTemp.uniformElectricField = box.uniformElectricField;

	cudaMemcpy(boxTemp.compounds, box.compounds.data(), sizeof(Compound) * box.boxparams.n_compounds, cudaMemcpyHostToDevice);
	boxTemp.compoundcoordsCircularQueue = box.compoundcoordsCircularQueue->CopyToDevice();
	cudaMemcpy(boxTemp.solvents, box.solvents.data(), sizeof(Solvent) * box.boxparams.n_solvents, cudaMemcpyHostToDevice);
	boxTemp.solventblockgrid_circularqueue = box.solventblockgrid_circularqueue->CopyToDevice();
	cudaMemcpy(boxTemp.bridge_bundle, box.bridge_bundle.get(), sizeof(CompoundBridgeBundleCompact), cudaMemcpyHostToDevice);
	boxTemp.bpLUTs = GenericCopyToDevice(box.bpLutCollection);

	BoxDevice* devPtr;
	cudaMallocManaged(&devPtr, sizeof(BoxDevice));
	cudaMemcpy(devPtr, &boxTemp, sizeof(BoxDevice), cudaMemcpyHostToDevice);

	LIMA_UTILS::genericErrorCheck("Error during MakeBox");


	return devPtr;
}
















SimulationDevice::SimulationDevice(const SimParams& params_host, Box* box_host)
{
	// Allocate structures for keeping track of solvents and compounds
	compound_grid = BoxGrid::MallocOnDevice<CompoundGridNode>(box_host->boxparams.boxSize);
	cudaMallocManaged(&compound_neighborlists, sizeof(NeighborList) * MAX_COMPOUNDS);

	cudaMallocManaged(&transfermodule_array, sizeof(SolventBlockTransfermodule) * BoxGrid::BlocksTotal(BoxGrid::NodesPerDim(box_host->boxparams.boxSize)));

	{
		SimSignals temp{};
		genericCopyToDevice(temp, &signals, 1);
	}

	genericCopyToDevice(params_host, &params, 1);

	{
		DatabuffersDevice databuffersTemp(box_host->boxparams.total_particles_upperbound, box_host->boxparams.n_compounds, params_host.data_logging_interval);
		databuffers = GenericCopyToDevice(&databuffersTemp, 1);
	}	

	box = MakeBox(*box_host);

	chargeGrid = BoxGrid::MallocOnDevice<Electrostatics::ChargeNode>(box_host->boxparams.boxSize);
	chargeGridChargeSums = BoxGrid::MallocOnDevice<float>(box_host->boxparams.boxSize);
	chargeGridOutputForceAndPot = BoxGrid::MallocOnDevice<ForceAndPotential>(box_host->boxparams.boxSize);
}

void SimulationDevice::deleteMembers() {
	box->DeleteBox();
	cudaFree(box);

	databuffers->freeMembers();
	cudaFree(databuffers);

	cudaFree(params);

	cudaFree(compound_grid);
	cudaFree(compound_neighborlists);

	//cudaFree(charge_octtree);
	cudaFree(chargeGrid);
	cudaFree(chargeGridChargeSums);
	cudaFree(chargeGridOutputForceAndPot);
}
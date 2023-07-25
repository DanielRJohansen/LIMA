#include "LIMA_ENGINE/include/Engine.cuh"



__host__ static TemperaturPackage getBoxTemperature(Simulation* simulation, ForceField_NB& forcefield_host) {
	TemperaturPackage package{};

	const uint64_t step = simulation->getStep() - 1;
	const uint64_t step_offset_a = step * simulation->boxparams_host.total_particles_upperbound;
	const uint64_t step_offset_b = (step - 1) * simulation->boxparams_host.total_particles_upperbound;
	const uint64_t solvent_offset = MAX_COMPOUND_PARTICLES * simulation->boxparams_host.n_compounds;

	const auto dt = simulation->simparams_host.constparams.dt;

	long double sum_kinE_compound = 0.;	// [J/mol]
	for (int compound_id = 0; compound_id < simulation->boxparams_host.n_compounds; compound_id++) {
		uint64_t compound_offset = compound_id * MAX_COMPOUND_PARTICLES;
		for (int pid = 0; pid < simulation->compounds_host[compound_id].n_particles; pid++) {	// i gotta move this somewhere else....
			const Float3 posa = simulation->traj_buffer->getCompoundparticleDatapoint(compound_id, pid, step);
			const Float3 posb = simulation->traj_buffer->getCompoundparticleDatapoint(compound_id, pid, step - 1);
			const float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[simulation->compounds_host[compound_id].atom_types[pid]].mass, dt / NANO_TO_LIMA);

			package.max_kinE_compound = std::max(package.max_kinE_compound, kinE);
			sum_kinE_compound += kinE;
		}
	}

	long double sum_kinE_solvents = 0.;	// [J/mol]
	for (int i = 0; i < simulation->boxparams_host.n_solvents; i++) {
		const Float3 posa = simulation->traj_buffer->getSolventparticleDatapoint(i, step);
		const Float3 posb = simulation->traj_buffer->getSolventparticleDatapoint(i, step-1);
		const float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[0].mass, dt / NANO_TO_LIMA);

		package.max_kinE_solvent = std::max(package.max_kinE_solvent, kinE);
		sum_kinE_solvents += static_cast<float>(kinE);
	}
	package.avg_kinE_solvent = static_cast<float>(sum_kinE_solvents / static_cast<long double>(simulation->boxparams_host.n_solvents));

	const long double total_kinE = (sum_kinE_compound + sum_kinE_solvents) / AVOGADROSNUMBER;
	package.temperature = EngineUtils::kineticEnergyToTemperature(total_kinE, simulation->extraparams.total_particles);

	return package;
}



void Engine::handleBoxtemp() {
	const float target_temp = 310.f;				// [k]
	const TemperaturPackage temp_package = getBoxTemperature(simulation, forcefield_host);
	const float temp = temp_package.temperature;

	simulation->temperature_buffer.push_back(temp_package.temperature);

	simulation->temperature = temp;	// For display :)

	
	if constexpr(!APPLY_THERMOSTAT) {
		// So we avoid dividing by 0
		const float temp_safe = temp == 0.f ? 1 : temp;
		float temp_scalar = target_temp / temp_safe;

		// I just added this to not change any temperatures too rapidly
		temp_scalar = std::clamp(temp_scalar, 1.f - MAX_THERMOSTAT_SCALER, 1.f + MAX_THERMOSTAT_SCALER);
		// Apply 1/n scalar for n steps.

		simulation->sim_dev->params->thermostat_scalar = temp_scalar;
	}	
}



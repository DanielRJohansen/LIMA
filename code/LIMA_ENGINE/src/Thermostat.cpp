#include "Engine.cuh"



__host__ static TemperaturPackage getBoxTemperature(Simulation* simulation, ForceField_NB& forcefield_host) {
	TemperaturPackage package{};

	const uint64_t step = simulation->getStep() - 1;
	const uint64_t step_offset_a = step * simulation->total_particles_upperbound;
	const uint64_t step_offset_b = (step - 1) * simulation->total_particles_upperbound;
	const uint64_t solvent_offset = MAX_COMPOUND_PARTICLES * simulation->n_compounds;

	long double sum_kinE_compound = 0.;
	for (uint64_t c = 0; c < simulation->n_compounds; c++) {
		uint64_t compound_offset = c * MAX_COMPOUND_PARTICLES;
		for (uint64_t i = 0; i < simulation->compounds_host[c].n_particles; i++) {	// i gotta move this somewhere else....
			Float3 posa = simulation->traj_buffer[i + compound_offset + step_offset_a];
			Float3 posb = simulation->traj_buffer[i + compound_offset + step_offset_b];
			float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[simulation->compounds_host[c].atom_types[i]].mass, simulation->dt * 1.f);			// Doesnt work, use forcefield_host!!

			package.max_kinE_compound = std::max(package.max_kinE_compound, kinE);
			sum_kinE_compound += kinE;
		}
	}
	package.avg_kinE_compound = static_cast<float>(sum_kinE_compound / static_cast<long double>(simulation->total_compound_particles));

	long double sum_kinE_solvents = 0.;
	for (int i = 0; i < simulation->n_solvents; i++) {
		Float3 posa = simulation->traj_buffer[i + solvent_offset + step_offset_a];
		Float3 posb = simulation->traj_buffer[i + solvent_offset + step_offset_b];
		float kinE = EngineUtils::calcKineticEnergy(&posa, &posb, forcefield_host.particle_parameters[0].mass, simulation->dt * 1.f);

		package.max_kinE_solvent = std::max(package.max_kinE_solvent, kinE);
		sum_kinE_solvents += static_cast<float>(kinE);
	}
	package.avg_kinE_solvent = static_cast<float>(sum_kinE_solvents / static_cast<long double>(simulation->n_solvents));


	float avg_kinE = static_cast<float>((sum_kinE_compound + sum_kinE_solvents) / static_cast<long double>(simulation->total_particles));
	package.temperature = avg_kinE * 2.f / (3.f * 8.3145f);

	return package;
}



void Engine::handleBoxtemp() {
	const float target_temp = 310.f;				// [k]
	const TemperaturPackage temp_package = getBoxTemperature(simulation, forcefield_host);
	float temp = temp_package.temperature;


	simulation->temperature_buffer[simulation->n_temp_values++] = temp_package.temperature;

	// So we avoid dividing by 0
	float temp_safe = temp == 0.f ? 1 : temp;
	float temp_scalar = target_temp / temp_safe;

	// I just added this to not change any temperatures too rapidly
	temp_scalar = std::clamp(temp_scalar, 1.f-MAX_THERMOSTAT_SCALER, 1.f + MAX_THERMOSTAT_SCALER);
	// Apply 1/n scalar for n steps.


	uint64_t step = simulation->getStep();
	if (step >= FIRST_TEMPERATURE_PRINT_STEP && PRINT_TEMP && std::abs(temp - target_temp) > 100.f) {
		LIMA_Printer::printNameValuePairs(
			"Temperature", temp,
			"Avg kinE sol", temp_package.avg_kinE_solvent,
			"Avg kinE comp", temp_package.avg_kinE_compound,
			"Max kinE sol", temp_package.max_kinE_solvent,
			"Max kinE comp", temp_package.max_kinE_compound
		);
	}
	
	if (APPLY_THERMOSTAT && step >= FIRST_THERMOSTAT_APPLICATION_STEP) {
		simulation->box->thermostat_scalar = temp_scalar;

		if (temp_scalar != temp_scalar) {//} || abs(temp_scalar) == "inf") {
			printf("Scalar: %f\n", simulation->box->thermostat_scalar);
			exit(0);
		}
	}
}



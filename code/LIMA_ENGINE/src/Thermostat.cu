#include "Engine.cuh"
#include "EngineUtils.cuh"
#include "SimulationDevice.cuh"

#include <algorithm>


static TemperaturPackage getBoxTemperature(Simulation* simulation) {
	TemperaturPackage package{};

	const uint64_t step = simulation->getStep() - 1;	// We haven't loaded data for current step onto host yet.
	const auto entryindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(step);


	long double sum_kinE_compound = 0.;	// [J/mol]
	for (int compound_id = 0; compound_id < simulation->boxparams_host.n_compounds; compound_id++) {
		for (int pid = 0; pid < simulation->compounds_host[compound_id].n_particles; pid++) {	// i gotta move this somewhere else....
			const float mass = simulation->forcefield->getNBForcefieldRef().particle_parameters[simulation->compounds_host[compound_id].atom_types[pid]].mass;
			const float velocity = simulation->vel_buffer->getCompoundparticleDatapointAtIndex(compound_id, pid, entryindex);
			const float kinE = EngineUtils::calcKineticEnergy(velocity, mass);

			package.max_kinE_compound = std::max(package.max_kinE_compound, kinE);
			sum_kinE_compound += kinE;
		}
	}

	long double sum_kinE_solvents = 0.;	// [J/mol]
	for (int solvent_id = 0; solvent_id < simulation->boxparams_host.n_solvents; solvent_id++) {
		const float mass = simulation->forcefield->getNBForcefieldRef().particle_parameters[ATOMTYPE_SOLVENT].mass;
		const float velocity = simulation->vel_buffer->getSolventparticleDatapointAtIndex(solvent_id, entryindex);
		const float kinE = EngineUtils::calcKineticEnergy(velocity, mass);

		package.max_kinE_solvent = std::max(package.max_kinE_solvent, kinE);
		sum_kinE_solvents += static_cast<float>(kinE);
	}
	package.avg_kinE_solvent = static_cast<float>(sum_kinE_solvents / static_cast<long double>(simulation->boxparams_host.n_solvents));

	const long double total_kinE = (sum_kinE_compound + sum_kinE_solvents) / AVOGADROSNUMBER;
	package.temperature = EngineUtils::kineticEnergyToTemperature(total_kinE, simulation->boxparams_host.total_particles);

	return package;
}



void Engine::handleBoxtemp() {
	const float target_temp = 310.f;				// [k]
	const TemperaturPackage temp_package = getBoxTemperature(simulation.get());
	const float temp = temp_package.temperature;

	simulation->temperature_buffer.push_back(temp_package.temperature);

	simulation->temperature = temp;	// For display :)

	
	if constexpr(APPLY_THERMOSTAT) {
		// So we avoid dividing by 0
		const float temp_safe = temp == 0.f ? 1 : temp;
		float temp_scalar = target_temp / temp_safe;

		
		// I just added this to not change any temperatures too rapidly. However in EM we can go faster, and should so we reach goal temperature before sim starts
		const float max_scalar = simulation->simparams_host.em_variant ? MAX_THERMOSTAT_SCALER * 10.f : MAX_THERMOSTAT_SCALER;
		temp_scalar = std::clamp(temp_scalar, 1.f - max_scalar, 1.f + max_scalar);
		
		
		// Apply 1/n scalar for n steps.

		sim_dev->signals->thermostat_scalar = temp_scalar;	// UNSAFE
	}	
}



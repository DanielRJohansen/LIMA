#include "Engine.cuh"
#include "PhysicsUtils.cuh"
#include <algorithm>
#include <numeric>

float GetBoxTemperature(Simulation* simulation) {

	const uint64_t step = simulation->getStep() - 1;	// We haven't loaded data for current step onto host yet.
	const auto entryindex = LIMALOGSYSTEM::getMostRecentDataentryIndex(step, simulation->simparams_host.data_logging_interval);


	long double sum_kinE_compound = 0.;	// [J/mol]
	for (int compound_id = 0; compound_id < simulation->box_host->boxparams.n_compounds; compound_id++) {
		for (int pid = 0; pid < simulation->box_host->compounds[compound_id].n_particles; pid++) {	// i gotta move this somewhere else....
			const float mass = simulation->forcefield.particle_parameters[simulation->box_host->compounds[compound_id].atom_types[pid]].mass;
			const float velocity = simulation->vel_buffer->getCompoundparticleDatapointAtIndex(compound_id, pid, entryindex);
			const float kinE = PhysicsUtils::calcKineticEnergy(velocity, mass);

			sum_kinE_compound += kinE;
		}
	}

	const float solventMass = simulation->forcefield.particle_parameters[ATOMTYPE_SOLVENT].mass;
	double sumKineSolvents = std::transform_reduce(
		&simulation->vel_buffer->getSolventparticleDatapointAtIndex(0, entryindex), 
		&simulation->vel_buffer->getSolventparticleDatapointAtIndex(0, entryindex) + simulation->box_host->boxparams.n_solvents,
		0.0, std::plus<>(), [solventMass](const float& velocity) {
		return PhysicsUtils::calcKineticEnergy(velocity, solventMass);
	});

	const long double totalKinEnergyAveraged = (sum_kinE_compound + sumKineSolvents) / static_cast<long double>(simulation->box_host->boxparams.total_particles);
	return PhysicsUtils::kineticEnergyToTemperature(totalKinEnergyAveraged);
}



float Engine::HandleBoxtemp() {
	const float temperature = GetBoxTemperature(simulation.get());

	simulation->temperature_buffer.push_back(temperature);
	runstatus.current_temperature = temperature;

	if (simulation->simparams_host.apply_thermostat) {		
		const float target_temp = simulation->simparams_host.em_variant ? 150.f : 310.f;				// [k]
		const float temp_safe = temperature == 0.f ? 1 : temperature;// So we avoid dividing by 0
		float temp_scalar = target_temp / temp_safe;

		
		// I just added this to not change any temperatures too rapidly. However in EM we can go faster, and should so we reach goal temperature before sim starts
		const float MAX_THERMOSTAT_SCALER = 0.001f / static_cast<float>(simulation->simparams_host.steps_per_temperature_measurement);	// change vel by 0.1% over NSTEPS
		const float max_scalar = simulation->simparams_host.em_variant ? MAX_THERMOSTAT_SCALER * 10.f : MAX_THERMOSTAT_SCALER;
		temp_scalar = std::clamp(temp_scalar, 1.f - max_scalar, 1.f + max_scalar);
		
		
		// Apply 1/n scalar for n steps.
		return temp_scalar;
		//sim_dev->signals->thermostat_scalar = temp_scalar;	// UNSAFE
	}
	return 1.f;
}



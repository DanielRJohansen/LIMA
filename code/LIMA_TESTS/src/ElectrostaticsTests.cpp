#include "Tests.h"
#include "Format.h"

#include "TestUtils.h"
#include "Environment.h"
#include "LimaTypes.cuh"
#include "LimaPositionSystem.cuh"
#include "Statistics.h"
#include "PhysicsUtils.cuh"
#include "PlotUtils.h"

#include <iostream>
#include <string>
#include <map>

namespace ElectrostaticsTests {
	using namespace TestUtils;


	TestRoutine CoulombForceSanityCheck(Environment&, EnvMode envmode) {
		const float calcedForce = PhysicsUtils::CalcCoulumbForce(1.f*elementaryChargeToKiloCoulombPerMole, 1.f*elementaryChargeToKiloCoulombPerMole, Float3{ 1.f, 0.f, 0.f }).len(); // [1/l N / mol]
		const float expectedForce = 2.307078e-10 * AVOGADROSNUMBER * NANO;  // [J/mol/nm] https://www.omnicalculator.com/physics/coulombs-law

		if (std::abs(calcedForce - expectedForce) / expectedForce >= 0.0001f)
			co_return LimaUnittestResult{ false, Lima::Format("Expected {:.2e} Actual {:.2e}", expectedForce, calcedForce), envmode == Full };
		// TODO: add potE to this also
		co_return LimaUnittestResult{ true, "Success", envmode == Full};
	}

	//static ForceEnergy CalcImmediateMirrorForceEnergy(const Float3& diff, const float chargeProduct, Float3 boxSize) {
	//	ForceEnergy forceEnergy{};
	//	for (int zOffset = -1; zOffset <= 1; zOffset += 1) {
	//		for (int yOffset = -1; yOffset <= 1; yOffset += 1) {
	//			for (int xOffset = -1; xOffset <= 1; xOffset += 1) {
	//				if (xOffset == 0 && yOffset == 0 && zOffset == 0)
	//					continue;

	//				const Float3 diffMirror = diff + Float3{ xOffset, yOffset, zOffset } *boxSize;
	//				forceEnergy.force += PhysicsUtils::CalcCoulumbForce(chargeProduct, 1, diffMirror);
	//				forceEnergy.potE += PhysicsUtils::CalcCoulumbPotential(chargeProduct, 1, diffMirror.len()) * 0.5f;
	//			}
	//		}
	//	}
	//		return forceEnergy;
	//}




	Float3 GetPositionOfParticleRelativeToSelfUsingTheWierdLogicOfTheKernel(const Float3& posOtherAbs, const NodeIndex nodeindexSelf) {
		// We cant use hyperdist, since the engine will hyperpos the block, and not individual particles
		//BoundaryConditionPublic::applyHyperposNM(posSelf, posOther, sim->simParams.box_size, PBC);
		// Instead we do this bullshit. Figure out the relative nodeindex compared to the self nodeindex
		// 
		// The final Coulomb force is calculated using the position from the second-to-last step, thus -2 not -1
		const NodeIndex nodeindexOther = LIMAPOSITIONSYSTEM::PositionToNodeIndexNM(posOtherAbs);
		const Float3 posOtherRel = posOtherAbs - LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindexOther);

		NodeIndex nodeindexOtherHyper = nodeindexOther;
		BoundaryConditionPublic::applyHyperpos(nodeindexSelf, nodeindexOtherHyper, Int3(3, 3, 3), PBC);
		const NodeIndex nodeindexOfOtherRelativeToSelf = nodeindexOtherHyper - nodeindexSelf;


		if (nodeindexOfOtherRelativeToSelf.largestMagnitudeElement() > 1)
			throw std::runtime_error("Expected to get a relative nodeindex that was immediately adjacent to nodeindexSelf");		

		const Float3 posOtherRelativeToSelf = (posOtherRel + LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindexOfOtherRelativeToSelf));
		return posOtherRelativeToSelf;
	}


	LimaUnittestResult TestAttractiveParticlesInteractingWithESandLJ(EnvMode envmode) {
		const fs::path work_folder = AutomatedTestsDir() / "Pool/";
		const int nSteps = 1000;



		SimParams params{};
		params.n_steps = nSteps;
		params.enable_electrostatics = true;
		params.data_logging_interval = 1;
		params.cutoff_nm = 2.f;
		Environment& environment = Environment::Get();
		SimulationJob job;
		job.workDir = work_folder;
		job.simParams = params;
		job.mode = envmode;
		job.postprocess = SimAnalysis::AnalyzeEnergy;
		job.preprocess = [](GroFile& grofile, TopologyFile&, SimParams&) {
			grofile.box_size = Float3{ 8.f, 4.f, 4.f };
			grofile.atoms[0].position = Float3{ 1.f, 1.5f, 1.5f };
			grofile.atoms[1].position = Float3{ 2.f, 1.5f, 1.5f };
		};
		auto result = environment.Submit(std::move(job)).Get();

		const float actualVC = result.analysis->variance_coefficient;
		const float maxVC = 1e-3;
		ASSERT(actualVC < maxVC, Lima::Format("VC {:.3e} / {:.3e}", actualVC, maxVC));

		return LimaUnittestResult{ true, "", envmode == Full };
	}

	void MakeChargeParticlesSim(
		GroFile& grofile, TopologyFile& topfile, const fs::path& workDir,
		const float boxLen, const AtomsSelection& atomsSelection, float particlesPerNm3) {
		grofile.m_path = workDir / "conf.gro";
		grofile.box_size = Float3{ boxLen };
		topfile.SetSystem("MySystem");
		topfile.path = workDir / "topol.top";

		//MDFiles::SimulationFilesCollection simfiles(env.getWorkdir());
		for (const auto& atom : atomsSelection) {
			auto moltype = std::make_shared<TopologyFile::Moleculetype>( atom.atomtype.atomname, 3);
			moltype->atoms.push_back(atom.atomtype);
			topfile.moleculetypes.insert({ atom.atomtype.atomname, moltype });
		}
		//simfiles.topfile->SetSystem("ElectroStatic Field Test");
		SimulationBuilder::DistributeParticlesInBox(grofile, topfile, atomsSelection, 0.24f, particlesPerNm3);

		// Overwrite the forcefield
		topfile.forcefieldInclude.emplace("lima_custom_forcefield.itp");
		topfile.forcefieldInclude->contents = std::move(GenericItpFile(FileUtils::GetLimaDir() / "resources" / "forcefields" / "lima_custom_forcefield.itp"));

		grofile.title = "ElectroStatic Field Test";
		topfile.title = "ElectroStatic Field Test";
		grofile.printToFile();
		topfile.printToFile();
	}

	TestRoutine TestChargedParticlesVelocityInUniformElectricField(
		Environment& environment, EnvMode envmode) {
		const fs::path workDir = AutomatedTestsDir() / "ElectrostaticField";
		AtomsSelection atoms{
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx1", 0, -1.f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx2", 0, -.5f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx3", 0, -0.f, 10.f}, 40},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx4", 0, 0.5f, 10.f}, 15},
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt2", 0, "lxx", "lx5", 0, 1.f, 10.f},  15}
			};

		SimParams simparams;
		simparams.dt = 0.2f * FEMTO_TO_NANO;
		simparams.coloring_method = ColoringMethod::Charge;
		simparams.data_logging_interval = 1;
		simparams.snf_select.insert(HorizontalChargeField);
		SimulationJob job;
		job.workDir = workDir;
		job.grofile.emplace();
		job.topfile.emplace();
		job.simParams = simparams;
		job.mode = envmode;
		job.preprocess = [workDir, atoms = std::move(atoms)](
			GroFile& grofile, TopologyFile& topfile, SimParams&) {
			MakeChargeParticlesSim(grofile, topfile, workDir, 7.f, atoms, 5.f);
		};
		job.configureSimulation = [](Simulation& simulation) {
			simulation.box->uniformElectricField =
				UniformElectricField{ Float3{-1.f, 0.f, 0.f }, 12.f };
		};
		auto result = co_await environment.Submit(std::move(job));

		if (envmode == Full)
			TestUtils::CompareForces1To1(workDir, *result.simulation, false);

		auto& sim = result.simulation;

		std::map<float, std::vector<float>> velDistributions;

		// Go through each particle in each compound, and assert that their velocities are as we expect in this horizontal electric field
		// TODO!!!
		//for (int cid = 0; cid < sim->box->boxparams.n_compounds; cid++) {
		//	const auto& compound = sim->box->compounds[cid];
		//	const auto& compoundInterimState = sim->box->pclusterInterimStates[cid];
		for (int pcId = 0; pcId < sim->box->persistentClusters.size(); pcId++){
			for (int pid = 0; pid < PersistentCluster::maxParticles; pid++) {
				if (!sim->box->persistentClusters[pcId].pqd[pid].Valid())
					continue;
				const float charge = sim->box->persistentClusters[pcId].pqd[pid].params.charge;
				const float velHorizontal = sim->box->pclusterInterimStates[pcId].vels_prev[pid].x;

				//const float velHorizontal = compoundInterimState.vels_prev[pid].x;

				velDistributions[charge].push_back(velHorizontal);
			}
		}

		if (envmode == Full) {
			for (const auto& pair : velDistributions) {
				const int charge = pair.first;
				const auto& velocities = pair.second;
				std::cout << "Charge: " << charge << " | Mean Velocity: " << Statistics::Mean(velocities) << " | Standard Deviation: " << Statistics::StdDev(velocities) << std::endl;
			}
		}

		std::vector<float> x;
		std::vector<float> y;
		for (const auto& pair : velDistributions) {
			x.insert(x.end(), pair.second.size(), pair.first);
			y.insert(y.end(), pair.second.begin(), pair.second.end());
		}

		const auto [slope, intercept] = Statistics::linearFit(x, y);

		if (slope >= 0.f) {
			std::string errorMsg = Lima::Format("Slope of velocity distribution should be negative, but got {:.4f} ", slope);
			co_return LimaUnittestResult{ false, errorMsg, envmode == Full };
		}
		if (std::abs(intercept) > 50.f) {
			std::string errorMsg = Lima::Format("Intercept of velocity distribution should be close to 0, but got {:.2f}",intercept);
			co_return LimaUnittestResult{ false, errorMsg, envmode == Full };
		}

		const float r2 = Statistics::calculateR2(x, y, slope, intercept);
		if (std::isnan(r2))
			co_return LimaUnittestResult{ false, "R2 value is nan", envmode == Full };
		if (r2 < 0.5f) {
			//std::string errorMsg = "R2 value " + std::to_string(r2) + " of velocity distribution should be close to 1";
			std::string errorMsg = Lima::Format("R2 value {:.2f} of velocity distribution should be close to 1", r2);
			co_return LimaUnittestResult{ false, errorMsg, envmode == Full };
		}

		co_return LimaUnittestResult{ true, Lima::Format("R2 Value: {:.2f}", r2), envmode == Full};
	}

	//static LimaUnittestResult TestElectrostaticsManyParticles(EnvMode envmode) {
	//	MakeChargeParticlesSim("ShortrangeElectrostaticsCompoundOnly", 5.f,
	//		AtomsSelection{
	//			{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "lxx", "lxx", 0, 5.f, 12.011}, 100}, // by naming the residue lxx we let these particles be full molecules, instead of tinymols, is that ideal? Does it matter? 
	//		},
	//		2.f // TODO: If we set this density to 32 as it should be, the result diverge too much. I should look into that later. And do a similar stresstest for a simple LJ system
	//	);

	//	const int nSteps = 1000;

	//	SimParams simparams;
	//	simparams.n_steps = nSteps;
	//	simparams.dt = 1.f * FEMTO_TO_NANO;
	//	simparams.coloring_method = ColoringMethod::Charge;
	//	simparams.data_logging_interval = 1;
	//	simparams.stepsPerNlistupdate = 1;
	//	simparams.enable_electrostatics = true;
	//	simparams.cutoff_nm = 2.f;
	//	auto env = basicSetup("ShortrangeElectrostaticsCompoundOnly", { simparams }, envmode);

	//	env->run();

	//	auto sim = env->GetSim();

	//	//LIMA_Print::plotEnergies(env->getAnalyzedPackage()->pot_energy, env->getAnalyzedPackage()->kin_energy, env->getAnalyzedPackage()->total_energy);

	//	
	//	//ASSERT(sim->boxparams_host.boxSize == BoxGrid::blocksizeNM * 3, "This test assumes entire BoxGrid is in Shortrange range");

	//	// First check that the potential energy is calculated as we would expect if we do it the simple way
	//	float maxForceError = 0.f;
	//	for (int cidSelf = 0; cidSelf < sim->box->boxparams.n_compounds; cidSelf++) {
	//		double potESum{};
	//		Float3 forceSum{};

	//		const Compound& compoundSelf = sim->box->compounds[cidSelf];
	//		const CompoundInterimState& compoundInterimSelf = sim->box->compoundInterimStates[cidSelf];
	//		const float chargeSelf = compoundSelf.atom_charges[0];

	//		// The final Coulomb force is calculated using the position from the second-to-last step, thus -2 not -1
	//		//const Float3 posSelfAbs = sim->traj_buffer->GetMostRecentCompoundparticleDatapoint(cidSelf, 0, simparams.n_steps - 2);
	//		const Float3 posSelfAbs = sim->traj_buffer->getCompoundparticleDatapointAtIndex(cidSelf, 0, simparams.n_steps - 2);	// We MUST have the position at this index, to get accurate forces
	//		const NodeIndex nodeindexSelf = LIMAPOSITIONSYSTEM::PositionToNodeIndexNM(posSelfAbs);
	//		const Float3 posSelfRel = posSelfAbs - LIMAPOSITIONSYSTEM::nodeIndexToAbsolutePosition(nodeindexSelf);

	//		for (int cidOther = 0; cidOther < sim->box->boxparams.n_compounds; cidOther++) {
	//			if (cidSelf == cidOther)
	//				continue;
	//			
	//			const auto& compoundOther = sim->box->compounds[cidOther];
	//			const float chargeOther = compoundOther.atom_charges[0];
	//			const Float3 posOtherAbs = sim->traj_buffer->getCompoundparticleDatapointAtIndex(cidOther, 0, simparams.n_steps - 2);
	//			const Float3 posOtherRelativeToSelf = GetPositionOfParticleRelativeToSelfUsingTheWierdLogicOfTheKernel(posOtherAbs, nodeindexSelf);

	//			const Float3 diff = posSelfRel - posOtherRelativeToSelf;

	//			potESum += PhysicsUtils::CalcCoulumbPotential(chargeSelf, chargeOther, diff.len()) * 0.5f;
	//			forceSum += PhysicsUtils::CalcCoulumbForce(chargeSelf, chargeOther, diff);
	//		}
	//		
	//		// Need a expected error because in the test we do true hyperdist, but in sim we do no hyperdist
	//		// The error arises because a particle is moved 1 boxlen, not when it is correct for hyperPos, but when it moves into the next node in the boxgrid
	//		// Thus this error arises only when the box is so small that a particle go directly from nodes such as (-1, 0 0) to (1,0,0)
	//		//const float potEError = std::abs(compoundInterimSelf.sumPotentialenergy(0) - potESum) / potESum;
	//		const float forceError = std::abs((compoundInterimSelf.forces_prev[0] - forceSum).len()) / forceSum.len();
	//		maxForceError = std::max(maxForceError, forceError);

	//		//ASSERT(potEError < 1e-4, Lima::Format("Actual PotE {:.7e} Expected potE: {:.7e} Error {:.7e}", compoundSelf.potE_interim[0], potESum, potEError));
	//		//ASSERT(forceError < 1e-4, Lima::Format("Actual Force {:.7e} Expected force {:.7e} Error {:.7e}", compoundSelf.forces_interim[0].len(), forceSum.len(), forceError));
	//	}

	//	// Now do the normal VC check
	//	const float targetVarCoeff = 8.66e-3f;
	//	auto analytics = SimAnalysis::analyzeEnergy(sim.get());


	//	ASSERT(analytics.variance_coefficient < targetVarCoeff, Lima::Format("VC {:.3e} / {:.3e}", analytics.variance_coefficient, targetVarCoeff));

	//	return LimaUnittestResult{ 
	//		true, 
	//		Lima::Format("VC {:.3e} / {:.3e} Max F error {:.3e}", analytics.variance_coefficient, targetVarCoeff, maxForceError),
	//		envmode == Full };
	//}


	// The electrostatic forces and energy of random +-1e ions (no LJ) against an exact Ewald sum in double precision: the realspace sum within
	// the cutoff, and the reciprocal sum over every wavevector where it is not negligible. SPME matches it to ~0.1% RMS.
	// Many charges and an RMS error, rather than pair forces: SPME gives each charge a small force from its own spread charge
	// (~0.7 kJ/mol/nm for 1e at 0.125 nm spacing), which dominates the error of a lone pair but averages out in a system
	TestRoutine TestPmeMatchesEwaldSum(Environment& environment, EnvMode envmode) {
		const float boxLen = 5.f;
		const fs::path workDir = HeavyTestsDir() / "PmeMatchesEwaldSum";
		const AtomsSelection atoms{
			{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "lxx", "lp", 0, 1.f, 12.011f}, 50},
			{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "lxx", "ln", 0, -1.f, 12.011f}, 50},
		};

		SimParams params{};
		params.n_steps = 1;
		params.data_logging_interval = 1;
		SimulationJob job;
		job.workDir = workDir;
		job.grofile.emplace();
		job.topfile.emplace();
		job.simParams = params;
		job.mode = envmode;
		auto generatedGrofile = std::make_shared<GroFile>();
		job.preprocess = [workDir, boxLen, atoms, generatedGrofile](GroFile& grofile, TopologyFile& topfile, SimParams&) {
			fs::create_directories(workDir);
			MakeChargeParticlesSim(grofile, topfile, workDir, boxLen, atoms, 1.f);
			*generatedGrofile = grofile;
		};
		auto completed = co_await environment.Submit(std::move(job));
		const GroFile& grofile = *generatedGrofile;
		const int n = static_cast<int>(grofile.atoms.size());

		const double kappa = PhysicsUtils::CalcEwaldkappa(params.cutoff_nm);	// [nm^-1]
		const double cutoff = params.cutoff_nm;									// [nm]
		const double coulomb = PhysicsUtils::modifiedCoulombConstant;
		const double volume = double(boxLen) * boxLen * boxLen;
		std::vector<double> charges(n);
		for (int i = 0; i < n; i++)
			charges[i] = (grofile.atoms[i].atomName == "lp" ? 1. : -1.) * elementaryChargeToKiloCoulombPerMole;
		auto position = [&](int i, int dim) { return double(grofile.atoms[i].position[dim]); };
		std::vector<std::array<double, 3>> expected(n, std::array<double, 3>{});
		double expectedEnergy = 0.;

		// Realspace part, minimum image (cutoff < boxLen / 2)
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) continue;
				double diff[3];
				for (int d = 0; d < 3; d++) {
					diff[d] = position(i, d) - position(j, d);
					diff[d] -= boxLen * std::round(diff[d] / boxLen);
				}
				const double dist = std::sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
				if (dist >= cutoff) continue;
				expectedEnergy += 0.5 * coulomb * charges[i] * charges[j] * std::erfc(kappa * dist) / dist;
				const double magnitude = coulomb * charges[i] * charges[j]
					* (std::erfc(kappa * dist) / (dist * dist) + 2. * kappa / std::sqrt(PI) * std::exp(-kappa * kappa * dist * dist) / dist);
				for (int d = 0; d < 3; d++)
					expected[i][d] += magnitude * diff[d] / dist;
			}
		}

		// Reciprocal part: wavevectors up to where exp(-k^2 / 4 kappa^2) < 1e-12
		const double kMax = 2. * kappa * std::sqrt(12. * std::log(10.));
		const int nMax = static_cast<int>(std::ceil(kMax * boxLen / (2. * PI)));
		for (int nx = -nMax; nx <= nMax; nx++) {
			for (int ny = -nMax; ny <= nMax; ny++) {
				for (int nz = -nMax; nz <= nMax; nz++) {
					if (nx == 0 && ny == 0 && nz == 0) continue;
					const double k[3]{ 2. * PI * nx / boxLen, 2. * PI * ny / boxLen, 2. * PI * nz / boxLen };
					const double kSquared = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
					if (kSquared > kMax * kMax) continue;
					const double greens = 4. * PI / kSquared * std::exp(-kSquared / (4. * kappa * kappa));

					// Structure factor S(k) = sum_j q_j exp(-i k.r_j)
					double structureRe = 0., structureIm = 0.;
					for (int j = 0; j < n; j++) {
						const double phase = k[0] * position(j, 0) + k[1] * position(j, 1) + k[2] * position(j, 2);
						structureRe += charges[j] * std::cos(phase);
						structureIm -= charges[j] * std::sin(phase);
					}
					expectedEnergy += coulomb / (2. * volume) * greens * (structureRe * structureRe + structureIm * structureIm);
					for (int i = 0; i < n; i++) {
						const double phase = k[0] * position(i, 0) + k[1] * position(i, 1) + k[2] * position(i, 2);
						// Im(exp(i k.r_i) S(k)) = sum_j q_j sin(k.(r_i - r_j))
						const double im = std::sin(phase) * structureRe + std::cos(phase) * structureIm;
						for (int d = 0; d < 3; d++)
							expected[i][d] += coulomb / volume * charges[i] * greens * im * k[d];
					}
				}
			}
		}

		// Self-energy: the reciprocal sum includes each charge's interaction with its own Gaussian
		for (int i = 0; i < n; i++)
			expectedEnergy -= kappa / std::sqrt(PI) * coulomb * charges[i] * charges[i];

		double errorSquaredSum = 0., forceSquaredSum = 0., actualEnergy = 0.;
		for (int i = 0; i < n; i++) {
			const Float3 actual = completed.simulation->forceBuffer->GetDatapoint(i, 0, 0);
			actualEnergy += completed.simulation->potE_buffer->GetDatapoint(i, 0, 0);
			double errorSquared = 0.;
			for (int d = 0; d < 3; d++) {
				errorSquared += (actual[d] - expected[i][d]) * (actual[d] - expected[i][d]);
				forceSquaredSum += expected[i][d] * expected[i][d];
			}
			errorSquaredSum += errorSquared;
		}
		const double relativeRmsError = std::sqrt(errorSquaredSum / forceSquaredSum);
		const double relativeEnergyError = std::abs(actualEnergy - expectedEnergy) / std::abs(expectedEnergy);
		co_return LimaUnittestResult{ relativeRmsError < 0.005 && relativeEnergyError < 0.005,
			Lima::Format("Force err {:.2f}% Energy err {:.2f}%", relativeRmsError * 100., relativeEnergyError * 100.), envmode == Full };
	}

	LimaUnittestResult PlotPmePotAsFactorOfDistance(EnvMode envmode) {
		const fs::path work_folder = AutomatedTestsDir() / "Pool/";
		Environment& environment = Environment::Get();

		SimParams params{};
		params.n_steps = 1;
		params.data_logging_interval = 1;

		const float c0 = -1.f * elementaryChargeToKiloCoulombPerMole;
		const float c1 = 1.f * elementaryChargeToKiloCoulombPerMole;

		std::vector<Float3> actualForce, expectedForce;
		std::vector<float> actualPot, expectedPot;
		std::vector<float> distances;

		for (float dist = 15.f; dist > 0.4; dist -= 1.13f) {
			const Float3 p0{ 15.025f, 10.025f, 10.025f };
			const Float3 p1 = p0 - Float3{ dist, 0.f, 0.f };

			SimulationJob job;
			job.workDir = work_folder;
			job.simParams = params;
			job.mode = envmode;
			job.preprocess = [p0, p1](GroFile& grofile, TopologyFile&, SimParams&) {
				grofile.box_size = Float3{ 30.f };
				grofile.atoms[0].position = p0;
				grofile.atoms[1].position = p1;
			};
			job.configureSimulation = [c0, c1](Simulation& simulation) {
				simulation.box->persistentClusters[0].pqd[0].params.charge = c0;
				simulation.box->persistentClusters[1].pqd[0].params.charge = c1;
			};


			Float3 hyperposOther = p1;
			BoundaryConditionPublic::applyHyperposNM(p0, hyperposOther, Float3{ 30.f }, PBC);
			const Float3 diff = p0 - hyperposOther;


			//ForceEnergy mirrorFE = CalcImmediateMirrorForceEnergy(diff, c0 * c1, grofile.box_size);
			Float3 mirrorForce = PhysicsUtils::CalcCoulumbForce(c0, c1, Float3{ diff.x - 30.f, 0.f, 0.f });
			float mirrorPotential = PhysicsUtils::CalcCoulumbPotential(c0, c1, (Float3{ diff.x - 30.f, 0.f, 0.f }).len()) * 0.5f;
			Float3 force = PhysicsUtils::CalcCoulumbForce(c0, c1, diff);
			float pot = PhysicsUtils::CalcCoulumbPotential(c0, c1, diff.len()) * 0.5f;

			expectedPot.push_back(pot);
			expectedForce.push_back(force + mirrorForce);

			auto sim = environment.Submit(std::move(job)).Get().simulation;

			actualPot.push_back(sim->potE_buffer->GetDatapoint(0, 0, 0));
			actualForce.push_back(sim->forceBuffer->GetDatapoint(0, 0, 0));

			distances.push_back(dist);
		}

		//PlotUtils::PlotData({ actual, expected }, {"Actual Pot", "Expected Pot"}, distances);
		std::vector<float> potError(actualPot.size());
		std::vector<float> forceError(actualPot.size());
		for (int i = 0; i < potError.size(); i++) {
			potError[i] = std::abs((actualPot[i] - expectedPot[i]) / expectedPot[i]);
			forceError[i] = std::abs((actualForce[i] - expectedForce[i]).len() / expectedForce[i].len());
		}

		//PlotUtils::PlotData({ actual, expected, errors }, { "Actual Force mag", "Expected Force mag", "Error"}, distances);
		PlotUtils::PlotData({ potError, forceError}, { "PotE Error", "Force Error"}, distances);
		return LimaUnittestResult{ true, "", envmode == Full };
	}

	LimaUnittestResult TestConsistentEnergyWhenGoingFromLresToSres(EnvMode envmode) {
		const fs::path work_folder = AutomatedTestsDir() / "Pool/";
		Environment& environment = Environment::Get();


		SimParams params{};
		params.n_steps = 500;
		params.data_logging_interval = 1;
		params.dt = 1.f * FEMTO_TO_NANO;

		const Float3 p0{ 7.0f, 10.f, 10.f };
		const Float3 p1{ 9.f, 10.f, 10.f };

		const float c0 = 1.f * elementaryChargeToKiloCoulombPerMole;
		const float c1 = -c0;

		SimulationJob job;
		job.workDir = work_folder;
		job.simParams = params;
		job.mode = envmode;
		job.postprocess = SimAnalysis::AnalyzeEnergy;
		job.preprocess = [p0, p1](GroFile& grofile, TopologyFile&, SimParams&) {
			grofile.box_size = Float3{ 20.f };
			grofile.atoms[0].position = p0;
			grofile.atoms[1].position = p1;
		};

		Float3 hyperposOther = p1;
		BoundaryConditionPublic::applyHyperposNM(p0, hyperposOther, Float3{ 20.f }, PBC);
		const Float3 diff = p0 - hyperposOther;
		const float expectedPotential = PhysicsUtils::CalcCoulumbPotential(c0, c1, diff.len()) * 0.5f;
		const Float3 expectedForce = PhysicsUtils::CalcCoulumbForce(c0, c1, diff);



		auto result = environment.Submit(std::move(job)).Get();

		// First go trough the traj data and find the step where the particles are less than 0.5 nm apart
		int step = -1;
		for (int i = 0; i < params.n_steps; i++)
		{
			const Float3 pos0 = result.simulation->traj_buffer->GetDatapoint(0, 0, i);
			const Float3 pos1 = result.simulation->traj_buffer->GetDatapoint(1, 0, i);
			const float dist = (pos0 - pos1).len();
			if (dist < 0.4f) {
				step = i;
				break;
			}
		}

		// Only use the energies up untill the step found above
		const auto& anal = *result.analysis;
		//LIMA_Print::plotEnergies(std::span(anal.pot_energy).subspan(0, step), std::span(anal.kin_energy).subspan(0, step), std::span(anal.total_energy).subspan(0, step));

		return LimaUnittestResult{ anal.variance_coefficient < 1e-3f, "", envmode == Full };
	}



	// Create many pos charged Ions as compounds. Set all LJ to 0. Compute exact SR and LR interactions between all particles. Run simulation 1 step, and compare the errors
	TestRoutine TestLongrangeEsNoLJManyParticles(
		Environment& environment, EnvMode envmode) {
		const Float3 boxlen{ 20.f };
		const float chargeExtern = 1.f;
		const float charge = chargeExtern * elementaryChargeToKiloCoulombPerMole;
		AtomsSelection atoms{
				{TopologyFile::AtomsEntry{";residue_X", 0, "lt1", 0, "lxx", "lxx", 0, chargeExtern, 12.011f}, 100}, // by naming the residue lxx we let these particles be full molecules, instead of tinymols, is that ideal? Does it matter? 
			};

		const fs::path work_folder = HeavyTestsDir() / "ShortrangeElectrostaticsCompoundOnly/";

		
		
		SimParams params{};
		params.n_steps = 2;
		params.data_logging_interval = 1;
		SimulationJob job;
		job.workDir = work_folder;
		job.grofile.emplace();
		job.topfile.emplace();
		job.simParams = params;
		job.mode = envmode;
		auto generatedGrofile = std::make_shared<GroFile>();
		job.preprocess = [work_folder, boxlen, atoms = std::move(atoms), generatedGrofile](
			GroFile& grofile, TopologyFile& topfile, SimParams&) {
			MakeChargeParticlesSim(grofile, topfile, work_folder, boxlen.x, atoms, 1.f);
			*generatedGrofile = grofile;
		};
		auto result = co_await environment.Submit(std::move(job));
		const GroFile& grofile = *generatedGrofile;


		// Now compute all expected forces and potentials
		std::vector<float> expectedPotentials(grofile.atoms.size());
		std::vector<Float3> expectedForces(grofile.atoms.size());

		for (int i = 0; i < grofile.atoms.size(); i++) {
			const auto& atom = grofile.atoms[i];
			float potential{};
			Float3 force{};
			for (int j = 0; j < grofile.atoms.size(); j++) {
				if (i == j)
					continue;

				Float3 hyperposOther = grofile.atoms[j].position;
				BoundaryConditionPublic::applyHyperposNM(atom.position, hyperposOther, boxlen, PBC);
				const Float3 diff = atom.position - hyperposOther;
				const float expectedPotential = PhysicsUtils::CalcCoulumbPotential(charge, charge, diff.len()) * 0.5f;
				const Float3 expectedForce = PhysicsUtils::CalcCoulumbForce(charge, charge, diff);

				potential += expectedPotential;
				force += expectedForce;
			}
			expectedPotentials[i] = potential;
			expectedForces[i] = force;
		}

		const auto& sim = result.simulation;
		const Float3 actualForce = sim->forceBuffer->GetDatapoint(0, 0, 0);

		std::vector<float> potErrors(grofile.atoms.size());
		std::vector<float> forceErrors(grofile.atoms.size());

		for (int i = 0; i < grofile.atoms.size(); i++) {
			const float potEError = std::abs(sim->potE_buffer->GetDatapoint(i, 0, 0) - expectedPotentials[i]) / expectedPotentials[i];
			Float3 actualForce = sim->forceBuffer->GetDatapoint(i, 0, 0);
			Float3 expectedForce = expectedForces[i];
			Float3 position = grofile.atoms[i].position;
			const float forceError = (sim->forceBuffer->GetDatapoint(i, 0, 0) - expectedForces[i]).len() / expectedForces[i].len();

			if (expectedForces[i].len() < 50'000.f) // [J/mol/nm
				continue; // Force is quite small, hard to be relative accurate here

			if (forceError > 4.f)
				int a = 0;

			potErrors[i] = potEError;
			forceErrors[i] = forceError;			
		}

		/*const float maxPotError = *std::max_element(potErrors.begin(), potErrors.end());
		const float meanPotError = Statistics::Mean(potErrors);
		ASSERT(meanPotError < 5e-2, Lima::Format("Mean PotE Error {:.3e}", meanPotError));
		ASSERT(maxPotError < 1, Lima::Format("Max PotE Error {:.3e}", maxPotError));*/
		
		const float maxForceError = *std::max_element(forceErrors.begin(), forceErrors.end());
		const float meanForceError = Statistics::Mean(forceErrors);
		if (meanForceError >= 0.18f)
			co_return LimaUnittestResult{ false, Lima::Format("Mean Force Error {:.3f}", meanForceError), envmode == Full };
		//ASSERT(maxForceError < 0.8f, Lima::Format("Max Force Error {:.3e}", maxForceError));

		co_return LimaUnittestResult{ true, "", envmode == Full };
	}
}


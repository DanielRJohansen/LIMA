#include "EnergyMinimizationTypes.h"
#include "Bodies.cuh"
#include "Simulation.cuh"
#include "ParallelFor.h"

#include <algorithm>
#include <cmath>

namespace EM {
	namespace {
		Float3 MinimumImage(Float3 difference, Float3 boxSize) {
			difference.x -= boxSize.x * std::round(difference.x / boxSize.x);
			difference.y -= boxSize.y * std::round(difference.y / boxSize.y);
			difference.z -= boxSize.z * std::round(difference.z / boxSize.z);
			return difference;
		}
	}

	std::vector<float> ComputeInverseStiffness(const std::vector<PersistentCluster>& pclusters, const std::vector<PersistentClusterMeta>& pclusterMeta,
		const BondGroups& bonds, Float3 boxSize, float nonbondedStiffness) {
		constexpr float minDistance = 0.08f; // [nm] Guards against overlapping atoms in unminimized input
		std::vector<double> stiffness(pclusters.size() * PersistentCluster::maxParticles, 0.);

		for (const BondGroup& group : bonds.groups) {
			const auto Slot = [&](int idInGroup) {
				const BondGroup::ParticleRef ref = bonds.particles[group.indexOfFirstParticle + idInGroup];
				return ref.pcid * PersistentCluster::maxParticles + ref.pid;
			};
			const auto Position = [&](int idInGroup) {
				const BondGroup::ParticleRef ref = bonds.particles[group.indexOfFirstParticle + idInGroup];
				return pclusters[ref.pcid].pqd[ref.pid].position;
			};
			const auto Distance = [&](int a, int b) {
				return std::max(MinimumImage(Position(a) - Position(b), boxSize).len(), minDistance);
			};
			// Distance from an outer atom to the axis of a torsion
			const auto AxisDistance = [&](int outer, int axis0, int axis1) {
				const Float3 axis = MinimumImage(Position(axis1) - Position(axis0), boxSize).norm();
				const Float3 arm = MinimumImage(Position(outer) - Position(axis0), boxSize);
				return std::max((arm - axis * arm.dot(axis)).len(), minDistance);
			};
			const auto AddTorsion = [&](const uint8_t* ids, double k) {
				const double outer0 = k / std::pow(AxisDistance(ids[0], ids[1], ids[2]), 2);
				const double outer1 = k / std::pow(AxisDistance(ids[3], ids[2], ids[1]), 2);
				stiffness[Slot(ids[0])] += outer0;
				stiffness[Slot(ids[1])] += outer0;
				stiffness[Slot(ids[2])] += outer1;
				stiffness[Slot(ids[3])] += outer1;
			};

			for (int i = 0; i < group.nSinglebonds; i++) {
				const SingleBond& bond = bonds.singlebonds[group.indexOfFirstSinglebond + i];
				stiffness[Slot(bond.idInBondgroup[0])] += bond.params.kb;
				stiffness[Slot(bond.idInBondgroup[1])] += bond.params.kb;
			}
			for (int i = 0; i < group.nAnglebonds; i++) {
				const AngleUreyBradleyBond& angle = bonds.anglebonds[group.indexOfFirstAnglebond + i];
				const uint8_t* ids = angle.atom_indexes;
				const double left = Distance(ids[0], ids[1]);
				const double right = Distance(ids[2], ids[1]);
				stiffness[Slot(ids[0])] += angle.params.kTheta / (left * left) + angle.params.kUB;
				stiffness[Slot(ids[2])] += angle.params.kTheta / (right * right) + angle.params.kUB;
				stiffness[Slot(ids[1])] += angle.params.kTheta * std::pow(1. / left + 1. / right, 2);
			}
			for (int i = 0; i < group.nDihedralbonds; i++) {
				const DihedralBond& dihedral = bonds.dihedralbonds[group.indexOfFirstDihedralbond + i];
				AddTorsion(dihedral.atom_indexes, std::abs(dihedral.params.k_phi) * dihedral.params.n * dihedral.params.n);
			}
			for (int i = 0; i < group.nImproperdihedralbonds; i++) {
				const ImproperDihedralBond& improper = bonds.improperdihedralbonds[group.indexOfFirstImproperdihedralbond + i];
				AddTorsion(improper.atom_indexes, improper.params.k_psi);
			}
		}

		std::vector<float> inverseStiffness(stiffness.size(), 0.f);
		for (size_t pc = 0; pc < pclusterMeta.size(); pc++)
			for (int lane = 0; lane < PersistentCluster::maxParticles; lane++)
				if (pclusterMeta[pc].particleIdsGlobal[lane] != -1) {
					const size_t slot = pc * PersistentCluster::maxParticles + lane;
					inverseStiffness[slot] = static_cast<float>(1. / (stiffness[slot] + nonbondedStiffness));
				}
		return inverseStiffness;
	}

	std::vector<uint8_t> FindWholeMoleculePclusters(size_t nPclusters, const BondGroups& bonds) {
		std::vector<uint8_t> wholeMolecule(nPclusters, 1);
		for (const BondGroup& group : bonds.groups) {
			const auto Pcluster = [&](int idInGroup) { return bonds.particles[group.indexOfFirstParticle + idInGroup].pcid; };
			const auto Mark = [&](const uint8_t* ids, int nIds) {
				for (int i = 1; i < nIds; i++)
					if (Pcluster(ids[i]) != Pcluster(ids[0])) {
						for (int j = 0; j < nIds; j++) wholeMolecule[Pcluster(ids[j])] = 0;
						return;
					}
			};
			for (int i = 0; i < group.nSinglebonds; i++) Mark(bonds.singlebonds[group.indexOfFirstSinglebond + i].idInBondgroup, 2);
			for (int i = 0; i < group.nAnglebonds; i++) Mark(bonds.anglebonds[group.indexOfFirstAnglebond + i].atom_indexes, 3);
			for (int i = 0; i < group.nDihedralbonds; i++) Mark(bonds.dihedralbonds[group.indexOfFirstDihedralbond + i].atom_indexes, 4);
			for (int i = 0; i < group.nImproperdihedralbonds; i++) Mark(bonds.improperdihedralbonds[group.indexOfFirstImproperdihedralbond + i].atom_indexes, 4);
		}
		return wholeMolecule;
	}

	bool Preconditioner::Fits(const Box& box) const {
		return inverseStiffness.size() == box.persistentClusters.size() * PersistentCluster::maxParticles
			&& wholeMolecule.size() == box.persistentClusters.size();
	}

	Preconditioner MakePreconditioner(const Box& box, float nonbondedStiffness) {
		return {
			ComputeInverseStiffness(box.persistentClusters, box.persistentClustersMetadata, box.bondgroups,
				box.boxparams.BoxSizeFloat(), nonbondedStiffness),
			FindWholeMoleculePclusters(box.persistentClusters.size(), box.bondgroups)
		};
	}

	void MakeMissingPreconditioners(const std::vector<const Box*>& boxes, const std::vector<Preconditioner*>& preconditioners,
		float nonbondedStiffness) {
		if (boxes.size() != preconditioners.size())
			throw std::invalid_argument("Preconditioner count does not match box count");
		ParallelUtils::ParallelFor(boxes.size(), [&](size_t i) {
			if (!preconditioners[i]->Fits(*boxes[i]))
				*preconditioners[i] = MakePreconditioner(*boxes[i], nonbondedStiffness);
			});
	}

}

#include "LIMA_ENGINE/include/EngineUtils.cuh"






namespace DEBUGUTILS {

	void DEBUGUTILS::findAllNearestSolventSolvent(SolventBlockGrid* solventblockgrid, size_t n_solvents, std::vector<float>& closestNeighbor)
	{
		printf("Finding nearest neighbor-solvent of each solvent..");
		closestNeighbor.resize(n_solvents);
		std::fill(closestNeighbor.begin(), closestNeighbor.end(), FLT_MAX);

		int solvent_index = 0;

		for (int sbi = 0; sbi < SolventBlockGrid::blocks_total; sbi++) {
			auto sb = solventblockgrid->getBlockPtr(sbi);
			for (int i = 0; i < sb->n_solvents; i++) {
				auto posi = sb->rel_pos[i];

				// Loop through all solvents at equal or greater index
				for (int sbj = 0; sbj < SolventBlockGrid::blocks_total; sbj++) {
					auto sb2 = solventblockgrid->getBlockPtr(sbj);
					for (int j = 0; j < sb2->n_solvents; j++) {

						if (sbi == sbj && i == j) { continue; }	// same solvent

						auto posj = sb2->rel_pos[j];
						auto dist = EngineUtils::calcDistance(sb->origo, posi, sb2->origo, posj);

						closestNeighbor[solvent_index] = std::min(closestNeighbor[solvent_index], dist);

					}
				}

				solvent_index++;
			}
		}
		printf(" Done!\n");
	}



}

#include "EngineUtils.cuh"






namespace DEBUGUTILS {
/*
	void DEBUGUTILS::findAllNearestSolventSolvent(SolventBlocksCircularQueue* queue, size_t n_solvents, std::vector<float>& closestNeighbor)
	{
		printf("Finding nearest neighbor-solvent of each solvent..");
		closestNeighbor.resize(n_solvents);
		std::fill(closestNeighbor.begin(), closestNeighbor.end(), FLT_MAX);

		int solvent_index = 0;

		for (int sbi = 0; sbi < queue->blocks_per_grid; sbi++) {
			auto sb = queue->getBlockPtr(sbi, 0);		// Look at step 0
			for (int i = 0; i < sb->n_solvents; i++) {
				auto posi = sb->rel_pos[i];

				// Loop through all solvents at equal or greater index
				for (int sbj = 0; sbj < queue->blocks_per_grid; sbj++) {
					auto sb2 = queue->getBlockPtr(sbj, 0);
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


*/
}

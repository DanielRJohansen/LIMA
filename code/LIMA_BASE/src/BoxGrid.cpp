#include "BoxGrid.cuh"
#include "BoundaryConditionPublic.h"


BoxGrid::TinymolBlockAdjacency::BlockRef* BoxGrid::TinymolBlockAdjacency::PrecomputeNeabyBlockIds(Int3 boxlenNM, float ljCutoffNm) {

    const float maxAllowedCutoff = sqrt(1*1+1*1+0);
    //if (ljCutoffNm > maxAllowedCutoff)
    //    throw std::invalid_argument(std::format("Due to harcoded optimizations of neighborlists, LIMA does not allow a cutoff above {}, was {})", maxAllowedCutoff, ljCutoffNm);

	const int blocksTotal = BoxGrid::BlocksTotal(boxlenNM);

	std::vector<BlockRef> nearbyBlockIds(blocksTotal * nNearbyBlocks);
	int globalIndex = 0;

	for (int i = 0; i < blocksTotal; i++) {
		NodeIndex index3d = Get3dIndex(i, boxlenNM);

		const int query_range = 2;
		for (int z = -query_range; z <= query_range; z++) {
			for (int y = -query_range; y <= query_range; y++) {
				for (int x = -query_range; x <= query_range; x++) {
					const NodeIndex dir{ x,y,z };

                    if (dir.largestMagnitudeElement() == 1 || (dir.largestMagnitudeElement() == 2 && dir.Magnitude() == 2)) {

                        NodeIndex nearbyIndex = NodeIndex{ index3d.x + x, index3d.y + y, index3d.z + z };
                        BoundaryConditionPublic::applyBC(nearbyIndex, boxlenNM);
						assert(Get1dIndex(nearbyIndex, boxlenNM) < UINT16_MAX);

						Float3 a = NodeIndex{ x,y,z }.toFloat3();
						assert(a.x >= -2 && a.x <= 2);
						assert(a.y >= -2 && a.y <= 2);
						assert(a.z >= -2 && a.z <= 2);
                        nearbyBlockIds[globalIndex++] = BlockRef{ static_cast<uint16_t>(Get1dIndex(nearbyIndex, boxlenNM)), Float3Compressed(NodeIndex{x,y,z}.toFloat3())};
                    }
				}
			}
		}
		assert(globalIndex % nNearbyBlocks == 0);
	}
	assert(globalIndex == nearbyBlockIds.size());
	return GenericCopyToDevice(nearbyBlockIds);
}


BoxGrid::TinymolBlockAdjacency::NearbyBlocksSequences* BoxGrid::TinymolBlockAdjacency::PrecomputeNearbyBlockSequences(Int3 boxlenNM) {
	const int blocksTotal = BoxGrid::BlocksTotal(boxlenNM);
	std::vector<NearbyBlocksSequences> nearbyBlockIds(blocksTotal);
	

	auto InsertSequence = [](NearbyBlocksSequences& nearbySequence, std::vector<int>& sequence) {
		if (nearbySequence.nSequences >= NearbyBlocksSequences::maxSequences)
			throw std::runtime_error("Exceeded max nearby block sequences");

		//printf("Sequence:");
		//for (const int bidx : sequence) {
		//	printf("%d ", bidx);
		//}
		//printf("\n");

		nearbySequence.sequences[nearbySequence.nSequences].blockIndexStart = sequence.front();
		nearbySequence.sequences[nearbySequence.nSequences].nBlocks = static_cast<int>(sequence.size());
		nearbySequence.nSequences++;
		sequence.clear();
		};


	for (int i = 0; i < blocksTotal; i++) {
		NodeIndex index3d = Get3dIndex(i, boxlenNM);

		const int query_range = 2;
		for (int z = -query_range; z <= query_range; z++) {
			for (int y = -query_range; y <= query_range; y++) {

				std::vector<int> sequence;

				for (int x = -query_range; x <= query_range; x++) {
					const NodeIndex dir{ x,y,z };


					// Blocks mustn't include self in sequence
					if (dir.largestMagnitudeElement() == 0) {
						if (sequence.empty())
							throw std::runtime_error("How is this even possible?");
						InsertSequence(nearbyBlockIds[i], sequence);
					}



					if (!(dir.largestMagnitudeElement() == 1 || (dir.largestMagnitudeElement() == 2 && dir.Magnitude() == 2)))
						continue;

					NodeIndex nearbyIndex = NodeIndex{ index3d.x + x, index3d.y + y, index3d.z + z };
					BoundaryConditionPublic::applyBC(nearbyIndex, boxlenNM);
					assert(Get1dIndex(nearbyIndex, boxlenNM) < UINT16_MAX);

					int globalIndex = Get1dIndex(nearbyIndex, boxlenNM);

					if (sequence.empty())
						sequence.push_back(globalIndex);
					else if (sequence.back() + 1 == globalIndex)
						sequence.push_back(globalIndex);
					else {
						// Store sequence
						InsertSequence(nearbyBlockIds[i], sequence);
						sequence.push_back(globalIndex);
					}
				}

				// Store last sequence
				if (!sequence.empty()) {
					InsertSequence(nearbyBlockIds[i], sequence);
				}
			}
		}
	}

	return GenericCopyToDevice(nearbyBlockIds);
}
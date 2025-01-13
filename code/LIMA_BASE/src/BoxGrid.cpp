#include "BoxGrid.cuh"
#include "BoundaryConditionPublic.h"


BoxGrid::TinymolBlockAdjacency::BlockRef* BoxGrid::TinymolBlockAdjacency::PrecomputeNeabyBlockIds(int boxlenNM, float ljCutoffNm) {

    const float maxAllowedCutoff = sqrt(1*1+1*1+0);
    //if (ljCutoffNm > maxAllowedCutoff)
    //    throw std::invalid_argument(std::format("Due to harcoded optimizations of neighborlists, LIMA does not allow a cutoff above {}, was {})", maxAllowedCutoff, ljCutoffNm);

	const int blocksPerDim = NodesPerDim(boxlenNM);
	const int blocksTotal = BlocksTotal(blocksPerDim);

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
                        nearbyBlockIds[globalIndex++] = BlockRef{ Get1dIndex(nearbyIndex, boxlenNM), NodeIndex{x,y,z}.toFloat3()};
                    }
				}
			}
		}
		assert(globalIndex % nNearbyBlocks == 0);
	}
	assert(globalIndex == nearbyBlockIds.size());
	return GenericCopyToDevice(nearbyBlockIds);
}



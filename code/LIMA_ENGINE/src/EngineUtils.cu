#include "EngineUtils.cuh"



__host__ void CompoundGridNode::addCompound(int16_t compound_id)
{
	if (n_nearby_compounds == max_elements) {
		throw std::exception("Failed to add compound to CompoundGridNode\n");
	}
	nearby_compound_ids[n_nearby_compounds++] = compound_id;
}

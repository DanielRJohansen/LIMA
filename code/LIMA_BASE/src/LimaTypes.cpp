#include "LimaTypes.cuh"
#include "RenderUtilities.cuh"
#include <limits>

RenderAtom::RenderAtom(Float3 pos, Float3 boxSize, char atomLetter) {
	const float radius = RenderUtilities::getRadius(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
	position = float4{ pos.x, pos.y, pos.z, radius };
	color = RenderUtilities::getColor(RenderUtilities::RAS_getTypeFromAtomletter(atomLetter));
}

BoundingBox::BoundingBox(const std::vector<Float3>& points) {
	min = Float3{std::numeric_limits<float>::max()};
	max = Float3{std::numeric_limits<float>::min()};
	for (const Float3& p : points) {
		for (int dim = 0; dim < 3; dim++) {
			min[dim] = std::min(min[dim], p[dim]);
			max[dim] = std::max(max[dim], p[dim]);
		}
	}
}


__host__ void BondedParticlesLUT::printMatrix(int n) const {
	// Print column indices
	std::cout << "     ";  // Space for row indices
	for (int j = 0; j < n; ++j) {
		std::string separator = j + 1 < 10 ? "  " : " ";
		std::cout << j + 1 << separator;
	}
	std::cout << '\n';

	// Print separator
	std::cout << "   +";
	for (int j = 0; j < n; ++j) {
		std::cout << "---";
	}
	std::cout << '\n';

	// Print rows with row indices
	for (int i = 0; i < n; ++i) {
		std::string separator = i + 1 < 10 ? "  | " : " | ";
		std::cout << i + 1 << separator;  // Row index
		for (int j = 0; j < n; ++j) {
			std::cout << (get(i, j) ? 'X' : 'O') << "  ";
		}
		std::cout << '\n';
	}
}
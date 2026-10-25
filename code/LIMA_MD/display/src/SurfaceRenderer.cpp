#include "SurfaceRenderer.h"

#include "Shaders.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace {

constexpr float isovalue = 0.5f;
// Each atom's density is exp(-k d^2 / r^2), with k chosen so a lone atom's isosurface lies at d = r
const float falloff = std::log(1.f / isovalue);
// Beyond this many radii an atom's density is below 1% of the isovalue, and is left out
constexpr float cutoffRadii = 2.6f;
// Larger systems get a coarser grid, to bound memory and time
constexpr std::size_t maxGridNodes = 48'000'000;

struct Grid {
	glm::vec3 origin{};
	float spacing = 0.f;
	glm::ivec3 size{};
	std::vector<float> density;
	std::vector<float> strongest;	// The largest single-atom density at each node
	std::vector<int> strongestAtom;	// Index into the atoms
	// Atom colors weighted by density^3: blends between atoms, but keeps the boundaries fairly sharp
	std::vector<glm::vec4> colorSum;	// {weighted rgb, total weight}

	std::size_t Index(int x, int y, int z) const {
		return static_cast<std::size_t>(x) + static_cast<std::size_t>(size.x) * (y + static_cast<std::size_t>(size.y) * z);
	}
};

Grid BuildDensity(const std::vector<SurfaceRenderer::Atom>& atoms, const std::vector<Float3>& positions)
{
	glm::vec3 low(std::numeric_limits<float>::max());
	glm::vec3 high(std::numeric_limits<float>::lowest());
	for (const auto& atom : atoms) {
		const glm::vec3 p{ positions[atom.id].x, positions[atom.id].y, positions[atom.id].z };
		low = glm::min(low, p - atom.radius * cutoffRadii);
		high = glm::max(high, p + atom.radius * cutoffRadii);
	}

	Grid grid;
	grid.spacing = SurfaceRenderer::gridSpacing;
	const auto SizeFor = [&](float spacing) { return glm::ivec3((high - low) / spacing) + 3; };
	while (static_cast<std::size_t>(SizeFor(grid.spacing).x) * SizeFor(grid.spacing).y * SizeFor(grid.spacing).z > maxGridNodes)
		grid.spacing *= 1.25f;
	grid.size = SizeFor(grid.spacing);
	grid.origin = low - grid.spacing;
	const std::size_t nodes = static_cast<std::size_t>(grid.size.x) * grid.size.y * grid.size.z;
	grid.density.assign(nodes, 0.f);
	grid.strongest.assign(nodes, 0.f);
	grid.strongestAtom.assign(nodes, -1);
	grid.colorSum.assign(nodes, glm::vec4(0.f));

	for (int a = 0; a < static_cast<int>(atoms.size()); ++a) {
		const auto& atom = atoms[a];
		const glm::vec3 p{ positions[atom.id].x, positions[atom.id].y, positions[atom.id].z };
		const glm::vec3 center = (p - grid.origin) / grid.spacing;
		const float reach = atom.radius * cutoffRadii / grid.spacing;
		const glm::ivec3 from = glm::max(glm::ivec3(glm::floor(center - reach)), glm::ivec3(0));
		const glm::ivec3 to = glm::min(glm::ivec3(glm::ceil(center + reach)), grid.size - 1);
		const float k = falloff / (atom.radius * atom.radius) * grid.spacing * grid.spacing;	// In grid units
		for (int z = from.z; z <= to.z; ++z) {
			for (int y = from.y; y <= to.y; ++y) {
				for (int x = from.x; x <= to.x; ++x) {
					const glm::vec3 d = glm::vec3(x, y, z) - center;
					const float d2 = glm::dot(d, d);
					if (d2 > reach * reach)
						continue;
					const float value = std::exp(-k * d2);
					const std::size_t i = grid.Index(x, y, z);
					grid.density[i] += value;
					const float weight = value * value * value;
					grid.colorSum[i] += glm::vec4(glm::vec3(atom.color) * weight, weight);
					if (value > grid.strongest[i]) {
						grid.strongest[i] = value;
						grid.strongestAtom[i] = a;
					}
				}
			}
		}
	}
	return grid;
}

// The 12 edges of a cell, as pairs of corners. Corner c is at offset (c & 1, c >> 1 & 1, c >> 2 & 1)
constexpr std::array<std::array<int, 2>, 12> cellEdges{ {
	{0, 1}, {2, 3}, {4, 5}, {6, 7},	// Along x
	{0, 2}, {1, 3}, {4, 6}, {5, 7},	// Along y
	{0, 4}, {1, 5}, {2, 6}, {3, 7},	// Along z
} };

glm::ivec3 Corner(int c) { return { c & 1, c >> 1 & 1, c >> 2 & 1 }; }

} // namespace

// Surface nets: one vertex per cell the isosurface passes through, at the mean of where it crosses the cell's edges,
// and one quad per grid edge it crosses, joining the four cells around that edge
SurfaceRenderer::Mesh SurfaceRenderer::BuildMesh(const std::vector<Atom>& atoms, const std::vector<Float3>& positions)
{
	Mesh mesh;
	if (atoms.empty())
		return mesh;
	const Grid grid = BuildDensity(atoms, positions);
	const glm::ivec3 cells = grid.size - 1;
	std::vector<int> cellVertex(static_cast<std::size_t>(cells.x) * cells.y * cells.z, -1);
	const auto CellIndex = [&cells](int x, int y, int z) {
		return static_cast<std::size_t>(x) + static_cast<std::size_t>(cells.x) * (y + static_cast<std::size_t>(cells.y) * z);
	};

	for (int z = 0; z < cells.z; ++z) {
		for (int y = 0; y < cells.y; ++y) {
			for (int x = 0; x < cells.x; ++x) {
				std::array<float, 8> values{};
				int inside = 0;
				for (int c = 0; c < 8; ++c) {
					const glm::ivec3 o = Corner(c);
					values[c] = grid.density[grid.Index(x + o.x, y + o.y, z + o.z)];
					inside += values[c] > isovalue;
				}
				if (inside == 0 || inside == 8)
					continue;

				glm::vec3 sum{};
				int crossings = 0;
				for (const auto& [c0, c1] : cellEdges) {
					if ((values[c0] > isovalue) == (values[c1] > isovalue))
						continue;
					const float t = (isovalue - values[c0]) / (values[c1] - values[c0]);
					sum += glm::mix(glm::vec3(Corner(c0)), glm::vec3(Corner(c1)), t);
					++crossings;
				}
				const glm::vec3 local = sum / static_cast<float>(crossings);

				// The density falls outwards, so the outward normal is against its gradient
				const glm::vec3 gradient{
					(values[1] - values[0]) + (values[3] - values[2]) + (values[5] - values[4]) + (values[7] - values[6]),
					(values[2] - values[0]) + (values[3] - values[1]) + (values[6] - values[4]) + (values[7] - values[5]),
					(values[4] - values[0]) + (values[5] - values[1]) + (values[6] - values[2]) + (values[7] - values[3]),
				};
				// Color interpolated from the corners, so color boundaries are smooth rather than following the grid.
				// Picking uses the atom with the most density at the nearest corner
				glm::vec4 colorSum{};
				for (int c = 0; c < 8; ++c) {
					const glm::ivec3 o = Corner(c);
					const glm::vec3 w = glm::mix(glm::vec3(1.f) - local, local, glm::vec3(o));
					colorSum += grid.colorSum[grid.Index(x + o.x, y + o.y, z + o.z)] * (w.x * w.y * w.z);
				}
				const glm::ivec3 nearest = glm::ivec3(x, y, z) + glm::ivec3(glm::round(local));
				const int atom = grid.strongestAtom[grid.Index(nearest.x, nearest.y, nearest.z)];
				const Atom& source = atoms[std::max(atom, 0)];
				const glm::vec4 color = colorSum.w > 0.f ? glm::vec4(glm::vec3(colorSum) / colorSum.w, 1.f) : source.color;

				cellVertex[CellIndex(x, y, z)] = static_cast<int>(mesh.vertices.size());
				mesh.vertices.push_back({ grid.origin + (glm::vec3(x, y, z) + local) * grid.spacing,
					glm::dot(gradient, gradient) > 0.f ? -glm::normalize(gradient) : glm::vec3(0.f, 0.f, 1.f),
					color, source.id });
			}
		}
	}

	// A quad for each grid edge crossing the surface, between the four cells sharing that edge
	for (int z = 1; z < cells.z; ++z) {
		for (int y = 1; y < cells.y; ++y) {
			for (int x = 1; x < cells.x; ++x) {
				const bool here = grid.density[grid.Index(x, y, z)] > isovalue;
				for (int axis = 0; axis < 3; ++axis) {
					glm::ivec3 next(x, y, z);
					next[axis] += 1;
					if (here == (grid.density[grid.Index(next.x, next.y, next.z)] > isovalue))
						continue;
					// The other two axes span the quad
					const int u = (axis + 1) % 3;
					const int v = (axis + 2) % 3;
					std::array<int, 4> quad{};
					bool complete = true;
					for (int corner = 0; corner < 4; ++corner) {
						glm::ivec3 cell(x, y, z);
						cell[u] -= (corner == 0 || corner == 3) ? 1 : 0;
						cell[v] -= (corner == 0 || corner == 1) ? 1 : 0;
						quad[corner] = cellVertex[CellIndex(cell.x, cell.y, cell.z)];
						complete &= quad[corner] >= 0;
					}
					if (!complete)
						continue;
					// Wind the quad to face outwards
					if (!here)
						std::swap(quad[1], quad[3]);
					for (const int i : { quad[0], quad[1], quad[2], quad[0], quad[2], quad[3] })
						mesh.indices.push_back(static_cast<unsigned int>(i));
				}
			}
		}
	}
	return mesh;
}

SurfaceRenderer::SurfaceRenderer()
{
	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vertexBuffer);
	glGenBuffers(1, &indexBuffer);
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, position)));
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, normal)));
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, color)));
	glEnableVertexAttribArray(3);
	glVertexAttribIPointer(3, 1, GL_INT, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, atomId)));
	glBindVertexArray(0);
}

SurfaceRenderer::~SurfaceRenderer()
{
	glDeleteBuffers(1, &indexBuffer);
	glDeleteBuffers(1, &vertexBuffer);
	glDeleteVertexArrays(1, &vao);
}

void SurfaceRenderer::Prepare(std::vector<Atom> newAtoms, const std::vector<Float3>& positions)
{
	atoms = std::move(newAtoms);
	Rebuild(positions);
}

void SurfaceRenderer::Update(const std::vector<Float3>& positions)
{
	// Building takes ~1 s for 100k+ atoms, so during a simulation it may only take a fifth of the render thread's time.
	// Small systems still rebuild for every new frame
	if (std::chrono::steady_clock::now() - lastBuildEnd < lastBuildDuration * 4)
		return;
	Rebuild(positions);
}

void SurfaceRenderer::Rebuild(const std::vector<Float3>& positions)
{
	const auto start = std::chrono::steady_clock::now();
	const Mesh mesh = BuildMesh(atoms, positions);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, mesh.vertices.size() * sizeof(Vertex), mesh.vertices.data(), GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	// The element buffer binding is part of the vertex array
	glBindVertexArray(vao);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.indices.size() * sizeof(unsigned int), mesh.indices.data(), GL_DYNAMIC_DRAW);
	glBindVertexArray(0);
	indexCount = static_cast<int>(mesh.indices.size());
	lastBuildEnd = std::chrono::steady_clock::now();
	lastBuildDuration = lastBuildEnd - start;
}

void SurfaceRenderer::Clear()
{
	atoms.clear();
	indexCount = 0;
}

void SurfaceRenderer::Draw(DrawSurfaceShader& shader, const glm::mat4& view, const glm::mat4& projection) const
{
	if (indexCount == 0)
		return;
	shader.Draw(vao, indexCount, view, projection);
}

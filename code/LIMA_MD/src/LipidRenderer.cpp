#include "LipidRenderer.h"

#include "Shaders.h"
#include "TubeMesh.h"

#include <algorithm>

namespace {

constexpr std::size_t controlPointStride = 2;
constexpr int samplesPerSegment = 4;
constexpr int tubeSides = 10;
constexpr int capRings = 4;
// Brighter than the cartoon's, keeps the pastel tails from going muddy on their shaded side
constexpr float ambient = 0.5f;
constexpr float tailRadius = 0.06f;	// [nm] Chunky rather than to scale, so the tails read as soft noodles
constexpr float sterolRadius = 0.11f;	// [nm] The ring system is flat and wide, a thick tube reads better than a thin one

// Tails take a pale version of their head's color, so each lipid reads as one piece
glm::vec4 TailColor(LipidHeadgroup headgroup)
{
	const glm::vec4 head = LipidRenderer::HeadColor(headgroup);
	const float whiteness = headgroup == LipidHeadgroup::Sterol ? 0.35f : 0.6f;
	return glm::mix(head, glm::vec4(1.f), whiteness);
}

} // namespace

glm::vec4 LipidRenderer::HeadColor(LipidHeadgroup headgroup)
{
	switch (headgroup) {
	case LipidHeadgroup::Choline: return { 1.00f, 0.71f, 0.52f, 1.f };	// peach
	case LipidHeadgroup::Ethanolamine: return { 0.52f, 0.76f, 1.00f, 1.f };	// sky
	case LipidHeadgroup::Serine: return { 0.76f, 0.62f, 1.00f, 1.f };	// lavender
	case LipidHeadgroup::Anionic: return { 0.56f, 0.88f, 0.62f, 1.f };	// mint
	case LipidHeadgroup::Sterol: return { 1.00f, 0.58f, 0.70f, 1.f };	// pink
	}
	return { 0.8f, 0.8f, 0.8f, 1.f };
}


void LipidRenderer::Prepare(const LipidModels& newLipids, const std::vector<Float3>& positions, Float3 newBoxSize)
{
	boxSize = newBoxSize;
	lipids = newLipids;
	meshes.assign(lipids.size(), {});
	for (std::size_t i = 0; i < lipids.size(); ++i) {
		meshes[i].color = TailColor(lipids[i].headgroup);
		// Picking a tail selects the lipid, the same as picking its head
		meshes[i].objectId = lipids[i].headParticleId;
	}
	Update(positions);
}

void LipidRenderer::Update(const std::vector<Float3>& positions)
{
	vertices.clear();
	std::vector<glm::vec3> path;
	for (std::size_t lipidIndex = 0; lipidIndex < lipids.size(); ++lipidIndex) {
		const LipidModel& lipid = lipids[lipidIndex];
		Mesh& mesh = meshes[lipidIndex];
		mesh.firstVertex = vertices.size();
		mesh.vertexCount = 0;
		if (lipid.headParticleId < 0 || static_cast<std::size_t>(lipid.headParticleId) >= positions.size())
			continue;

		// Unwrap around the head, as the head bead is drawn at its wrapped position
		const glm::vec3 head = TubeMesh::ToVec3(positions[lipid.headParticleId]);
		const float radius = lipid.headgroup == LipidHeadgroup::Sterol ? sterolRadius : tailRadius;
		for (const std::vector<int>& tail : lipid.tails) {
			if (std::ranges::any_of(tail, [&positions](int id) { return id < 0 || static_cast<std::size_t>(id) >= positions.size(); }))
				continue;
			const auto AddPoint = [&](int particleId) {
				path.push_back(head + TubeMesh::MinimumImage(TubeMesh::ToVec3(positions[particleId]) - head, boxSize));
			};
			// Every other atom, so the spline follows the chain rather than its zigzag
			path.clear();
			for (std::size_t i = 0; i < tail.size(); i += controlPointStride)
				AddPoint(tail[i]);
			if ((tail.size() - 1) % controlPointStride != 0)
				AddPoint(tail.back());
			if (path.size() < 2)
				continue;
			const std::vector<TubeMesh::Sample> samples = TubeMesh::SamplePath(path, 0, path.size() - 1, samplesPerSegment);
			TubeMesh::AppendTube(vertices, samples, radius, tubeSides);
			TubeMesh::AppendRoundCap(vertices, samples.back(), radius, tubeSides, capRings);
		}
		mesh.vertexCount = vertices.size() - mesh.firstVertex;
	}
}

void LipidRenderer::Clear()
{
	lipids.clear();
	meshes.clear();
	vertices.clear();
}

void LipidRenderer::Draw(DrawTrianglesShader& shader, const glm::mat4& viewProjection) const
{
	// One upload for all lipids, but a draw per lipid, as each has its own color and picking id
	shader.Upload(vertices);
	const glm::mat4 identity(1.f);
	for (const Mesh& mesh : meshes)
		shader.DrawUploaded(mesh.firstVertex, mesh.vertexCount, viewProjection, identity, mesh.color, mesh.objectId, ambient);
}

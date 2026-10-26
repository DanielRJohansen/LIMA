#pragma once

#include <cstdint>
#include <vector>

// What kind of molecule a particle belongs to, for selecting eg. "all lipids" in the display
enum class ParticleCategory : uint8_t { Other, Protein, Lipid, Water, Ion };

// Indexed by particle id in Box numbering
using ParticleCategories = std::vector<ParticleCategory>;

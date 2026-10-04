#pragma once

// Precompiled header for LIMA's host (C++) code, see lima_pch() in the root CMakeLists.txt.
// Only stable, widely used headers belong here: touching anything included below recompiles every C++ file.

// Standard library
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Third party
#include <cuda_runtime.h>
#include <glm.hpp>

// LIMA headers included by most files. Changing one of these recompiles every C++ file, which is close to what
// changing them did anyway
#include "LimaTypes.cuh"
#include "Bodies.cuh"
#include "Simulation.cuh"
#include "Filehandling.h"
#include "MDFiles.h"
#include "Utilities.h"
#include "TimeIt.h"

#pragma once

// CUB's block-wide algorithms (cub::BlockReduce). Under HIP, cub:: maps to hipCUB. Separate from GpuCubDevice.cuh,
// whose headers would add a kernel to every file including this one
#if defined(LIMA_HIP)
#include <hipcub/hipcub.hpp>
namespace cub = hipcub;
#else
#include <cub/block/block_reduce.cuh>
#endif

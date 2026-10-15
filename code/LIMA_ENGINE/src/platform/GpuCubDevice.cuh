#pragma once

// CUB's device-wide algorithms (cub::DeviceReduce, cub::DeviceScan) and cuda::std::plus. Under HIP they map to hipCUB
// and std (clang treats constexpr functions as host+device, so std::plus works in kernels)
#if defined(LIMA_HIP)
#include <hipcub/hipcub.hpp>
#include <functional>
namespace cub = hipcub;
namespace cuda { namespace std = ::std; }
#else
#include <cub/device/device_reduce.cuh>
#include <cub/device/device_scan.cuh>
#include <cuda/std/functional>
#endif

// HIP go/no-go spike: does a plain kernel and PME's FFT setup run on this GPU? Not part of the LIMA build.
//   CUDA: nvcc -O2 spike.cu -o spike.exe -lcufft
//   HIP:  hipcc -x hip --offload-arch=gfx1103 -O2 spike.cu -o spike.exe -lhipfft
//         ROCm 7.2 on Windows: from a VS 2022 (MSVC 14.44) developer shell, and run with %HIP_PATH%bin on PATH
// Prints the device, the checks, and ends with SPIKE_PASS or SPIKE_FAIL (exit code 0 or 1)

#include "GpuCompat.h"

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

namespace {
	bool allPassed = true;

	void Check(bool ok, const char* what) {
		std::printf("%s %s\n", ok ? "PASS" : "FAIL", what);
		allPassed = allPassed && ok;
	}
	void CheckGpu(cudaError_t status, const char* what) {
		if (status != cudaSuccess) {
			std::printf("FAIL %s: %s\n", what, cudaGetErrorString(status));
			std::printf("SPIKE_FAIL\n");
			std::exit(1);
		}
	}
	void CheckFft(cufftResult status, const char* what) {
		if (status != CUFFT_SUCCESS) {
			std::printf("FAIL %s: fft error %d\n", what, static_cast<int>(status));
			std::printf("SPIKE_FAIL\n");
			std::exit(1);
		}
	}
}

__global__ void IndexKernel(int* out, int* warpSizeOut) {
	const int id = blockIdx.x * blockDim.x + threadIdx.x;
	out[id] = threadIdx.x * blockIdx.x;
	if (id == 0) *warpSizeOut = warpSize;
}

static void TestKernel() {
	constexpr int nBlocks = 4, nThreads = 256, n = nBlocks * nThreads;
	int* out = nullptr;
	int* warpSizeDevice = nullptr;
	CheckGpu(cudaMalloc(&out, n * sizeof(int)), "malloc");
	CheckGpu(cudaMalloc(&warpSizeDevice, sizeof(int)), "malloc");
	IndexKernel<<<nBlocks, nThreads>>>(out, warpSizeDevice);
	CheckGpu(cudaGetLastError(), "kernel launch");
	CheckGpu(cudaDeviceSynchronize(), "kernel run");

	std::vector<int> host(n);
	int warpSizeHost = 0;
	CheckGpu(cudaMemcpy(host.data(), out, n * sizeof(int), cudaMemcpyDeviceToHost), "copy back");
	CheckGpu(cudaMemcpy(&warpSizeHost, warpSizeDevice, sizeof(int), cudaMemcpyDeviceToHost), "copy back");
	bool ok = true;
	for (int i = 0; i < n; i++)
		ok = ok && host[i] == (i % nThreads) * (i / nThreads);
	Check(ok, "kernel writes threadIdx.x * blockIdx.x");
	std::printf("warpSize %d\n", warpSizeHost);
	Check(warpSizeHost == 32, "warpSize is 32 (LIMA assumes wave32)");
	cudaFree(out);
	cudaFree(warpSizeDevice);
}

// Batched 3D R2C forward and C2R inverse, planned exactly as PME::Controller does. A rectangular grid, as LIMA's boxes
// may be, with x/2+1 odd
static void TestFft() {
	constexpr int nx = 32, ny = 28, nz = 24, batch = 2;
	constexpr int nReal = nx * ny * nz;
	constexpr int nComplex = (nx / 2 + 1) * ny * nz;
	int dimensions[3]{ nz, ny, nx };
	int complexDimensions[3]{ nz, ny, nx / 2 + 1 };

	cudaStream_t stream;
	CheckGpu(cudaStreamCreate(&stream), "stream");
	cufftHandle planForward, planInverse;
	CheckFft(cufftPlanMany(&planForward, 3, dimensions, dimensions, 1, nReal, complexDimensions, 1, nComplex, CUFFT_R2C, batch), "plan R2C");
	CheckFft(cufftPlanMany(&planInverse, 3, dimensions, complexDimensions, 1, nComplex, dimensions, 1, nReal, CUFFT_C2R, batch), "plan C2R");
	CheckFft(cufftSetStream(planForward, stream), "set stream");
	CheckFft(cufftSetStream(planInverse, stream), "set stream");

	std::vector<float> input(size_t(nReal) * batch);
	std::mt19937 rng(1234);
	std::uniform_real_distribution<float> dist(-1.f, 1.f);
	for (float& v : input) v = dist(rng);

	cufftReal* real = nullptr;
	cufftComplex* fourier = nullptr;
	CheckGpu(cudaMalloc(&real, input.size() * sizeof(float)), "malloc");
	CheckGpu(cudaMalloc(&fourier, size_t(nComplex) * batch * sizeof(cufftComplex)), "malloc");
	CheckGpu(cudaMemcpy(real, input.data(), input.size() * sizeof(float), cudaMemcpyHostToDevice), "upload");

	CheckFft(cufftExecR2C(planForward, real, fourier), "exec R2C");
	CheckGpu(cudaStreamSynchronize(stream), "R2C run");
	std::vector<cufftComplex> spectrum(size_t(nComplex) * batch);
	CheckGpu(cudaMemcpy(spectrum.data(), fourier, spectrum.size() * sizeof(cufftComplex), cudaMemcpyDeviceToHost), "download");

	// One mode of the second batch member against a direct DFT on the host
	const int kx = 3, ky = 5, kz = 7, b = 1;
	std::complex<double> reference = 0;
	for (int z = 0; z < nz; z++)
		for (int y = 0; y < ny; y++)
			for (int x = 0; x < nx; x++) {
				const double phase = -2.0 * 3.141592653589793 * (double(kx) * x / nx + double(ky) * y / ny + double(kz) * z / nz);
				reference += double(input[size_t(b) * nReal + (z * ny + y) * nx + x]) * std::polar(1.0, phase);
			}
	const cufftComplex gpu = spectrum[size_t(b) * nComplex + (kz * ny + ky) * (nx / 2 + 1) + kx];
	const double modeError = std::abs(std::complex<double>(gpu.x, gpu.y) - reference) / std::abs(reference);
	std::printf("mode (%d,%d,%d) relative error %.2e\n", kx, ky, kz, modeError);
	Check(modeError < 1e-4, "R2C mode matches host DFT");

	CheckFft(cufftExecC2R(planInverse, fourier, real), "exec C2R");
	CheckGpu(cudaStreamSynchronize(stream), "C2R run");
	std::vector<float> roundTrip(input.size());
	CheckGpu(cudaMemcpy(roundTrip.data(), real, roundTrip.size() * sizeof(float), cudaMemcpyDeviceToHost), "download");
	double maxError = 0;
	for (size_t i = 0; i < input.size(); i++)
		maxError = std::fmax(maxError, std::fabs(roundTrip[i] / nReal - input[i]));
	std::printf("round trip max abs error %.2e\n", maxError);
	Check(maxError < 1e-5, "R2C + C2R round trip returns the input (scaled by N)");

	cufftDestroy(planForward);
	cufftDestroy(planInverse);
	cudaFree(real);
	cudaFree(fourier);
	cudaStreamDestroy(stream);
}

int main() {
	cudaDeviceProp properties{};
	CheckGpu(cudaGetDeviceProperties(&properties, 0), "device properties");
	int driver = 0, runtime = 0;
	cudaDriverGetVersion(&driver);
	cudaRuntimeGetVersion(&runtime);
#if defined(__HIP__)
	std::printf("%s %s (%s), driver %d, runtime %d\n", LIMA_GPU_VENDOR, properties.name, properties.gcnArchName, driver, runtime);
#else
	std::printf("%s %s (sm_%d%d), driver %d, runtime %d\n", LIMA_GPU_VENDOR, properties.name, properties.major, properties.minor, driver, runtime);
#endif

	TestKernel();
	TestFft();

	std::printf(allPassed ? "SPIKE_PASS\n" : "SPIKE_FAIL\n");
	return allPassed ? 0 : 1;
}

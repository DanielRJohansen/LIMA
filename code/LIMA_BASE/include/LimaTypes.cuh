#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <span>
#include "Constants.h"

#include <array>


enum ATOM_TYPE { NONE, O, C, P, N, H, SOL, S, LIMA_CUSTOM};


struct Int3 {
	__host__ __device__ constexpr Int3() {}
	__host__ __device__ constexpr Int3(const int& x, const int& y, const int& z) : x(x), y(y), z(z) {}

	__host__ __device__ inline Int3 operator + (const Int3 a) const { return Int3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Int3 operator - (const Int3 a) const { return Int3(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ inline Int3 operator * (const int a) const { return Int3(x * a, y * a, z * a); }
	__host__ __device__ inline Int3 operator / (const int a) const { return Int3(x / a, y / a, z / a); }

	__host__ __device__ void operator += (const Int3& a) { x += a.x; y += a.y; z += a.z; }
	__host__ __device__ void operator /= (const int a) { x /= a; y /= a; z /= a; }

	__host__ __device__ bool operator== (const Int3& a) const { return (x == a.x && y == a.y && z == a.z); }
	__host__ __device__ bool operator!= (const Int3& a) const { return (x != a.x || y != a.y || z != a.z); }

	__host__ __device__ int manhattanLen() const { return std::abs(x) + std::abs(y) + std::abs(z); }
	__device__ int MaxAbsElement() const { return std::max(std::abs(x), std::max(std::abs(y), std::abs(z))); }
	__host__ Int3 abs() const { return Int3{ std::abs(x), std::abs(y), std::abs(z) }; }

	__host__ __device__ void print(char c = '_', bool prefix_newline = false) const {
		char nl = prefix_newline ? '\n' : ' ';
		printf("%c %c %d\t %d\t %d\n", nl, c, x, y, z);
	}

	std::string toString() const {
		return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
	}

	int x = 0, y = 0, z = 0;
};



struct Float3 {
	float x = 0, y = 0, z = 0;

	__host__ __device__ Float3() {}
	__host__ __device__ Float3(float a) : x(a), y(a), z(a) {}
	__host__ __device__ Float3(float x, float y, float z) : x(x), y(y), z(z) {}
	__host__ __device__ Float3(float* a) { x = a[0]; y = a[1]; z = a[2]; }
	__host__ __device__ Float3(int a) : x(static_cast<float>(a)), y(static_cast<float>(a)), z(static_cast<float>(a)) {}
	__host__ __device__ Float3(const int& x, const int& y, const int& z) : x(static_cast<float>(x)), y(static_cast<float>(y)), z(static_cast<float>(z)) {}
	__host__ Float3(const double& x, const double& y, const double& z) : x(static_cast<float>(x)), y(static_cast<float>(y)), z(static_cast<float>(z)) {}

	__host__ __device__ inline Float3 operator - () const { return Float3(-x, -y, -z); }
	__host__ __device__ inline Float3 operator * (const float a) const { return Float3(x * a, y * a, z * a); }
	__host__ __device__ inline Float3 operator * (const Float3& a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	__host__ __device__ inline Float3 operator / (const float a) const { return Float3(x / a, y / a, z / a); }
	__host__ __device__ inline Float3 operator / (const Float3& a) const { return Float3(x / a.x, y / a.y, z / a.z); }
	__host__ __device__ inline Float3 operator + (const Float3& a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Float3 operator - (const Float3& a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ inline bool operator == (const Float3& a) const { return (a.x == x && a.y == y && a.z == z); }
	__host__ __device__ inline void operator += (const Float3& a) { x += a.x; y += a.y; z += a.z; }
	__host__ __device__ inline void operator -= (const Float3& a) { x -= a.x; y -= a.y; z -= a.z; }
	__host__ __device__ inline void operator *= (const float a) { x *= a; y *= a; z *= a; }

	__host__ __device__ inline bool operator < (const Float3 a) const { return x < a.x&& y < a.y&& z < a.z; }
	__host__ __device__ inline bool operator > (const Float3 a) const { return x > a.x && y > a.y && z > a.z; }

	float* begin() { return &x; }
	const float* begin() const { return &x; }
	float* end() { return &x + 3; }
	const float* end() const { return &x + 3; }

	__host__ __device__ inline float operator[] (int index) const {
		switch (index) {
			case 0:
				return x;
			case 1:
				return y;
			case 2:
				return z;
			default:
				return -404;
		}
	}

	__host__ __device__ inline float& operator[] (int index) {
		switch (index) {
		default:			// Sadly it is not reasonable to catch this in release
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		}
	}

	__host__ __device__ Float3 norm_d() const {
		double l = len_d();
		if (l)
			return *this * (1. / l);
		return Float3{};
	}
	__host__ __device__ Float3 norm() const {
		const float l = len();
		if (l)
			return *this * (1.f / l);
		return Float3{};
	}
	__device__ Float3 norm_fast() const {		// Unsafe, may divide by 0
		return *this * (1.f / len());
	}

	__device__ bool isNan() const {
		return isnan(x) || isnan(y) || isnan(z);
	}
	__host__ __device__ Float3 round() const { return Float3{ roundf(x), roundf(y), roundf(z) }; }
	__host__ __device__ Float3 square() const { return Float3(x * x, y * y, z * z); }
	__host__ __device__ inline float len() const { return (float)sqrtf(x * x + y * y + z * z); }
	__host__ __device__ inline double len_d() const { return sqrt((double)x * x + (double)y * y + (double)z * z); }
	__host__ __device__ inline float lenSquared() const { return (float)(x * x + y * y + z * z); }
	__host__ __device__ Float3 zeroIfAbove(float a) { return Float3(x * (x < a), y * (y < a), z * (z < a)); }
	__host__ __device__ Float3 zeroIfBelow(float a) { return Float3(x * (x > a), y * (y > a), z * (z > a)); }

	__host__ __device__ Float3 Floor() { return Float3(floorf(x), floorf(y), floorf(z));}

	__host__ __device__ inline static float getAngle(const Float3& v1, const Float3& v2) {
		float val = (v1.dot(v2)) / (v1.len() * v2.len());	// If i make this float, we get values over 1, even with the statements below! :(
		//if (val > 1.f || val < -1.f) { printf("Val1 %f !!\n", val);}
		val = val > 1.f ? 1.f : val;
		val = val < -1.f ? -1.f : val;
		return acos(val);
	}
	__device__ static float getAngleOfNormVectors(const Float3& v1, const Float3& v2) {
		float val = v1.dot(v2);
		val = val > 1.f ? 1.f : val;
		val = val < -1.f ? -1.f : val;
		return acosf(val);
	}

	__host__ __device__ static float getAngle(Float3 a, Float3 middle, Float3 b) {
		return getAngle(a - middle, b - middle);
	}


	__host__ __device__ Float3 cross(const Float3 a) const { 
		return Float3(
			y * a.z - z * a.y, 
			z * a.x - x * a.z, 
			x * a.y - y * a.x); }
	__host__ __device__ float dot(Float3 a) const { return (x * a.x + y * a.y + z * a.z); }
	__host__ __device__ Float3 abs() const {
		return Float3(
			std::abs(x),
			std::abs(y),
			std::abs(z)
		);
	}


	__host__ __device__ void print(char c = '_', bool prefix_newline=false) const {
		char nl = prefix_newline ? '\n' : ' ';
		if (len() < 100)
			printf("%c %c %.12f %.12f %.12f\n",nl, c, x, y, z);
		else if (len() < 500'000) {
			printf("%c %c %.4f %.4f %.4f\n", nl, c, x, y, z);
		}
		else
			printf("%c %c %.0f\t %.0f\t %.0f\n",nl, c, x, y, z);
	}

	std::string toString() const {
		//return std::format("{} {} {}", x, y, z);
		return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
	}

	__host__ __device__ void rotateAroundOrigo(Float3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// pitch and yaw is relative to global coordinates. 

		*this = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);
		*this = rodriguesRotatation(*this, Float3(0, 1, 0), pitch_yaw_roll.y);
		*this = rodriguesRotatation(*this, Float3(0, 0, 1), pitch_yaw_roll.z);
	}

	__host__ __device__ static Float3 rodriguesRotatation(const Float3 v, const Float3 k, const float theta) {
		return v * cos(theta) + k.cross(v) * sin(theta) + k * (k.dot(v)) * (1.f - cos(theta));
	}

	__device__ __host__ inline float LargestMagnitudeElement() const {
		return std::max(
			std::max(std::abs(x), std::abs(y)),
			std::abs(z)
		);
	}


	// Not used right now!
	__host__ __device__ static Float3 centerOfMass(Float3* arr_ptr, uint32_t arr_size) {	// Only run before sim, so we can cast to double without slowing sim
		Float3 sum = Float3(0, 0, 0);
		for (uint32_t i = 0; i < arr_size; i++) {
			sum = sum + arr_ptr[i];
		}
		return sum * (1.f / static_cast<float>(arr_size));
	}

	// Assumes *this is a point, and arguments are two points on a line
	__host__ __device__ float distToLine(const Float3& p1, const Float3& p2) const {
		return ((p2 - p1).cross(p1 - (*this))).len() / (p2 - p1).len();
	}

	// Concat to string operation
	friend std::string& operator<<(std::string& str, const Float3& f) {
		str = str + std::to_string(f.x) + " " + std::to_string(f.y) + " " + std::to_string(f.z);
		return str;
	}

};


struct Double3 {
	__host__ __device__ Double3() {}
	__host__ __device__ Double3(double a) : x(a), y(a), z(a) {}
	__host__ __device__ Double3(double x, double y, double z) : x(x), y(y), z(z) {}
	__host__ __device__ Double3(Float3 a) : x((double)a.x), y((double)a.y), z((double)a.z) {}

	__host__ __device__ inline Double3 operator + (const Float3 a) const {
		return Double3(x + (double)a.x, y + (double)a.y, z + (double)a.z);
	}
	__host__ __device__ inline Double3 operator + (const Double3 a) const { return Double3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Double3 operator / (const double a) const { return Double3(x / a, y / a, z / a); }
	__host__ __device__ inline void operator += (const Float3 a) { x += (double)a.x; y += (double)a.y; z += (double)a.z; }
	__host__ __device__ inline void operator += (const Double3 a) { x += a.x; y += a.y; z += a.z; }
	__host__ __device__ inline Double3 operator - (const Double3 a) const { return Double3(x - a.x, y - a.y, z - a.z); }

	__host__ __device__ inline double len() { return (double)sqrt(x * x + y * y + z * z); }

	__host__ __device__ Float3 toFloat3() const {
		return Float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
	}

	__host__ __device__ void print(char c = '_') {
		printf("%c %.10f %.10f %.10f\n", c, x, y, z);
	}

	double x = 0, y = 0, z = 0;
};

// LIMA Coordinate3
struct Coord {
	int32_t x = 0, y = 0, z = 0;	// [lm]

	__device__ __host__ Coord() {};
	__device__ __host__ explicit Coord(Float3 pos_abs) {
		x = static_cast<int32_t>(pos_abs.x);
		y = static_cast<int32_t>(pos_abs.y);
		z = static_cast<int32_t>(pos_abs.z);
	}

	__host__ __device__ Coord(int32_t a) : x(a), y(a), z(a) {}
	__host__ __device__ Coord(int32_t x, int32_t y, int32_t z) : x(x), y(y), z(z) {}

	__host__ __device__ Coord operator + (const Coord& a) const { return Coord(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ Coord operator - (const Coord& a) const { return Coord(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ Coord operator / (const int32_t& a) const { return Coord(x / a, y / a, z / a); }
	__host__ __device__ Coord operator - () const { return Coord(-x, -y, -z); }
	__host__ __device__ void operator += (const Coord& a) { x += a.x; y += a.y; z += a.z; };
	__host__ __device__ void operator -= (const Coord& a) { x -= a.x; y -= a.y; z -= a.z; };
	__host__ __device__ Coord operator * (const int32_t a) const { return Coord{ x * a, y * a, z * a }; }
	__host__ __device__ Coord operator * (const float) const = delete;
	__host__ __device__ bool operator == (const Coord& a) const { return x == a.x && y == a.y && z == a.z; }
	__host__ __device__ bool operator != (const Coord& a) const { return !(*this == a); }
	//__device__ Coord operator >> (const uint32_t a) const { return Coord(x >> a, y >> a, z >> a); }
	
	__host__ __device__ int32_t dot(const Coord& a) const { return (x * a.x + y * a.y + z * a.z); }
	__host__ __device__ void print(char c = '_', bool nl=1) const { 
		if (nl) printf(" %c %d %d %d\n", c, x, y, z);
		else printf(" %c %d %d %d", c, x, y, z);
	}
	// Print in pico, assuming baseline is lima
	__host__ __device__ void printS(char c = '_') const { 
		printf(" %c %d %d %d [pico]\n", c, x / 100000, y / 100000, z / 100000); }
	__host__ __device__ bool isZero() const { return (x == 0 && y == 0 && z == 0); }

	__device__ __host__ static int32_t max(int l, int r) { return l > r ? l : r; }

	__device__ __host__ int32_t maxElement() const { return max(std::abs(x), max(std::abs(y), std::abs(z))); }

	__host__ int32_t* get(int dim) {
		switch (dim)
		{
		case 0:
			return &x;
		case 1:
			return &y;
		case 2:
			return &z;
		default:
			throw std::runtime_error("Requested bad dimension");
		}
	}

	__host__ __device__ Float3 toFloat3() const { 
		return Float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)); 
	}

};

struct NodeIndex : public Int3 {
	__host__ __device__ NodeIndex() : Int3() {}
	__host__ __device__ NodeIndex(const int& x, const int& y, const int& z) : Int3(x,y,z) {}
	__host__ __device__ NodeIndex(const Int3& a) : Int3(a) {}

	//__host__ __device__ int32_t dot(const NodeIndex& a) const { return (x * a.x + y * a.y + z * a.z); }

	// This function does NOT return anything position related, only distance related
	__host__ __device__ Float3 toFloat3() const {
		return Float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
	}

	__device__ bool isZero() const { return (x == 0 && y == 0 && z == 0); }

	__device__ int sum() const { return std::abs(x) + std::abs(y) + std::abs(z); }

	__host__ int largestMagnitudeElement() const {
		const Int3 m = this->abs();
		return std::max(
			std::max(m.x, m.y),
			m.z
		);
	}

	__device__ __host__ bool isInBox(int nodes_per_dim) const {
		if (x < 0 || y < 0 || z < 0 || x >= nodes_per_dim || y >= nodes_per_dim || z >= nodes_per_dim)
			return false;
		return true;
	}
};





struct BoundingBox {
	BoundingBox(
		Float3 min = Float3{ std::numeric_limits<float>::max() },
		Float3 max = Float3{ std::numeric_limits<float>::lowest() }) 
		: min(min), max(max) {}


	Float3 min, max;

	bool intersects(BoundingBox b) {
		return
			min.x <= b.max.x && max.x >= b.min.x &&
			min.y <= b.max.y && max.y >= b.min.y &&
			min.z <= b.max.z && max.z >= b.min.z;
	}
	bool pointIsInBox(Float3 point) {
		return (min < point) && (point < max);
	}
	void addPadding(float margin) {
		min += Float3(-margin);
		max += Float3(margin);
	}
};

template <int len>
class FixedSizeMatrix{
public:
	__device__ FixedSizeMatrix() {}
	__host__ FixedSizeMatrix(bool val) {
		//uint8_t bit = val ? 1 : 0;
		for (int i = 0; i < m_size; i++) {
			matrix[i] = val ? 0xFF : 0;
		}
	}

	__host__ __device__ bool get(int i1, int i2) const {
		int index = i1 + i2 * m_len;
		int byteIndex = index / 8;
		int bitIndex = index % 8;
		return (matrix[byteIndex] >> bitIndex) & 1U;
	}

	__host__ void set(int i1, int i2, bool val) {
		int index = i1 + i2 * m_len;
		int byteIndex = index / 8;
		int bitIndex = index % 8;
		if (val)
			matrix[byteIndex] |= (1U << bitIndex);
		else
			matrix[byteIndex] &= ~(1U << bitIndex);
	}

	__device__ void load(const FixedSizeMatrix<len>& src) {
		for (int i = threadIdx.x; i < m_size; i += blockDim.x) {
			matrix[i] = src.matrix[i];
		}
	}


	__host__ void printMatrix(int n) const {
		// Print column indices
		std::cout << "     ";  // Space for row indices
		for (int j = 0; j < n; ++j) { 
			std::string separator = j + 1 < 10 ? "  " : " ";
			std::cout << j+1 << separator;
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

private:
	const static int m_len = len;
	const static int m_size = (m_len * m_len + 7) / 8; // Ceil division
	uint8_t matrix[m_size]{};
};

using BondedParticlesLUT = FixedSizeMatrix<MAX_COMPOUND_PARTICLES>;

class BondedParticlesLUTManager {
	static const int max_bonded_compounds = 5;	// first 3: self, res-1 and res+1. The rest are various h bonds i think
	static const int n_elements = MAX_COMPOUNDS * max_bonded_compounds;

	BondedParticlesLUT luts[n_elements];
	uint32_t connected_compound_ids_masks[MAX_COMPOUNDS];

	__device__ __host__ int getLocalIndex(int id_self, int id_other) {
		return (max_bonded_compounds/2) + (id_other - id_self);
	}
	__device__ __host__ int getGlobalIndex(int local_index, int id_self) {
		return id_self*max_bonded_compounds + local_index;
	}
	__device__ __host__ uint32_t getMask(int index) {
		return 1u << index;
	}

public:
	BondedParticlesLUTManager() {
		for (int i = 0; i < n_elements; i++) {
			luts[i] = BondedParticlesLUT(false);
		}
	}

	__device__ BondedParticlesLUT* get(int id_self, int id_other){
		// The around around when this function is called on device, should ensure 
		// that there is always an entry in the table for the 2 compounds 
		return &luts[getGlobalIndex(getLocalIndex(id_self, id_other), id_self)];
	}

	__host__ BondedParticlesLUT* get(int id_self, int id_other, bool) {
		if (std::abs(id_self - id_other > 2)) {
			throw std::runtime_error("Cannot get BPLUT for compounds with distanecs > 2 in id-space");
		}

		const int local_index = getLocalIndex(id_self, id_other);
		if (connected_compound_ids_masks[id_self] & getMask(local_index)) {
			return &luts[getGlobalIndex(local_index, id_self)];
		}
		return nullptr;
	}

	__host__ void addNewConnectedCompoundIfNotAlreadyConnected(int id_self, int id_other) {

		if (std::abs(id_self - id_other > 2)) {
			throw std::runtime_error("Cannot connect compounds that are too far apart in id space");
		}

		connected_compound_ids_masks[id_self] |= getMask(getLocalIndex(id_self, id_other));
	}
};








template<typename T>
T* genericMoveToDevice(T* data_ptr, int n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
	if (n_elements == 0) { return nullptr; }

	T* gpu_ptr = nullptr;
	size_t bytesize = n_elements * sizeof(T);

	cudaMallocManaged(&gpu_ptr, bytesize);
	auto cuda_status = cudaMemcpy(gpu_ptr, data_ptr, bytesize, cudaMemcpyHostToDevice);

	if (cuda_status != cudaSuccess) {
		std::cout << "\nCuda error code: " << cuda_status << " - " << cudaGetErrorString(cuda_status) << std::endl;
		throw std::runtime_error("Move to device failed");
	}

	cudaDeviceSynchronize();

	if (n_elements == 1)
		delete data_ptr;
	else
		delete[] data_ptr;

	//data_ptr = gpu_ptr;

	//printf("Moved %.2f MB to device\n", bytesize*1e-6);
	return gpu_ptr;
}

// TODO MEMLEAK: this function is not freeing the memory of the original data
// Assumes data is a ptr to device memory. Will copy what was in data, allocate new memory for data and move the copied data into that
template<typename T>
void genericCopyToHost(T** data, uint32_t n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
	T* data_host = new T[n_elements];

	size_t bytesize = sizeof(T) * n_elements;
	cudaMemcpy(data_host, *data, bytesize, cudaMemcpyDeviceToHost);

	*data = data_host;
}

template<typename T>
void GenericCopyToHost(T* srcDevice, std::vector<T>& destHost, size_t nElements) {
	destHost.resize(nElements);
	cudaMemcpy(destHost.data(), srcDevice, nElements * sizeof(T), cudaMemcpyDeviceToHost);
}

template<typename T>
void genericCopyToDevice(const T& src, T** dest, int n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
	size_t bytesize = n_elements * sizeof(T);

	cudaMallocManaged(dest, bytesize);
	cudaMemcpy(*dest, &src, bytesize, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}



 
// Order is critical, as many times something like "bool verbose = vl > V1" occurs
enum VerbosityLevel {
	SILENT,
	CRITICAL_INFO, 
	V1,
	V2,
	V3
};

enum EnvMode { Full, ConsoleOnly, Headless };

struct RenderAtom {
	float4 position = Disabled(); // {posX, posY, posZ, radius} [nm]
	float4 color{};					// {r, g, b, a} [0-1]	

	bool IsDisabled() const { return position.x == std::numeric_limits<float>::max() && position.y == std::numeric_limits<float>::max() && position.z == std::numeric_limits<float>::max(); }
	__device__ __host__ static constexpr float4 Disabled() { return float4{ std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() }; }
};
#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <math.h>
#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <span>
#include <optional>
#include "Constants.h"
#include <array>
#include <ranges>
//#include <generator>

// TODO: EASY: LARGE: Its a huge waste that the boxsize is a Int3, when it really should be a packed into a single 32 bit DWORD..
// 1024 nm boxsize is a reasonable limitation. However we cant use the same type for PME grid obviously
struct Int3 {
	constexpr Int3() {}
	constexpr Int3(const int& x, const int& y, const int& z) : x(x), y(y), z(z) {}

	constexpr Int3 operator + (const Int3 a) const { return Int3(x + a.x, y + a.y, z + a.z); }
	constexpr Int3 operator - (const Int3 a) const { return Int3(x - a.x, y - a.y, z - a.z); }
	constexpr Int3 operator * (const int a) const { return Int3(x * a, y * a, z * a); }
	constexpr Int3 operator / (const int a) const { return Int3(x / a, y / a, z / a); }
	constexpr Int3 operator - () const { return Int3(-x, -y, -z); }

	constexpr void operator += (const Int3& a) { x += a.x; y += a.y; z += a.z; }
	constexpr void operator /= (const int a) { x /= a; y /= a; z /= a; }

	constexpr bool operator== (const Int3& a) const { return (x == a.x && y == a.y && z == a.z); }
	constexpr bool operator!= (const Int3& a) const { return (x != a.x || y != a.y || z != a.z); }

	__host__ __device__ int manhattanLen() const { return std::abs(x) + std::abs(y) + std::abs(z); }
	__device__ int MaxAbsElement() const { return std::max(std::abs(x), std::max(std::abs(y), std::abs(z))); }
	__device__ __host__ Int3 abs() const { return Int3{ std::abs(x), std::abs(y), std::abs(z) }; }
	constexpr int InnerProduct() const { return x * y * z; }
	constexpr int Min() const { return std::min(x, std::min(y, z)); }

	__device__ __host__ void print(char c = '_', bool prefix_newline = false) const {
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

	constexpr Float3() {}
	constexpr explicit Float3(float a) : x(a), y(a), z(a) {}
	constexpr Float3(float x, float y, float z) : x(x), y(y), z(z) {}
	constexpr explicit Float3(int a) : x(static_cast<float>(a)), y(static_cast<float>(a)), z(static_cast<float>(a)) {}
	constexpr explicit Float3(const int& x, const int& y, const int& z) : x(static_cast<float>(x)), y(static_cast<float>(y)), z(static_cast<float>(z)) {}
	constexpr explicit Float3(const double& x, const double& y, const double& z) : x(static_cast<float>(x)), y(static_cast<float>(y)), z(static_cast<float>(z)) {}
	constexpr explicit Float3(const float4& a) : x(a.x), y(a.y), z(a.z) {}


	constexpr Float3 operator - () const { return Float3(-x, -y, -z); }
	constexpr Float3 operator * (const float a) const { return Float3(x * a, y * a, z * a); }
	constexpr Float3 operator * (const Float3& a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	constexpr Float3 operator / (const float a) const { return Float3(x / a, y / a, z / a); }
	constexpr Float3 operator / (const Float3& a) const { return Float3(x / a.x, y / a.y, z / a.z); }
	constexpr Float3 operator + (const Float3& a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	constexpr Float3 operator - (const Float3& a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	constexpr bool operator == (const Float3& a) const { return (a.x == x && a.y == y && a.z == z); }
	constexpr bool operator != (const Float3& a) const { return !(*this == a); }
	constexpr void operator += (const Float3& a) { x += a.x; y += a.y; z += a.z; }
	constexpr void operator -= (const Float3& a) { x -= a.x; y -= a.y; z -= a.z; }
	constexpr void operator *= (const float a) { x *= a; y *= a; z *= a; }
	constexpr bool operator < (const Float3 a) const { return x < a.x&& y < a.y&& z < a.z; }
	constexpr bool operator > (const Float3 a) const { return x > a.x && y > a.y && z > a.z; }

	constexpr float3 Tofloat3() const { return float3{ x, y, z }; }
	constexpr float4 Tofloat4(float w=0) const { return float4{ x, y, z, w }; }
	constexpr Int3 ToInt3() const { return Int3{ static_cast<int>(x), static_cast<int>(y), static_cast<int>(z) }; }
	__host__ static Float3 FromInt3(const Int3& a) { return Float3{ static_cast<float>(a.x), static_cast<float>(a.y), static_cast<float>(a.z) }; }


	__host__ inline float operator[] (int index) const {
		switch (index) {
		default:
			throw std::runtime_error("Illegal index in operator");
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
		}
	}

	__host__ inline float& operator[] (int index) {
		switch (index) {
		default:
			throw std::runtime_error("Illegal index in operator");
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
		//norm3d(x, y, z): TODO: OPTIM: Use the cuda math api for this and similar functions
		const float l = len();
		if (l)
			return *this * (1.f / l);
		return Float3{};
	}
	__device__ constexpr Float3 norm_fast() const {		// Unsafe, may divide by 0
		return *this * (1.f / len());
	}

	__device__ bool isNan() const {
		return isnan(x) || isnan(y) || isnan(z);
	}
	constexpr Float3 round() const { return Float3{ roundf(x), roundf(y), roundf(z) }; }
	constexpr Float3 square() const { return Float3(x * x, y * y, z * z); }
    __host__ __device__ inline float len() const { return std::sqrt(x * x + y * y + z * z); }
	__host__ __device__ inline double len_d() const { return sqrt((double)x * x + (double)y * y + (double)z * z); }
	constexpr float lenSquared() const { return (x * x + y * y + z * z); }
	constexpr Float3 zeroIfAbove(float a) { return Float3(x * (x < a), y * (y < a), z * (z < a)); }
	constexpr Float3 zeroIfBelow(float a) { return Float3(x * (x > a), y * (y > a), z * (z > a)); }
	constexpr Float3 sqrtElementwise() const { return Float3{ sqrtf(x), sqrtf(y), sqrtf(z) }; }


    constexpr Float3 Floor() const { return Float3(std::floor(x), std::floor(y), std::floor(z));}

	__host__ __device__ inline static float getAngle(const Float3& v1, const Float3& v2) {
		float val = (v1.dot(v2)) / (v1.len() * v2.len());	// If i make this float, we get values over 1, even with the statements below! :(
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


	constexpr Float3 cross(const Float3 a) const {
		return Float3(
			y * a.z - z * a.y, 
			z * a.x - x * a.z, 
			x * a.y - y * a.x); }
	constexpr float dot(Float3 a) const { return (x * a.x + y * a.y + z * a.z); }
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
	__host__ void print(const std::string& str) const {
		printf("%s %.6f %.6f %.6f\n", str.c_str(), x, y, z);
	}
	__host__ __device__ void print(int i) const {
		printf("%d %.6f %.6f %.6f\n", i, x, y, z);
	}

	std::string toString() const {
		//return std::format("{} {} {}", x, y, z);
		return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
	}

	// TODO: This is now proven wrong, and should be replaced, by a dedicated grometry namespace, using rotationmatrices instead
	__host__ __device__ void rotateAroundOrigo(Float3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// pitch and yaw is relative to global coordinates. 

		*this = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);
		*this = rodriguesRotatation(*this, Float3(0, 1, 0), pitch_yaw_roll.y);
		*this = rodriguesRotatation(*this, Float3(0, 0, 1), pitch_yaw_roll.z);
	}

	constexpr Float3 RotateAroundOrigo(Float3 pitch_yaw_roll) const {	//pitch around x, yaw around z, tilt around y
		// pitch and yaw is relative to global coordinates. 
		Float3 point = *this;
		point = rodriguesRotatation(point, Float3(1, 0, 0), pitch_yaw_roll.x);
		point = rodriguesRotatation(point, Float3(0, 1, 0), pitch_yaw_roll.y);
		point = rodriguesRotatation(point, Float3(0, 0, 1), pitch_yaw_roll.z);
		return point;
	}

	constexpr static Float3 rodriguesRotatation(const Float3 v, const Float3 k, const float theta) {
		return v * cos(theta) + k.cross(v) * sin(theta) + k * (k.dot(v)) * (1.f - cos(theta));
	}

	__device__ __host__ inline float LargestMagnitudeElement() const {
		return fmaxf(
			fmaxf(fabsf(x), fabsf(y)),
			fabsf(z)
		);
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

	static Float3 ElementwiseMin(const Float3& a, const Float3& b) {
		return Float3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
	}
	static Float3 ElementwiseMax(const Float3& a, const Float3& b) {
		return Float3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
	}

};

__device__ inline float4 Add(const float4& a, const float4& b) {
	return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

// Can only present integer values
struct Float3Compressed {

	Float3Compressed(){}
	Float3Compressed(const Float3& a) {
		const std::uint16_t x = (static_cast<int>(a.x) + 2) & 0x7u;
		const std::uint16_t y = (static_cast<int>(a.y) + 2) & 0x7u;
		const std::uint16_t z = (static_cast<int>(a.z) + 2) & 0x7u;
		data = static_cast<std::uint16_t>(x | (y << 3) | (z << 6));
	}

	constexpr Float3 Decode() const {
		const int x = static_cast<int>((data & 0x7u)) - 2;
		const int y = static_cast<int>(((data >> 3) & 0x7u)) - 2;
		const int z = static_cast<int>(((data >> 6) & 0x7u)) - 2;
		return Float3{ static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) };
	}

	uint16_t data;
};

struct ForceEnergy {
	Float3 force{};	// [J/mol/nm]
	float potE{};		// [J/mol]

	constexpr ForceEnergy operator+ (const ForceEnergy& a) const {
		return ForceEnergy{ force + a.force, potE + a.potE };
	}		
	constexpr void operator+= (const ForceEnergy& a) {
		force += a.force;
		potE += a.potE;
	}
	
	__host__ bool operator != (const ForceEnergy& a) const {
		return (force != a.force) || (potE != a.potE);
	}
};


struct Double3 {
	__host__ __device__ Double3() {}
	__host__ __device__ Double3(double a) : x(a), y(a), z(a) {}
	__host__ __device__ Double3(double x, double y, double z) : x(x), y(y), z(z) {}
	__host__ __device__ Double3(const Float3& a) : x((double)a.x), y((double)a.y), z((double)a.z) {}

	__host__ __device__ inline Double3 operator + (const Float3 a) const {
		return Double3(x + (double)a.x, y + (double)a.y, z + (double)a.z);
	}
	__host__ __device__ inline Double3 operator + (const Double3 a) const { return Double3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Double3 operator / (const double a) const { return Double3(x / a, y / a, z / a); }
	__host__ __device__ inline void operator += (const Float3 a) { x += (double)a.x; y += (double)a.y; z += (double)a.z; }
	__host__ __device__ inline void operator += (const Double3 a) { x += a.x; y += a.y; z += a.z; }
	__host__ __device__ inline Double3 operator - (const Double3 a) const { return Double3(x - a.x, y - a.y, z - a.z); }

	__host__ __device__ inline double len() const { return (double)sqrt(x * x + y * y + z * z); }

	__host__ __device__ Float3 toFloat3() const {
		return Float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
	}

	__host__ __device__ void print(char c = '_') const {
		printf("%c %.10f %.10f %.10f\n", c, x, y, z);
	}

	double x = 0, y = 0, z = 0;
};

struct NodeIndex : public Int3 {
	constexpr NodeIndex() : Int3() {}
	constexpr NodeIndex(const int& x, const int& y, const int& z) : Int3(x, y, z) {}
	constexpr NodeIndex(const Int3& a) : Int3(a) {}

	//constexpr NodeIndex operator+(const NodeIndex& a) const { return NodeIndex(x + a.x, y + a.y, z + a.z); }

	// This function does NOT return anything position related, only distance related
	constexpr Float3 toFloat3() const {
		return Float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
	}

	constexpr bool isZero() const { return (x == 0 && y == 0 && z == 0); }

	__device__ __host__ int Magnitude() const { return std::abs(x) + std::abs(y) + std::abs(z); }

	__device__ __host__ int largestMagnitudeElement() const {
		const Int3 m = this->abs();
		return std::max(
			std::max(m.x, m.y),
			m.z
		);
	}

	constexpr bool isInBox(const Int3& gridDim) const {
		if (x < 0 || y < 0 || z < 0 || x >= gridDim.x || y >= gridDim.y|| z >= gridDim.z)
			return false;
		return true;
	}
};



// Very fine-grained integer position, counted in [lm]
struct Coord {

	static const int nanoToLima_i = 100'000'000;
	static constexpr float nanoToLima_f = static_cast<float>(nanoToLima_i);
	static constexpr float limaToNano_f = static_cast<float>(1. / static_cast<double>(nanoToLima_f));

	int32_t x = 0, y = 0, z = 0;	// [lm]

	constexpr Coord() {};
	constexpr Coord(int32_t a) : x(a), y(a), z(a) {}
	constexpr Coord(int32_t x, int32_t y, int32_t z) : x(x), y(y), z(z) {}

	__device__ __host__ explicit Coord(Float3 pos) :
		x(static_cast<int32_t>(pos.x * nanoToLima_f)),
		y(static_cast<int32_t>(pos.y * nanoToLima_f)),
		z(static_cast<int32_t>(pos.z * nanoToLima_f))
		{
		if constexpr (POSITION_CHECKS) {
			if (std::abs(pos.x) > 1000.f || std::abs(pos.y) > 1000.f || std::abs(pos.z) > 1000.f) {// TODO magic nr,relates to intmax, fix!
				printf("pos %f %f %f\n", pos.x, pos.y, pos.z);
			}
		}
	}

	__device__ __host__ explicit Coord(const NodeIndex& nodeIndex) 
		: x(nodeIndex.x*nanoToLima_i), y(nodeIndex.y* nanoToLima_i), z(nodeIndex.z* nanoToLima_i) {
		if constexpr (POSITION_CHECKS) {
			if (std::abs(nodeIndex.x) > 1000 || std::abs(nodeIndex.y) > 1000 || std::abs(nodeIndex.z) > 1000) {// TODO magic nr,relates to intmax, fix!
				printf("NodeIndex %d %d %d\n", nodeIndex.x, nodeIndex.y, nodeIndex.z);
			}
		}
	}

	constexpr Float3 ToRelpos() const {
		return Float3{ static_cast<float>(x), static_cast<float>(y), static_cast<float>(z) } * limaToNano_f;
	}





	constexpr Coord operator + (const Coord& a) const { return Coord(x + a.x, y + a.y, z + a.z); }
	constexpr Coord operator - (const Coord& a) const { return Coord(x - a.x, y - a.y, z - a.z); }
	constexpr Coord operator / (const int32_t& a) const { return Coord(x / a, y / a, z / a); }
	constexpr Coord operator - () const { return Coord(-x, -y, -z); }
	constexpr void operator += (const Coord& a) { x += a.x; y += a.y; z += a.z; };
	constexpr void operator -= (const Coord& a) { x -= a.x; y -= a.y; z -= a.z; };
	constexpr Coord operator * (const int32_t a) const { return Coord{ x * a, y * a, z * a }; }
	constexpr Coord operator * (const float) const = delete;
	constexpr bool operator == (const Coord& a) const { return x == a.x && y == a.y && z == a.z; }
	constexpr bool operator != (const Coord& a) const { return !(*this == a); }
	
	constexpr int32_t dot(const Coord& a) const { return (x * a.x + y * a.y + z * a.z); }
	__host__ __device__ void print(char c = '_', bool nl=1) const { 
		if (nl) printf(" %c %d %d %d\n", c, x, y, z);
		else printf(" %c %d %d %d", c, x, y, z);
	}
	// Print in pico, assuming baseline is lima
	__host__ __device__ void printS(char c = '_') const { 
		printf(" %c %d %d %d [pico]\n", c, x / 100000, y / 100000, z / 100000); }
	constexpr bool isZero() const { return (x == 0 && y == 0 && z == 0); }

	__device__ __host__ int32_t maxElement() const { return std::max(std::abs(x), std::max(std::abs(y), std::abs(z))); }

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

	//__host__ __device__ Float3 toFloat3() const { 
	//	return Float3(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)); 
	//}
};

struct Rotation {
	Float3 center{};
	Float3 rotation{}; // rotX, rotY, rotZ
	constexpr bool Valid() const { return rotation != Float3{0.f}; }
};



struct BoundingBox {
	constexpr BoundingBox(
		Float3 min = Float3{ std::numeric_limits<float>::max() },
		Float3 max = Float3{ std::numeric_limits<float>::lowest() }) 
		: min(min), max(max) {}

	constexpr BoundingBox(const std::vector<Float3>& points);
	//BoundingBox(std::generator<Float3> generator); 
	BoundingBox(std::ranges::input_range auto&& range) {
		min = Float3{ std::numeric_limits<float>::max() };
		max = Float3{ std::numeric_limits<float>::min() };
		for (const Float3& p : range) {
			min.x = std::min(min.x, p.x);
			min.y = std::min(min.y, p.y);
			min.z = std::min(min.z, p.z);
			max.x = std::max(max.x, p.x);
			max.y = std::max(max.y, p.y);
			max.z = std::max(max.z, p.z);
		}
	}
	Float3 min, max;

	constexpr Float3 Center() const {
		return (min + max) * 0.5f;
	}
	constexpr Float3 Dimensions() const {
		return max - min;
	}
	constexpr bool intersects(BoundingBox b) const {
		return
			min.x <= b.max.x && max.x >= b.min.x &&
			min.y <= b.max.y && max.y >= b.min.y &&
			min.z <= b.max.z && max.z >= b.min.z;
	}
	constexpr bool pointIsInBox(Float3 point) const {
		return (min < point) && (point < max);
	}
	constexpr void addPadding(float padding) {
		min += Float3(-padding);
		max += Float3(padding);
	}
};

template<typename T>
void GenericCopyToHost(T* srcDevice, std::vector<T>& destHost, size_t nElements) {
	destHost.resize(nElements);
	cudaMemcpy(destHost.data(), srcDevice, nElements * sizeof(T), cudaMemcpyDeviceToHost);
}

template<typename T>
std::vector<T> GenericCopyToHost(T* srcDevice, size_t nElements) {
	std::vector<T> destHost(nElements);
	cudaMemcpy(destHost.data(), srcDevice, nElements * sizeof(T), cudaMemcpyDeviceToHost);
	return destHost;
}

template<typename T>
T GenericCopyToHost(T* srcDevice) {
	static_assert(std::is_trivially_copyable<T>::value, "GenericCopyToHost can only be used with trivially copyable types");
	T destHost;
	cudaMemcpy(&destHost, srcDevice, sizeof(T), cudaMemcpyDeviceToHost);
	return destHost;
}

template<typename T>
void genericCopyToDevice(const T& src, T** dest, int n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
	size_t bytesize = n_elements * sizeof(T);

	cudaMallocManaged(dest, bytesize);  // optim: THis shouldnt be managed
	cudaMemcpy(*dest, &src, bytesize, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}

template<typename T>
T* GenericCopyToDevice(const T* src, int n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
	T* dest;
	size_t bytesize = n_elements * sizeof(T);

	cudaMallocManaged(&dest, bytesize);
	cudaMemcpy(dest, src, bytesize, cudaMemcpyHostToDevice);

	return dest;
}

template<typename T>
T* GenericCopyToDevice(const std::vector<T>& src) {
	T* dest;
	cudaMalloc(&dest, src.size() * sizeof(T));
	cudaMemcpy(dest, src.data(), src.size() * sizeof(T), cudaMemcpyHostToDevice);
	return dest;
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

	RenderAtom() {}
	RenderAtom(Float3 positionNM, Float3 boxSize, char atomLetter);
	RenderAtom(float4 pos, float4 color) : position(pos), color(color) {};

	float4 position = Disabled(); // {posX, posY, posZ, radius} [normalized]
	float4 color{};					// {r, g, b, a} [0-1]	
	uint4 flags{};

	void HighLight(bool highLight) {
		flags.x = highLight ? 1 : 0;
	}

	bool IsDisabled() const { return position.x == std::numeric_limits<float>::max() && position.y == std::numeric_limits<float>::max() && position.z == std::numeric_limits<float>::max(); }
	__device__ __host__ static constexpr float4 Disabled() { return float4{ std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() }; }
};
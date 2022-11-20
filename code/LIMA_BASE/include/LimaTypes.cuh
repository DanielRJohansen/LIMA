#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"



#include "math.h"
#include <iostream>

#include <string>
#include <fstream>
#include <sstream>


struct TestType {
	static int get();
};


constexpr float PI = 3.14159f;

enum ATOM_TYPE { NONE, O, C, P, N, H, SOL, S };


struct Int3 {
	__host__ __device__ Int3() {}
	__host__ __device__ Int3(int x, int y, int z) : x(x), y(y), z(z) {}


	__host__ __device__ inline Int3 operator + (const Int3 a) const { return Int3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Int3 operator - (const Int3 a) const { return Int3(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ inline Int3 operator * (const int a) const { return Int3(x * a, y * a, z * a); }
	__host__ __device__ inline Int3 operator * (const float a) const { return Int3((int)floor((float)x * a), (int)floor((float)y * a), (int)floor((float)z * a)); }


	int x = 0, y = 0, z = 0;
};



struct Float3 {
	__host__ __device__ Float3() {}
	__host__ __device__ Float3(float a) : x(a), y(a), z(a) {}
	__host__ __device__ Float3(float x, float y, float z) : x(x), y(y), z(z) {}
	__host__ __device__ Float3(float* a) { x = a[0]; y = a[1]; z = a[2]; }

	__host__ __device__ inline Float3 operator * (const float a) const { return Float3(x * a, y * a, z * a); }
	//__host__ __device__ inline Float3 operator * (const double a) const { return Float3((float) (x * a), (float) (y * a), (float) (z * a)); }
	__host__ __device__ inline Float3 operator * (const Float3 a) const { return Float3(x * a.x, y * a.y, z * a.z); }
	__host__ __device__ inline Float3 operator + (const Float3 a) const { return Float3(x + a.x, y + a.y, z + a.z); }
	__host__ __device__ inline Float3 operator - (const Float3 a) const { return Float3(x - a.x, y - a.y, z - a.z); }
	__host__ __device__ inline bool operator == (const Float3 a) const { return (a.x == x && a.y == y && a.z == z); }
	__host__ __device__ inline void operator += (const Float3 a) { x += a.x; y += a.y; z += a.z; }
	__host__ __device__ inline void operator -= (const Float3 a) { x -= a.x; y -= a.y; z -= a.z; }
	__host__ __device__ inline void operator *= (const float a) { x *= a; y *= a; z *= a; }

	__host__ __device__ inline bool operator < (const Float3 a) { return x < a.x&& y < a.y&& z < a.z; }
	__host__ __device__ inline bool operator > (const Float3 a) { return x > a.x && y > a.y && z > a.z; }

	__host__ __device__ Float3 norm() {
		float l = len();
		if (l)
			return *this * (1.f / l);
		return Float3(0, 0, 0);
	}
	__device__ Float3 norm_fast() {		// Unsafe, may divide by 0
		return *this * (1.f / len());
	}
	__host__ __device__ Float3 square() { return Float3(x * x, y * y, z * z); }
	__host__ __device__ inline float len() { return (float)sqrtf(x * x + y * y + z * z); }
	__host__ __device__ inline float lenSquared() { return (float)(x * x + y * y + z * z); }
	__host__ __device__ Float3 zeroIfAbove(float a) { return Float3(x * (x < a), y * (y < a), z * (z < a)); }
	__host__ __device__ Float3 zeroIfBelow(float a) { return Float3(x * (x > a), y * (y > a), z * (z > a)); }
	__host__ __device__ Float3 elementwiseModulus(float a) {
		while (x > a)
			x -= a;
		while (y > a)
			y -= a;
		while (z > a)
			z -= a;
		return *this;
	}

	__host__ __device__ inline static float getAngle(Float3 v1, Float3 v2) {
		float val = (v1.dot(v2)) / (v1.len() * v2.len());	// If i make this float, we get values over 1, even with the statements below! :(
		//if (val > 1.f || val < -1.f) { printf("Val1 %f !!\n", val);}
		val = val > 1.f ? 1.f : val;
		val = val < -1.f ? -1.f : val;
		if (val > 1.f || val < -1.f) {
			printf("Val2 %f !!\n", val);
		}
		return acos(val);
	}
	__host__ __device__ static float getAngle(Float3 a, Float3 middle, Float3 b) {
		return getAngle(a - middle, b - middle);
	}


	__host__ __device__ Float3 cross(Float3 a) const { return Float3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }
	__host__ __device__ float dot(Float3 a) const { return (x * a.x + y * a.y + z * a.z); }
	__host__ __device__ Float3 abs() const {
		return Float3(
			std::abs(x),
			std::abs(y),
			std::abs(z)
		);
	}


	__host__ __device__ void print(char c = '_') {
		if (len() < 100000)
			printf("%c %f %f %f\n", c, x, y, z);
		else
			printf("%c %.0f\t %.0f\t %.0f\n", c, x, y, z);
	}



	__host__ __device__ void rotateAroundOrigo(Float3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// pitch and yaw is relative to global coordinates. 

		*this = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);
		*this = rodriguesRotatation(*this, Float3(0, 1, 0), pitch_yaw_roll.y);
		*this = rodriguesRotatation(*this, Float3(0, 0, 1), pitch_yaw_roll.y);
	}

	// I think the one below is wrong...
	__host__ __device__ Float3 _rotateAroundOrigin(Float3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// pitch and yaw is relative to global coordinates. Tilt is relative to body direction

		Float3 v = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);

		// Yaw around z
		v = rodriguesRotatation(v, Float3(0, 0, 1), pitch_yaw_roll.y);

		return v;
	}
	__host__ __device__ Float3 rotateAroundVector(Float3 pitch_yaw_roll, Float3 k) {	// k=normal = z-pointing		
		Float3 v = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);

		v = rodriguesRotatation(v, Float3(0, 0, 1), pitch_yaw_roll.y);
		v = rodriguesRotatation(v, k, pitch_yaw_roll.z);
		return v;
	}
	__host__ __device__ static Float3 rodriguesRotatation(Float3 v, Float3 k, float theta) {
		return v * cos(theta) + k.cross(v) * sin(theta) + k * (k.dot(v)) * (1 - cos(theta));
	}

	__host__ __device__ float at(int index) const{
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

	__host__ __device__ float* placeAt(int index) {
		switch (index) {
		default:			// Sadly it is not reasonable to catch this in release
		case 0:
			return &x;
		case 1:
			return &y;
		case 2:
			return &z;
		}
	}

	// Not used right now!
	__host__ __device__ static Float3 centerOfMass(Float3* arr_ptr, uint32_t arr_size) {	// Only run before sim, so we can cast to double without slowing sim
		Float3 sum = Float3(0, 0, 0);
		for (uint32_t i = 0; i < arr_size; i++) {
			sum = sum + arr_ptr[i];
		}
		return sum * (1.f / arr_size);
	}



	float x = 0, y = 0, z = 0;
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
	__host__ __device__ inline void operator += (const Float3 a) { x += (double)a.x; y += (double)a.y; z += (double)a.z; }
	__host__ __device__ inline void operator += (const Double3 a) { x += a.x; y += a.y; z += a.z; }


	__host__ __device__ inline double len() { return (double)sqrt(x * x + y * y + z * z); }

	__host__ __device__ void print(char c = '_') {
		printf("%c %f %f %f\n", c, x, y, z);
	}

	double x = 0, y = 0, z = 0;
};

struct BoundingBox {
	BoundingBox() {}
	BoundingBox(Float3 min, Float3 max) : min(min), max(max) {}


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

// TODO make this
struct CompactBool {

};

// This struct is agnostic as to which compounds the particles in the table belongs to!!
template <typename T, int len>
class FixedSizeMatrix {
public:
	__device__ FixedSizeMatrix() {}	// Fast constructor for device
	__host__ FixedSizeMatrix(T val) {
		for (int i = 0; i < m_len * m_len; i++) { matrix[i] = val; }
	}
	__host__ FixedSizeMatrix(FixedSizeMatrix<T, len>& src) { *this = src; }

	//__host__ init()

	__host__ __device__ T* get(int i1, int i2) { return &matrix[i1 + i2 * m_len]; }

	__host__ void set(int i1, int i2, T object) { matrix[i1 + i2 * m_len] = object; }

	__device__ void load(FixedSizeMatrix<T, len>& src) {
		for (int i = threadIdx.x; i < m_size; i += blockDim.x) {
			matrix[i] = src.get(i);
		}
	}

private:		
	const static int m_len = len;
	const static int m_size = m_len * m_len;
	T matrix[m_len * m_len];

	// To be used for loading only
	__device__ T get(int index) { return matrix[index]; }
};









template<typename T>
T* genericMoveToDevice(T* data_ptr, int n_elements) {	// Currently uses MallocManaged, switch to unmanaged for safer operation
	T* gpu_ptr = nullptr;
	int bytesize = n_elements * sizeof(T);

	cudaMallocManaged(&gpu_ptr, bytesize);
	cudaMemcpy(gpu_ptr, data_ptr, bytesize, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	delete[] data_ptr;

	data_ptr = gpu_ptr;

	//printf("Moved %.2f MB to device\n", bytesize*1e-6);
	return gpu_ptr;
}




struct Trajectory {
	Trajectory() {}
	Trajectory(std::string path) {
		positions = new Float3[max_size];
		float buffer[3] = {};

		printf("Path: ");
		std::cout << path << std::endl;

		std::string file_contents = readFileIntoString(path);
		std::istringstream sstream(file_contents);
		std::string record;

		int counter = 0;

		int step_cnt = 0;
		while (std::getline(sstream, record)) {
			std::istringstream line(record);

			int dim = 0;
			while (getline(line, record, ';')) {
				buffer[dim++] = stof(record);

				if (dim == 3) {
					positions[counter++] = Float3(buffer);
					dim = 0;
					//positions[counter - 1].print();
				}
			}
			if (step_cnt == 0) {
				n_particles = counter;
			}
			step_cnt++;
		}

		particle_type = new int[n_particles];
		n_steps = step_cnt;
		printf("Loaded trajectory with %d particles and %d steps\n", n_particles, n_steps);
	}
	Trajectory(int size) {
		positions = new Float3[size];
	}

	void moveToDevice() {
		positions = genericMoveToDevice(positions, max_size);
		particle_type = genericMoveToDevice(particle_type, n_particles);
		cudaDeviceSynchronize();
		printf("Success\n");
	}



	std::string readFileIntoString(const std::string& path) {
		auto ss = std::ostringstream{};
		std::ifstream input_file(path);
		if (!input_file.is_open()) {
			std::cerr << "Could not open the file - '"
				<< path << "'" << std::endl;
			exit(EXIT_FAILURE);
		}
		ss << input_file.rdbuf();
		return ss.str();
	}

	int max_size = 100 * 10000;

	Float3* positions = nullptr;
	int* particle_type = nullptr;	//0 solvent, 1 compound, 2 virtual

	int n_particles = 0;
	int n_steps = 0;
};



class HashTable {
public:
	HashTable(uint16_t* keys, int n_keys, int ts) : table_size(ts) {
		table = new int[table_size]();
		for (int i = 0; i < table_size; i++) {
			table[i] = -1;
		}

		for (int i = 0; i < n_keys; i++)
			insert(keys[i]);
	}

	bool insert(uint16_t key, int offset = 1) {			// Returns true for sucessful insertion
		if (offset > 100000) {
			printf("Hashtable insertion failed\n");
			exit(1);
		}


		uint32_t hash = (getHash(key) + offset) % table_size;
		if (table[hash] == key) {		// Key already exists in table
			return false;
		}
		else if (table[hash] == -1) {	// Key doesn't exist in table
			table[hash] = key;
			return true;
		}
		else {							// Hash already taken, recurse
			return insert(key, offset * 2);
		}
	}
	~HashTable() {
		delete[] table;
	}

private:

	uint32_t getHash(uint16_t key) {
		return (uint32_t)floor(table_size * (fmod((double)key * k, 1.)));
	}


	int* table = nullptr;
	int table_size = 0;
	double k = 0.79026;
};




struct RenderBall {
	RenderBall() {}
	__host__ __device__ RenderBall(Float3 pos, float radius, Int3 color) :pos(pos), radius(radius), color(color) {}
	Float3 pos;	// only uses x and y
	float radius = 0.f;
	Int3 color;
	bool disable = false;
};

enum VerbosityLevel {
	SILENT,
	CRITICAL_INFO, 
	V1,
	V2,
	V3
};
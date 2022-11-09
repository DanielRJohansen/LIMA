#include "QuantomTypes.cuh"


//__host__ __device__ inline double Float3::len() { return (double) sqrtf(x * x + y * y + z * z); }

/*
__host__ __device__ Float3 rotateAroundOrigin(Float3 pitch_yaw_roll) {	//pitch around x, yaw around z, tilt around y
		// Pitch around x
	printf("\n\n\n Rotating %f %f %f\n", x, y, z);
	Float3 v = rodriguesRotatation(*this, Float3(1, 0, 0), pitch_yaw_roll.x);

	printf("After x: %f %f %f\n", v.x, v.y, v.z);

	// Yaw around y
	v = rodriguesRotatation(v, Float3(0, 1, 0), pitch_yaw_roll.y);
	printf("After y: %f %f %f\n", v.x, v.y, v.z);

	// Tilt around itself
	v = rodriguesRotatation(v, v, pitch_yaw_roll.z);
	printf("After self: %f %f %f\n", v.x, v.y, v.z);
	return v;


	/*double sin_pitch = sin(pitch_yaw_roll.x);
	double cos_pitch = cos(pitch_yaw_roll.x);

	double sin_yaw = sin(pitch_yaw_roll.z);
	double cos_yaw = cos(pitch_yaw_roll.z);

	double sin_tilt = sin(pitch_yaw_roll.y);
	double cos_tilt = cos(pitch_yaw_roll.y);

	//Rotate around x
	double y_x = cos_pitch * y + sin_pitch * z;
	double z_x = -sin_pitch * y + cos_pitch * z;

	// Rotate 

}*/
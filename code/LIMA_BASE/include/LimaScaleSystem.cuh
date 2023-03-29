#pragma once

#include <assert.h>

namespace LIMA_UTILS {

	/*template <typename T>
	assertEqual(T a, T b) */
	int roundUp(int numToRound, int multiple)
	{
		assert(multiple);
		return ((numToRound + multiple - 1) / multiple) * multiple;
	}
}
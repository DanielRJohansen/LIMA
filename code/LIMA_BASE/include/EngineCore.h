#pragma once

#include "BoundaryConditionPublic.h"


namespace LIMAPOSITIONSYSTEM {
	static float calcHyperDistNM(const Float3& p1, const Float3& p2, float boxlen_nm, BoundaryConditionSelect bc) {
		Float3 temp = p2;
		BoundaryConditionPublic::applyHyperposNM(p1, temp, boxlen_nm, bc);
		return (p1 - temp).len();
	}
}
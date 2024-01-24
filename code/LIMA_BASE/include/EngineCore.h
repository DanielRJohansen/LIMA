#pragma once

#include "BoundaryConditionPublic.h"


namespace LIMAPOSITIONSYSTEM {
	static float calcHyperDistNM(const Float3* const p1, const Float3* const p2, float boxlen_nm, BoundaryConditionSelect bc) {
		Float3 temp = *p2;
		//LIMAPOSITIONSYSTEM::applyHyperposNM<BoundaryCondition>(p1, &temp);
		BoundaryConditionPublic::applyHyperposNM(p1, &temp, boxlen_nm, bc);
		return (*p1 - temp).len();
	}
}
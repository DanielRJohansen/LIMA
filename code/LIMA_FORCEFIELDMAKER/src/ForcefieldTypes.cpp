#include "ForcefieldTypes.h"




using std::cout, std::endl, std::string, std::vector;



//void Singlebondtype::sort() {
//	if (!FTHelpers::isSorted(&bonded_typenames[0], &bonded_typenames[1])) {
//		swap(bonded_typenames[0], bonded_typenames[1]);
//		//std::swap(id1, id2);	// TODO: Should this not be like this???
//	}
//}

//void Anglebondtype::sort() {
//	// The middle type must always be in the middle, so we only sort the outer 2 types
//	// No need to sort if our edges are the same
//	if (bonded_typenames[0] != bonded_typenames[2]) {
//		if (!FTHelpers::isSorted(&bonded_typenames[0], &bonded_typenames[2])) {
//			flip()
//		}
//	}
//}

//void Dihedralbondtype::sort() {
//	if (bonded_typenames[0] != bonded_typenames[3]) {
//		if (!FTHelpers::isSorted(&bonded_typenames[0], &bonded_typenames[3])) {
//			flip();
//		}
//	}
//	else {			// In case the outer two is identical, we check the inner two.
//		if (!FTHelpers::isSorted(&bonded_typenames[1], &bonded_typenames[2])) {
//			flip();
//		}
//	}
//}

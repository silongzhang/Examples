#include"DataStructure.h"


// Whether data of this instance are consistent.
bool Instance::consistent() const {
	const int NIt = varInt.size, NCt = varCont.size;
	if (NIt == 0 || NCt == 0) return false;

	// Whether varInt is consistent.
	if (varInt.type != VariableType::Integer) return false;
	if (varInt.upperbounds.size() != NIt) return false;
	for (const auto num : varInt.upperbounds)
		if (lessThanReal(num, 0, PPM) && !equalToReal(num, INFUB, PPM))
			return false;

	// Whether varCont is consistent.
	if (varCont.type != VariableType::Continuous) return false;
	if (varCont.upperbounds.size() != NCt) return false;
	for(const auto num: varCont.upperbounds)
		if (lessThanReal(num, 0, PPM) && !equalToReal(num, INFUB, PPM))
			return false;

	// Whether obj is consistent.
	if (obj.coefInt.size() != NIt || obj.coefCont.size() != NCt) return false;

	// Whether cpCons is consistent.
	if (cpCons.empty()) return false;
	for (const auto& cons : cpCons) {
		if (cons.coefInt.size() != NIt || cons.coefCont.size() != NCt) return false;

		// There should be some nonzero(s) in cons.coefCont.
		int i = 0;
		for (; i < NCt; ++i) {
			if (!equalToReal(cons.coefCont[i], 0, PPM))
				break;
		}
		if (i == NCt) return false;
	}

	// Whether itCons is consistent.
	for (const auto& cons : itCons) {
		if (cons.coefInt.size() != NIt || cons.coefCont.size() != NCt) return false;

		// All elements in cons.coefCont must be zero.
		int i = 0;
		for (; i < NCt; ++i) {
			if (!equalToReal(cons.coefCont[i], 0, PPM))
				break;
		}
		if (i != NCt) return false;
	}

	return true;
}


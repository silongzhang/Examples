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


// Whether (i) it is consistent, (ii) it is a minimization problem and (iii) for any continuous variable y, y >= 0 is specified in the domain whereas y >= a and y <= b are specified in constraints.
bool Instance::standard() const {
	if (!consistent()) return false;

	if (obj.type == ObjectiveType::Maximization) return false;

	for (const auto elem : varCont.upperbounds)
		if (!equalToReal(elem, INFUB, PPM))
			return false;

	return true;
}


void Instance::standardize() {
	try {
		if (!consistent()) throw exception();

		if (obj.type == ObjectiveType::Maximization) {
			obj.type == ObjectiveType::Minimization;
			for (auto& num : obj.coefInt) num = -num;
			for (auto& num : obj.coefCont) num = -num;
			cout << "*****************************" << endl;
			cout << "The maximization problem is transformed into a minimization problem." << endl;
			cout << "*****************************" << endl;
		}

		const int NInt = varInt.size, NCont = varCont.size;
		vector<double> cfIt(NInt, 0);
		for (int i = 0; i < varCont.upperbounds.size(); ++i) {
			if (!equalToReal(varCont.upperbounds[i], INFUB, PPM)) {
				vector<double> cfCt(NCont, 0);
				cfCt[i] = 1;
				cpCons.push_back(Constraint(ConstraintType::Le, cfIt, cfCt, varCont.upperbounds[i]));

				varCont.upperbounds[i] = INFUB;
			}
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::standardize", exc);
	}
}


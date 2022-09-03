#include"DataStructure.h"


bool ParameterTest::valid() const {
	if (NInt < 1 || NCont < 1) return false;
	if (lessThanReal(MaxUBInt, 1, PPM) || lessThanReal(MaxUBCont, 1, PPM)) return false;
	if (lessThanReal(ratioInfUBInt, 0, PPM) || greaterThanReal(ratioInfUBInt, 1, PPM)) return false;
	if (lessThanReal(ratioInfUBCont, 0, PPM) || greaterThanReal(ratioInfUBCont, 1, PPM)) return false;
	if (!lessThanReal(rangeObjCoefInt.first, rangeObjCoefInt.second - 2, PPM)) return false;
	if (!lessThanReal(rangeObjCoefCont.first, rangeObjCoefCont.second - 2, PPM)) return false;
	if (!lessThanReal(rangeConsCoefInt.first, rangeConsCoefInt.second - 2, PPM)) return false;
	if (!lessThanReal(rangeConsCoefCont.first, rangeConsCoefCont.second - 2, PPM)) return false;
	if (lessThanReal(sparsity, 0, PPM) || !lessThanReal(sparsity, 1, PPM)) return false;
	if (lessThanReal(ratioEq, 0, PPM) || greaterThanReal(ratioEq, 1, PPM)) return false;
	if (lessThanReal(ratioLe, 0, PPM) || greaterThanReal(ratioLe, 1, PPM)) return false;
	if (lessThanReal(ratioGe, 0, PPM) || greaterThanReal(ratioGe, 1, PPM)) return false;
	if (!equalToReal(ratioEq + ratioLe + ratioGe, 1, PPM)) return false;
	if (NCpCons < 1) return false;

	return true;
}


vector<double> generateRandomVector(default_random_engine& engine, int N, pair<double, double> range, double sparsity, bool canEmpty) {
	vector<double> result;

	uniform_real_distribution<double> unifZeroOne(0, 1);
	uniform_int_distribution<int> unifRange(range.first, range.second);
	bool flag = !canEmpty;

	do {
		result = vector<double>(N, 0);
		for (int i = 0; i < N; ++i) {
			if (greaterThanReal(unifZeroOne(engine), sparsity, PPM)) {
				result[i] = unifRange(engine);
				if (!equalToReal(result[i], 0, PPM)) flag = false;
			}
		}
	} while (flag);

	return result;
}


void setRandomSignAndRhs(default_random_engine& engine, const vector<double>& midInt, const vector<double>& midCont, Constraint& cons, double ratioEq, double ratioLe, double ratioGe, double redundancy) {
	uniform_real_distribution<double> unifZeroOne(0, 1);
	try {
		if (midInt.size() != cons.coefInt.size() || midCont.size() != cons.coefCont.size()) throw exception();

		double expr = 0;
		for (int i = 0; i < cons.coefInt.size(); ++i) expr += midInt[i] * cons.coefInt[i];
		for (int i = 0; i < cons.coefCont.size(); ++i) expr += midCont[i] * cons.coefCont[i];

		const double sign = unifZeroOne(engine);
		if (lessThanReal(sign, ratioEq, PPM)) {
			cons.type = ConstraintType::Eq;
			cons.rhs = expr;
		}
		else if (lessThanReal(sign, ratioEq + ratioLe, PPM)) {
			cons.type = ConstraintType::Le;
			cons.rhs = expr + abs(expr) * redundancy;
		}
		else {
			cons.type = ConstraintType::Ge;
			cons.rhs = expr - abs(expr) * redundancy;
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("setRandomSignAndRhs", exc);
	}
}


Instance ParameterTest::generateInstance() const {
	Instance result;
	try {
		if (!valid()) throw exception();

		default_random_engine engine(seed);
		uniform_real_distribution<double> unifZeroOne(0, 1);

		// Set varInt and varCont.
		result.varInt.type = VariableType::Integer, result.varCont.type = VariableType::Continuous;
		result.varInt.size = NInt, result.varCont.size = NCont;
		result.varInt.upperbounds.resize(NInt), result.varCont.upperbounds.resize(NCont);

		uniform_int_distribution<int> unifUBInt(1, MaxUBInt);
		for (auto& num : result.varInt.upperbounds)
			num = lessThanReal(unifZeroOne(engine), ratioInfUBInt, PPM) ? INFUB : unifUBInt(engine);

		uniform_int_distribution<int> unifUBCont(1, MaxUBCont);
		for (auto& num : result.varCont.upperbounds)
			num = lessThanReal(unifZeroOne(engine), ratioInfUBCont, PPM) ? INFUB : unifUBCont(engine);

		// Set obj.
		result.obj.type = minimization ? ObjectiveType::Minimization : ObjectiveType::Maximization;
		result.obj.coefInt.resize(NInt), result.obj.coefCont.resize(NCont);

		uniform_int_distribution<int> unifObjInt(rangeObjCoefInt.first, rangeObjCoefInt.second);
		for (auto& num : result.obj.coefInt)
			num = unifObjInt(engine);

		uniform_int_distribution<int> unifObjCont(rangeObjCoefCont.first, rangeObjCoefCont.second);
		for (auto& num : result.obj.coefCont)
			num = unifObjCont(engine);

		// Treat the midpoint as a feasible solution.
		vector<double> midInt(NInt), midCont(NCont);
		for (int i = 0; i < NInt; ++i) midInt[i] = floor(result.varInt.upperbounds[i] / 2);
		for (int i = 0; i < NCont; ++i) midCont[i] = floor(result.varCont.upperbounds[i] / 2);

		// Set cpCons and itCons.
		uniform_int_distribution<int> unifConsInt(rangeConsCoefInt.first, rangeConsCoefInt.second);
		uniform_int_distribution<int> unifConsCont(rangeConsCoefCont.first, rangeConsCoefCont.second);

		result.cpCons.resize(NCpCons), result.itCons.resize(NItCons);
		for (auto& cons : result.cpCons) {
			cons.coefInt = generateRandomVector(engine, NInt, rangeConsCoefInt, sparsity, true);
			cons.coefCont = generateRandomVector(engine, NCont, rangeConsCoefCont, sparsity, false);
			setRandomSignAndRhs(engine, midInt, midCont, cons, ratioEq, ratioLe, ratioGe, redundancy);
		}
		for (auto& cons : result.itCons) {
			cons.coefInt = generateRandomVector(engine, NInt, rangeConsCoefInt, sparsity, false);
			cons.coefCont = vector<double>(NCont, 0);
			setRandomSignAndRhs(engine, midInt, midCont, cons, ratioEq, ratioLe, ratioGe, redundancy);
		}

		if (!result.consistent()) throw exception();
	}
	catch (const exception& exc) {
		printErrorAndExit("ParameterTest::generateInstance", exc);
	}

	return result;
}


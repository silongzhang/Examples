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
			obj.type = ObjectiveType::Minimization;
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


void Solution::clear() {
	status = SolutionStatus::Unkown;
	objective = InfinityPos;
	LB = InfinityNeg;
	valueInt.clear();
	valueEta = InfinityNeg;
	valueCont.clear();

	nOptCutLP = nOptCutIP = nFeasCutLP = nFeasCutIP = 0;
	elapsedTime = 0;
}


IloExpr product(IloEnv env, const vector<double>& coefs, IloNumVarArray vars) {
	IloExpr expr(env);
	try {
		if (coefs.size() != vars.getSize() || coefs.empty()) throw exception();

		for (int i = 0; i < coefs.size(); ++i)
			expr += coefs[i] * vars[i];
	}
	catch (const exception& exc) {
		printErrorAndExit("product", exc);
	}
	return expr;
}


IloExpr product(IloEnv env, const vector<double>& coefs, IloIntVarArray vars) {
	IloExpr expr(env);
	try {
		if (coefs.size() != vars.getSize() || coefs.empty()) throw exception();

		for (int i = 0; i < coefs.size(); ++i)
			expr += coefs[i] * vars[i];
	}
	catch (const exception& exc) {
		printErrorAndExit("product", exc);
	}
	return expr;
}


void addConstraint(IloModel model, const vector<double>& coefs, IloNumVarArray vars, ConstraintType type, double rhs) {
	try {
		if (coefs.size() != vars.getSize() || coefs.empty()) throw exception();

		IloExpr expr = product(model.getEnv(), coefs, vars);
		if (type == ConstraintType::Eq) model.add(expr == rhs);
		else if (type == ConstraintType::Le) model.add(expr <= rhs);
		else model.add(expr >= rhs);
	}
	catch (const exception& exc) {
		printErrorAndExit("addConstraint", exc);
	}
}


IloRange genCons(IloEnv env, const vector<double>& coefs, IloNumVarArray vars, ConstraintType type, double rhs) {
	IloRange result;
	try {
		if (coefs.size() != vars.getSize() || coefs.empty()) throw exception();

		IloExpr expr = product(env, coefs, vars);
		if (type == ConstraintType::Eq) result = IloRange(env, rhs, expr, rhs);
		else if (type == ConstraintType::Le) result = IloRange(env, -IloInfinity, expr, rhs);
		else result = IloRange(env, rhs, expr, IloInfinity);
	}
	catch (const exception& exc) {
		printErrorAndExit("genCons", exc);
	}
	return result;
}


vector<double> getValues(IloCplex cplex, IloNumVarArray vars) {
	vector<double> result;
	try {
		result.resize(vars.getSize());
		for (int i = 0; i < vars.getSize(); ++i) result[i] = cplex.getValue(vars[i]);
	}
	catch (const exception& exc) {
		printErrorAndExit("getValues", exc);
	}
	return result;
}


IloExpr Instance::exprRhs(IloEnv env, IloNumArray duals, IloNumVarArray vars) const {
	IloExpr expr(env);
	try {
		if (duals.getSize() != cpCons.size() || vars.getSize() != varInt.size) throw exception();

		for (int i = 0; i < cpCons.size(); ++i) {
			expr += duals[i] * (cpCons[i].rhs - product(env, cpCons[i].coefInt, vars));
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::exprRhs", exc);
	}
	return expr;
}


void setRhs(IloRange constraint, ConstraintType type, double rhs) {
	if (type == ConstraintType::Eq) constraint.setBounds(rhs, rhs);
	else if (type == ConstraintType::Le) constraint.setUB(rhs);
	else constraint.setLB(rhs);
}


void Solution::renew(IloCplex cplexRMP, IloNumVarArray X, IloNumVar eta, IloCplex cplexSP, IloNumVarArray Y) {
	try {
		objective = cplexRMP.getObjValue() - cplexRMP.getValue(eta) + cplexSP.getObjValue();
		LB = cplexRMP.getObjValue();
		valueInt = getValues(cplexRMP, X);
		valueEta = cplexRMP.getValue(eta);
		valueCont = getValues(cplexSP, Y);

		status = SolutionStatus::Feasible;
		for (const auto num : valueInt) {
			if (!isInteger(num, PPM)) {
				status = SolutionStatus::Infeasible;
				break;
			}
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("Solution::renew", exc);
	}
}


void Solution::print() const {
	cout << "status = ";
	::print(status);
	cout << '\t' << "objective = " << objective << '\t' << "LB = " << LB << '\t' << "nOptCutLP = " << nOptCutLP << '\t' << "nOptCutIP = " << nOptCutIP << '\t' << "nFeasCutLP = " << nFeasCutLP << '\t' << "nFeasCutIP = " << nFeasCutIP << '\t' << "elapsedTime = " << elapsedTime << endl;
}


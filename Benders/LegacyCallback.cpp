#include"DataStructure.h"

// Implement the benders decomposition with legacy callback.


void Instance::initiateModels(IloEnv env, IloModel modelRMP, IloModel modelSP, IloNumVarArray X, IloNumVarArray Y, IloNumVar eta, const vector<double>& currentValInt, IloRangeArray consSP) const {
	try {
		// Define variables.
		for (int i = 0; i < varInt.size; ++i)
			X.add(IloNumVar(env, 0, equalToReal(varInt.upperbounds[i], INFUB, PPM) ? IloInfinity : varInt.upperbounds[i]));

		// Set the objective functions.
		IloExpr expr = product(env, obj.coefInt, X);
		expr += eta;
		modelRMP.add(IloMinimize(env, expr));

		expr = product(env, obj.coefCont, Y);
		modelSP.add(IloMinimize(env, expr));

		// Set constraints other than Benders cuts in the RMP.
		for (const auto& cons : itCons) addConstraint(modelRMP, cons.coefInt, X, cons.type, cons.rhs);
		modelRMP.add(eta >= InfinityNeg);						// Avoid the RMP to be unbounded.

		// Set constraints in the SP.
		for (const auto& cons : cpCons) {
			double rhs = cons.rhs - inner_product(cons.coefInt.begin(), cons.coefInt.end(), currentValInt.begin(), 0.0);
			consSP.add(genCons(env, cons.coefCont, Y, cons.type, rhs));
		}
		modelSP.add(consSP);
	}
	catch (const exception& exc) {
		printErrorAndExit("initiateModels", exc);
	}
}


bool solveLRMP(IloCplex cplexRMP, IloNumVarArray X, IloNumVar eta, Solution& incumbent, vector<double>& currentValInt, double& currentValEta) {
	try {
		solveModel(cplexRMP);
		if (cplexRMP.getStatus() == IloAlgorithm::Status::Infeasible) {
			cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cout << "The instance is infeasible! Terminate!" << endl;
			incumbent.status = SolutionStatus::Infeasible;
			return false;
		}
		else if (cplexRMP.getStatus() == IloAlgorithm::Status::Unbounded) {
			cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cout << "The current RMP might be unbounded (Assumption 6 is not satisfied)! Terminate!" << endl;
			incumbent.status = SolutionStatus::Unkown;
			return false;
		}
		else if (cplexRMP.getStatus() == IloAlgorithm::Status::Optimal) {
			currentValInt = getValues(cplexRMP, X);
			currentValEta = cplexRMP.getValue(eta);
			incumbent.LB = cplexRMP.getObjValue();			// Renew LB.
			return true;
		}
		else throw exception();
	}
	catch (const exception& exc) {
		printErrorAndExit("solveLRMP", exc);
	}
	return false;
}


bool solveSP(IloCplex cplexSP, IloRangeArray consSP, IloNumArray& dualSP, Solution& incumbent, const vector<Constraint>& cpCons, const vector<double>& currentValInt) {
	try {
		// Set the right hand side of constraints in the SP.
		for (int i = 0; i < consSP.getSize(); ++i) {
			double rhs = cpCons[i].rhs - inner_product(cpCons[i].coefInt.begin(), cpCons[i].coefInt.end(), currentValInt.begin(), 0.0);
			setRhs(consSP[i], cpCons[i].type, rhs);
		}

		solveModel(cplexSP);								// Solve the subproblem.
		if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible || cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
			dualSP = getDuals(cplexSP, consSP);				// Get dual values.
			return true;
		}
		else {
			if (cplexSP.getStatus() == IloAlgorithm::Status::Unbounded) {
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				cout << "The instance might be unbounded! Terminate!" << endl;
			}
			incumbent.status = SolutionStatus::Unkown;
			return false;
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("solveSP", exc);
	}
	return false;
}


bool addBendersCuts(IloEnv env, IloCplex cplexSP, IloModel modelRMP, IloNumVarArray X, IloNumVar eta, IloNumArray dualSP, Solution& incumbent, const Instance& instance, double currentValEta, bool integral) {
	try {
		if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible) {
			modelRMP.add(0 >= instance.exprRhs(env, dualSP, X));			// Add feasibility cut.
			integral ? ++incumbent.nFeasCutIP : ++incumbent.nFeasCutLP;
		}
		else if (cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
			if (lessThanReal(currentValEta, cplexSP.getObjValue(), PPM)) {
				modelRMP.add(eta >= instance.exprRhs(env, dualSP, X));		// Add optimality cut.
				integral ? ++incumbent.nOptCutIP : ++incumbent.nOptCutLP;
			}
			else if (equalToReal(currentValEta, cplexSP.getObjValue(), PPM)) return false;
			else throw exception();
		}
		else throw exception();
	}
	catch (const exception& exc) {
		printErrorAndExit("addBendersCuts", exc);
	}
	return true;
}


ILOLAZYCONSTRAINTCALLBACK7(LazyBendersCuts, IloNumVarArray, X, IloNumVarArray, Y, IloNumVar, eta, IloCplex, cplexSP, IloRangeArray, consSP, Solution&, incumbent, const Instance&, instance) {
	try {
		incumbent.LB = getBestObjValue();														// Renew LB.
		if (greaterThanReal(getObjValue(), getIncumbentObjValue(), PPM)) return;

		IloEnv env = getEnv();
		IloNumArray valueX(env), dualSP(env);
		getValues(valueX, X);
		const vector<double> currentValInt = assign(valueX);
		const double currentValEta = getValue(eta);

		if (!solveSP(cplexSP, consSP, dualSP, incumbent, instance.cpCons, currentValInt)) throw exception();
		if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible)
			add(0 >= instance.exprRhs(env, dualSP, X)), ++incumbent.nFeasCutIP;					// Add feasibility cut.
		else if (cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
			if (greaterThanReal(incumbent.objective, getObjValue() - getValue(eta) + cplexSP.getObjValue(), PPM)) {
				incumbent.objective = getObjValue() - getValue(eta) + cplexSP.getObjValue();	// Renew UB.
				incumbent.valueInt = currentValInt;
				incumbent.valueEta = currentValEta;
				incumbent.valueCont = ::getValues(cplexSP, Y);
				incumbent.status = SolutionStatus::Feasible;
			}

			if (lessThanReal(currentValEta, cplexSP.getObjValue(), PPM))
				add(eta >= instance.exprRhs(env, dualSP, X)), ++incumbent.nOptCutIP;			// Add optimality cut.
			else if (greaterThanReal(currentValEta, cplexSP.getObjValue(), PPM)) throw exception();

			cout << "UB = " << incumbent.objective << '\t' << "LB = " << incumbent.LB << '\t' << "nOptCutLP = " << incumbent.nOptCutLP << '\t' << "nOptCutIP = " << incumbent.nOptCutIP << '\t' << "nFeasCutLP = " << incumbent.nFeasCutLP << '\t' << "nFeasCutIP = " << incumbent.nFeasCutIP << endl;
		}
		else throw exception();
	}
	catch (const exception& exc) {
		printErrorAndExit("LazyBendersCuts", exc);
	}
}


ILOMIPINFOCALLBACK4(StoppingCriteria, bool&, aborted, const ParameterAlgorithm&, parameter, const Solution&, incumbent, clock_t, start) {
	try {
		if (!aborted && parameter.stop(start, incumbent.objective, incumbent.LB)) {
			aborted = true;
			abort();
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("StoppingCriteria", exc);
	}
}


// Assumption 1. The MILP model has a finite optimal objective value.
// Assumption 2. # of integer variables > 0 and # of continuous variables > 0.
// Assumption 3. Variables are already transformed to be nonnegative.
// Assumption 4. Invoke the method "standardize" in advance if the instance is not "standard".
// Assumption 5. It gives the opposite number of the optimal objective value of the original instance if it is a maximization problem.
// Assumption 6. Relaxed master problems (RMPs) cannot be unbounded.
// Note 1. For a minimization problem, a sufficient condition for Assumption 6 is that integer variables are all nonnegative and their coefficients in the objective function are all nonnegative. A similar sufficient condition can be derived for a maximization problem.
// Note 2. Another sufficient condition for Assumption 6 is that the number of feasible combinations of values of integer variables is finite. For example, each integer variable is bounded from both below and above.
Solution Instance::solveBendersLegacyCallback(const ParameterAlgorithm& parameter) const {
	Solution incumbent;
	incumbent.clear();
	IloEnv env;
	try {
		const clock_t start = clock();
		cout << "Running Instance::solveBendersLegacyCallback ..." << endl;
		if (!standard()) throw exception();
		const int NInt = varInt.size, NCont = varCont.size;

		// Declare modeling components.
		IloNumVar eta(env, -IloInfinity, IloInfinity);
		IloNumVarArray X(env), Y(env, NCont, 0, IloInfinity);
		IloModel modelRMP(env), modelSP(env);
		vector<double> currentValInt(NInt, 0);
		double currentValEta = InfinityNeg;
		IloRangeArray consSP(env);
		IloNumArray dualSP(env);

		// Initiate models.
		initiateModels(env, modelRMP, modelSP, X, Y, eta, currentValInt, consSP);
		IloCplex cplexRMP(modelRMP), cplexSP(modelSP);
		cplexSP.setParam(IloCplex::RootAlg, IloCplex::Primal);

		cout << "#############################" << endl;
		cout << "Add Benders cuts based on solutions of the LP relaxation of the RMP." << endl;
		bool proceed = true;
		for (int iter = 1; true; ++iter) {
			cout << "+++++++++++++++++++++++++++++" << endl;
			cout << "Iter = " << iter << '\t' << "LB = " << incumbent.LB << '\t' << "nOptCutLP = " << incumbent.nOptCutLP << '\t' << "nFeasCutLP = " << incumbent.nFeasCutLP << endl;

			// Solve the LP relaxation of the RMP.
			proceed = solveLRMP(cplexRMP, X, eta, incumbent, currentValInt, currentValEta);
			if (!proceed) break;

			// Solve the SP.
			proceed = solveSP(cplexSP, consSP, dualSP, incumbent, cpCons, currentValInt);
			if (!proceed) break;

			// Add Benders cuts.
			if (!addBendersCuts(env, cplexSP, modelRMP, X, eta, dualSP, incumbent, *this, currentValEta, false)) break;
		}

		cplexSP.setOut(env.getNullStream());
		bool aborted = false;
		if (proceed) {
			// Restore integrality constraints.
			modelRMP.add(IloConversion(env, X, IloNumVar::Type::Int));

			// Use callback to add lazy constraints.
			cplexRMP.use(LazyBendersCuts(env, X, Y, eta, cplexSP, consSP, incumbent, *this));

			// Use callback to stop the procedure according to stopping criteria.
			cplexRMP.use(StoppingCriteria(env, aborted, parameter, incumbent, start));

			solveModel(cplexRMP);
			incumbent.LB = cplexRMP.getBestObjValue();
		}

		if (incumbent.status == SolutionStatus::Feasible && !aborted) incumbent.status = SolutionStatus::Optimal;
		incumbent.elapsedTime = runTime(start);
		cout << "elapsed time (Instance::solveBendersLegacyCallback): " << runTime(start) << endl;
		cout << "# of nodes processed = " << cplexRMP.getNnodes() << '\t' << "# of nodes remained = " << cplexRMP.getNnodesLeft() << endl;
		incumbent.print();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersLegacyCallback", exc);
	}

	env.end();
	return incumbent;
}


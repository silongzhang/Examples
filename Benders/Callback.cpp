#include"DataStructure.h"

// Implement the benders decomposition with callback.


void Instance::initiateModels(IloEnv env, IloModel modelRMP, IloModel modelSP, IloNumVarArray X, IloNumVarArray Y, IloNumVar eta, const vector<double>& currentValInt, double currentValEta, IloRangeArray consSP) const {
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
	return true;
}


// Assumption 1. The MILP model has a finite optimal objective value.
// Assumption 2. # of integer variables > 0 and # of continuous variables > 0.
// Assumption 3. Variables are already transformed to be nonnegative.
// Assumption 4. Invoke the method "standardize" in advance if the instance is not "standard".
// Assumption 5. It gives the opposite number of the optimal objective value of the original instance if it is a maximization problem.
// Assumption 6. Relaxed master problems (RMPs) cannot be unbounded.
// Note 1. For a minimization problem, a sufficient condition for Assumption 6 is that integer variables are all nonnegative and their coefficients in the objective function are all nonnegative. A similar sufficient condition can be derived for a maximization problem.
// Note 2. Another sufficient condition for Assumption 6 is that the number of feasible combinations of values of integer variables is finite. For example, each integer variable is bounded from both below and above.
Solution Instance::solveBendersCallback(const ParameterAlgorithm& parameter) const {
	Solution incumbent;
	incumbent.clear();
	IloEnv env;
	try {
		const clock_t start = clock();
		cout << "Running Instance::solveBendersCallback ..." << endl;
		if (!standard()) throw exception();
		const int NInt = varInt.size, NCont = varCont.size;

		// Declare modeling components.
		IloNumVar eta(env, -IloInfinity, IloInfinity);
		IloNumVarArray X(env), Y(env, NCont, 0, IloInfinity);
		IloModel modelRMP(env), modelSP(env);
		vector<double> currentValInt = vector<double>(NInt, 0);
		double currentValEta = InfinityNeg;
		IloRangeArray consSP(env);

		// Initiate models.
		initiateModels(env, modelRMP, modelSP, X, Y, eta, currentValInt, currentValEta, consSP);
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

			// Add Benders cuts.

		}

		if (proceed) {
			// Restore integrality constraints.

			// Use callback to add lazy constraints.

		}

		incumbent.elapsedTime = runTime(start);
		cout << "elapsed time (Instance::solveBendersCallback): " << runTime(start) << endl;
		incumbent.print();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersCallback", exc);
	}

	env.end();
	return incumbent;
}


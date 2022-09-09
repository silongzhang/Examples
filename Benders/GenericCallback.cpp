#include"DataStructure.h"

// Implement the benders decomposition with generic callback.


tuple<SolutionStatus, double, vector<double>> Instance::solveNewSP(const vector<double>& valInt) const {
	tuple<SolutionStatus, double, vector<double>> result;
	try {

	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveSP", exc);
	}

	return result;
}


void BendersGenericCallback::invoke(const IloCplex::Callback::Context& context) {
	try {
		IloEnv env = context.getEnv();
		IloNumArray valueX(env);
		vector<double> valInt;
		double valEta;

		// Get the current values of variables.
		if (context.inCandidate()) {
			if (context.isCandidatePoint()) {
				context.getCandidatePoint(X, valueX);
				valInt = assign(valueX);
				valEta = context.getCandidatePoint(eta);
			}
			else throw exception();						// The LP relaxation of the RMP is unbounded.
		}
		else throw exception();							// Unexpected context ID.







	}
	catch (const exception& exc) {
		printErrorAndExit("BendersGenericCallback::invoke", exc);
	}
}


Solution Instance::solveBendersGenericCallback(const ParameterAlgorithm& parameter) const {
	Solution incumbent;
	incumbent.clear();
	IloEnv env;
	try {
		const clock_t start = clock();
		cout << "Running Instance::solveBendersGenericCallback ..." << endl;
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

		if (proceed) {
			// Restore integrality constraints.
			modelRMP.add(IloConversion(env, X, IloNumVar::Type::Int));








		}








		incumbent.elapsedTime = runTime(start);
		cout << "elapsed time (Instance::solveBendersGenericCallback): " << runTime(start) << endl;
		incumbent.print();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersGenericCallback", exc);
	}

	env.end();
	return incumbent;
}


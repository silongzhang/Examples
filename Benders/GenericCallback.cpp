#include"DataStructure.h"

// Implement the benders decomposition with generic callback.


tuple<SolutionStatus, double, vector<double>> Instance::solveNewSP(const vector<double>& valInt, Solution& incumbent) const {
	tuple<SolutionStatus, double, vector<double>> result;
	IloEnv env;
	try {
		// Declare modeling components.
		IloNumVarArray Y(env, varCont.size, 0, IloInfinity);
		IloModel modelSP(env);
		IloExpr expr = product(env, obj.coefCont, Y);
		modelSP.add(IloMinimize(env, expr));
		IloRangeArray consSP(env);
		for (const auto& cons : cpCons) {
			double rhs = cons.rhs - inner_product(cons.coefInt.begin(), cons.coefInt.end(), valInt.begin(), 0.0);
			consSP.add(genCons(env, cons.coefCont, Y, cons.type, rhs));
		}
		modelSP.add(consSP);

		// Solve the SP.
		IloCplex cplexSP(modelSP);
		cplexSP.setParam(IloCplex::RootAlg, IloCplex::Primal);
		IloNumArray dualSP(env);
		if (!solveSP(cplexSP, consSP, dualSP, incumbent, cpCons, valInt)) throw exception();

		get<2>(result) = assign(dualSP);
		if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible) {
			get<0>(result) = SolutionStatus::Infeasible;
		}
		else if (cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
			get<0>(result) = SolutionStatus::Optimal;
			get<1>(result) = cplexSP.getObjValue();
		}
		else throw exception();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveNewSP", exc);
	}

	env.end();
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

		// Solve the SP.
		const auto res = instance.solveNewSP(valInt, incumbent);
		IloNumArray dualSP(env, get<2>(res).size());
		for (int i = 0; i < dualSP.getSize(); ++i) dualSP[i] = get<2>(res)[i];

		if (get<0>(res) == SolutionStatus::Infeasible) {
			IloRange cons(env, -IloInfinity, instance.exprRhs(env, dualSP, X), 0);
			context.rejectCandidate(cons);
			++incumbent.nFeasCutIP;
		}
		else if (get<0>(res) == SolutionStatus::Optimal) {

		}
		else throw exception();







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


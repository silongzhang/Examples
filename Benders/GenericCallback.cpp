#include"DataStructure.h"

// Implement the benders decomposition with generic callback.


void BendersGenericCallback::invoke(const IloCplex::Callback::Context& context) {
	try {
		if (context.inCandidate() && context.isCandidatePoint()) {
			// Avoid race conditions in parallel search.
			static mutex mtx;
			lock_guard<mutex> lock(mtx);

			// Get current values of variables.
			IloEnv env = context.getEnv();
			IloNumArray valueX(env);
			context.getCandidatePoint(X, valueX);
			vector<double> valInt = assign(valueX);
			double valEta = context.getCandidatePoint(eta);

			// Define the subproblem.
			IloEnv envSP;
			IloNumVarArray Y(envSP, instance.varCont.size, 0, IloInfinity);
			IloModel modelSP(envSP);
			modelSP.add(IloMinimize(envSP, product(envSP, instance.obj.coefCont, Y)));
			IloRangeArray consSP(envSP);
			for (const auto& cons : instance.cpCons) {
				double rhs = cons.rhs - inner_product(cons.coefInt.begin(), cons.coefInt.end(), valInt.begin(), 0.0);
				consSP.add(genCons(envSP, cons.coefCont, Y, cons.type, rhs));
			}
			modelSP.add(consSP);

			// Solve the subproblem.
			IloCplex cplexSP(modelSP);
			cplexSP.setParam(IloCplex::RootAlg, IloCplex::Primal);
			solveModel(cplexSP);
			if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible) {
				IloNumArray dualSP(getDuals(cplexSP, consSP));
				IloRange cons(env, -IloInfinity, instance.exprRhs(env, dualSP, X), 0);
				context.rejectCandidate(cons);												// Add feasibility cut.
				++incumbent.nFeasCutIP;
			}
			else if (cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
				if (greaterThanReal(incumbent.objective, context.getCandidateObjective() - valEta + cplexSP.getObjValue(), PPM)) {
					incumbent.objective = context.getCandidateObjective() - valEta + cplexSP.getObjValue();			// Renew UB.
					incumbent.valueInt = valInt;
					incumbent.valueEta = valEta;
					incumbent.valueCont = ::getValues(cplexSP, Y);
					incumbent.status = SolutionStatus::Feasible;
					cout << "Found a better solution with objective = " << incumbent.objective << endl;
				}

				if (lessThanReal(valEta, cplexSP.getObjValue(), PPM)) {
					IloNumArray dualSP(getDuals(cplexSP, consSP));
					IloRange cons(env, -IloInfinity, instance.exprRhs(env, dualSP, X) - eta, 0);
					context.rejectCandidate(cons);											// Add optimality cut.
					++incumbent.nOptCutIP;
				}
				else if (greaterThanReal(valEta, cplexSP.getObjValue(), PPM)) throw exception();
			}
			else throw exception();
			envSP.end();
		}
		else throw exception();							// Unexpected context ID or unbounded LP relaxation.
	}
	catch (const exception& exc) {
		printErrorAndExit("BendersGenericCallback::invoke", exc);
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

			// Use callback to add lazy constraints.
			BendersGenericCallback bendersGCB(*this, incumbent, X, eta);
			CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate;
			cplexRMP.use(&bendersGCB, contextmask);

			// Stopping criteria.
			cplexRMP.setParam(IloCplex::Param::TimeLimit, parameter.timeLimit);
			cplexRMP.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, parameter.optThreshold);

			solveModel(cplexRMP);
			incumbent.LB = cplexRMP.getBestObjValue();
		}

		if (equalToReal(incumbent.objective, incumbent.LB, PPM)) incumbent.status = SolutionStatus::Optimal;
		incumbent.elapsedTime = runTime(start);
		cout << "elapsed time (Instance::solveBendersGenericCallback): " << runTime(start) << endl;
		cout << "# of nodes processed = " << cplexRMP.getNnodes() << '\t' << "# of nodes remained = " << cplexRMP.getNnodesLeft() << endl;
		incumbent.print();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersGenericCallback", exc);
	}

	env.end();
	return incumbent;
}


#include"DataStructure.h"

// Implement the benders decomposition with recursive procedure.


// Assumption 1. The MILP model has a finite optimal objective value.
// Assumption 2. # of integer variables > 0 and # of continuous variables > 0.
// Assumption 3. Variables are already transformed to be nonnegative.
// Assumption 4. Invoke the method "standardize" in advance if the instance is not "standard".
// Assumption 5. It gives the opposite number of the optimal objective value of the original instance if it is a maximization problem.
// Assumption 6. Relaxed master problems (RMPs) cannot be unbounded.
// Note 1. For a minimization problem, a sufficient condition for Assumption 6 is that integer variables are all nonnegative and their coefficients in the objective function are all nonnegative. A similar sufficient condition can be derived for a maximization problem.
// Note 2. Another sufficient condition for Assumption 6 is that the number of feasible combinations of values of integer variables is finite. For example, each integer variable is bounded from both below and above.
Solution Instance::solveBendersRecursive(const ParameterAlgorithm& parameter) const {
	Solution incumbent;
	incumbent.clear();
	IloEnv env;
	try {
		clock_t last = clock();
		cout << "Running Instance::solveBendersRecursive ..." << endl;
		if (!standard()) throw exception();
		const int NInt = varInt.size, NCont = varCont.size;

		// Define variables.
		IloNumVar eta(env, -IloInfinity, IloInfinity);
		IloNumVarArray X(env), Y(env, NCont, 0, IloInfinity);
		for (int i = 0; i < NInt; ++i)
			X.add(IloNumVar(env, 0, equalToReal(varInt.upperbounds[i], INFUB, PPM) ? IloInfinity : varInt.upperbounds[i]));

		// Set the objective functions.
		IloModel modelRMP(env), modelSP(env);
		IloExpr expr = product(env, obj.coefInt, X);
		expr += eta;
		modelRMP.add(IloMinimize(env, expr));

		expr = product(env, obj.coefCont, Y);
		modelSP.add(IloMinimize(env, expr));

		// Set constraints other than Benders cuts in the RMP.
		for (const auto& cons : itCons) addConstraint(modelRMP, cons.coefInt, X, cons.type, cons.rhs);
		modelRMP.add(eta >= InfinityNeg);						// Avoid the RMP to be unbounded.

		// Set the initial values of X and eta.
		vector<double> currentValInt = vector<double>(NInt, 0);
		double currentValEta = InfinityNeg;

		// Set constraints in the SP.
		IloRangeArray consSP(env);
		for (const auto& cons : cpCons) {
			double rhs = cons.rhs - inner_product(cons.coefInt.begin(), cons.coefInt.end(), currentValInt.begin(), 0.0);
			consSP.add(genCons(env, cons.coefCont, Y, cons.type, rhs));
		}
		modelSP.add(consSP);

		// Solve the RMP and SP alternately.
		IloCplex cplexRMP(modelRMP), cplexSP(modelSP);
		cplexSP.setParam(IloCplex::RootAlg, IloCplex::Primal);
		bool integral = false;									// Whether integrality constraints are restored.
		for (int iter = 1; !parameter.stop(last, incumbent.objective, incumbent.LB); ++iter) {
			cout << "+++++++++++++++++++++++++++++" << endl;
			cout << "Iter = " << iter << '\t' << "Integral = " << integral << '\t' << "UB = " << incumbent.objective << '\t' << "LB = " << incumbent.LB << '\t' << "nOptCutLP = " << incumbent.nOptCutLP << '\t' << "nOptCutIP = " << incumbent.nOptCutIP << '\t' << "nFeasCutLP = " << incumbent.nFeasCutLP << '\t' << "nFeasCutIP = " << incumbent.nFeasCutIP << endl;

			solveModel(cplexRMP);								// Solve the RMP.
			if (cplexRMP.getStatus() == IloAlgorithm::Status::Infeasible) {
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				cout << "The instance is infeasible! Terminate!" << endl;
				incumbent.status = SolutionStatus::Infeasible;
				break;
			}
			else if (cplexRMP.getStatus() == IloAlgorithm::Status::Unbounded) {
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				cout << "The current RMP might be unbounded (Assumption 6 is not satisfied)! Terminate!" << endl;
				incumbent.status = SolutionStatus::Unkown;
				break;
			}
			else if (cplexRMP.getStatus() == IloAlgorithm::Status::Optimal) {
				currentValInt = getValues(cplexRMP, X);
				currentValEta = cplexRMP.getValue(eta);
				incumbent.LB = cplexRMP.getObjValue();			// Renew LB.

				for (int i = 0; i < consSP.getSize(); ++i) {
					double rhs = cpCons[i].rhs - inner_product(cpCons[i].coefInt.begin(), cpCons[i].coefInt.end(), currentValInt.begin(), 0.0);
					setRhs(consSP[i], cpCons[i].type, rhs);		// Set the right hand side of constraints in the SP.
				}
			}
			else throw exception();

			solveModel(cplexSP);								// Solve the subproblem.
			if (cplexSP.getStatus() == IloAlgorithm::Status::Unbounded) {
				cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
				cout << "The instance might be unbounded! Terminate!" << endl;
				incumbent.status = SolutionStatus::Unkown;
				break;
			}
			else if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible) {
				IloNumArray dualSP(getDuals(cplexSP, consSP));							// Get dual values.
				modelRMP.add(0 >= exprRhs(env, dualSP, X));								// Add feasibility cut.
				integral ? ++incumbent.nFeasCutIP : ++incumbent.nFeasCutLP;
			}
			else if (cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
				if (integral && greaterThanReal(incumbent.objective, cplexRMP.getObjValue() - cplexRMP.getValue(eta) + cplexSP.getObjValue(), PPM))
					incumbent.renew(cplexRMP, X, eta, cplexSP, Y);						// Renew UB.

				IloNumArray dualSP(getDuals(cplexSP, consSP));							// Get dual values.
				if (lessThanReal(currentValEta, cplexSP.getObjValue(), PPM)) {
					modelRMP.add(eta >= exprRhs(env, dualSP, X));						// Add optimality cut.
					integral ? ++incumbent.nOptCutIP : ++incumbent.nOptCutLP;
				}
				else if (equalToReal(currentValEta, cplexSP.getObjValue(), PPM)) {
					if (integral) {
						incumbent.status = SolutionStatus::Optimal;
						break;															// Solved to optimality.
					}
					else {
						integral = true;												// Restore integrality constraints.
						modelRMP.add(IloConversion(env, X, IloNumVar::Type::Int));
					}
				}
				else throw exception();
			}
			else throw exception();
		}

		incumbent.elapsedTime = runTime(last);
		cout << "elapsed time (Instance::solveBendersRecursive): " << runTime(last) << endl;
		incumbent.print();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersRecursive", exc);
	}

	env.end();
	return incumbent;
}


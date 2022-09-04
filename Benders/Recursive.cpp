#include"DataStructure.h"

// Implement the benders decomposition with recursive procedure.


// 1. # of integer variables > 0 and # of continuous variables > 0.
// 2. Variables are already transformed to be nonnegative.
// 3. Invoke the method "standardize" in advance if the instance is not "standard".
// 4. It gives the opposite number of the optimal objective value of the original instance if it is a maximization problem.
Solution Instance::solveBendersRecursive() const {
	Solution result;
	result.clear();
	IloEnv env;
	try {
		clock_t last = clock();
		cout << "Running Instance::solveBendersRecursive ..." << endl;
		if (!standard()) throw exception();
		const int NInt = varInt.size, NCont = varCont.size;
		result.valueInt.resize(NInt), result.valueCont.resize(NCont);

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

		// Set the initial values of X and eta.
		result.valueInt = vector<double>(NInt, 0);
		result.valueEta = -IloInfinity;

		// Set constraints in the SP.
		IloRangeArray consSP(env);
		for (const auto& cons : cpCons) {
			double rhs = cons.rhs - inner_product(cons.coefInt.begin(), cons.coefInt.end(), result.valueInt.begin(), 0.0);
			consSP.add(genCons(env, cons.coefCont, Y, cons.type, rhs));
		}
		modelSP.add(consSP);

		// Solve the SP and RMP alternately.
		IloCplex cplexRMP(modelRMP), cplexSP(modelSP);
		cplexSP.setParam(IloCplex::RootAlg, IloCplex::Primal);
		for (int iter = 1; true; ++iter) {
			cout << "Iter = " << iter << endl;
			cplexSP.solve();
		}

		cout << "elapsed time (Instance::solveBendersRecursive): " << runTime(last) << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersRecursive", exc);
	}

	env.end();
	return result;
}


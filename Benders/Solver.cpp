#include"DataStructure.h"

// Solve the model with the Cplex solver directly.


IloRange genConsSolver(IloEnv env, IloIntVarArray X, IloNumVarArray Y, const Constraint& cons) {
	IloRange result;
	try {
		if (X.getSize() != cons.coefInt.size() || Y.getSize() != cons.coefCont.size()) throw exception();

		IloExpr expr(env);
		for (int i = 0; i < X.getSize(); ++i) expr += cons.coefInt[i] * X[i];
		for (int i = 0; i < Y.getSize(); ++i) expr += cons.coefCont[i] * Y[i];

		if (cons.type == ConstraintType::Eq) result = IloRange(env, cons.rhs, expr, cons.rhs);
		else if (cons.type == ConstraintType::Le) result = IloRange(env, -IloInfinity, expr, cons.rhs);
		else result = IloRange(env, cons.rhs, expr, IloInfinity);
	}
	catch (const exception& exc) {
		printErrorAndExit("generateConstraint", exc);
	}
	return result;
}


bool Instance::solveSolver() const {
	bool result;
	IloEnv env;
	try {
		clock_t last = clock();
		cout << "Running Instance::solveSolver ..." << endl;
		if (!consistent()) throw exception();
		const int NInt = varInt.size, NCont = varCont.size;

		// Define variables.
		IloIntVarArray X(env);
		for (int i = 0; i < NInt; ++i)
			X.add(IloIntVar(env, 0, equalToReal(varInt.upperbounds[i], INFUB, PPM) ? IloIntMax : varInt.upperbounds[i]));

		IloNumVarArray Y(env);
		for (int i = 0; i < NCont; ++i)
			Y.add(IloNumVar(env, 0, equalToReal(varCont.upperbounds[i], INFUB, PPM) ? IloInfinity : varCont.upperbounds[i]));

		// Set the objective function.
		IloModel model(env);
		IloExpr expr(env);
		for (int i = 0; i < NInt; ++i) expr += obj.coefInt[i] * X[i];
		for (int i = 0; i < NCont; ++i) expr += obj.coefCont[i] * Y[i];
		obj.type == ObjectiveType::Minimization ? model.add(IloMinimize(env, expr)) : model.add(IloMaximize(env, expr));

		// Set constraints.
		IloRangeArray constraints(env);
		for (int i = 0; i < cpCons.size(); ++i) constraints.add(genConsSolver(env, X, Y, cpCons[i]));
		for (int i = 0; i < itCons.size(); ++i) constraints.add(genConsSolver(env, X, Y, itCons[i]));
		model.add(constraints);

		// Solve the model.
		IloCplex cplex(model);
		result = solveModel(cplex);

		cout << "elapsed time (Instance::solveSolver): " << runTime(last) << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveSolver", exc);
	}

	env.end();
	return result;
}


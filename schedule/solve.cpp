#include"input.h"

void Input::solve() {
	IloEnv env;
	try {
		time_point<system_clock> startTime = system_clock::now();
		cout << "Running Input::solve ...." << endl;

		// Define the model.
		IloModel model(env);

		// Define the variables.
		vector<vector<IloBoolVarArray>> X(NPerson);
		vector<IloBoolVarArray> Y(NPerson);
		for (int i = 0; i < NPerson; ++i) {
			X[i].resize(NDay);
			for (int j = 0; j < NDay; ++j) {
				X[i][j] = IloBoolVarArray(env, NShift);
			}

			Y[i] = IloBoolVarArray(env, NDay - 1);
		}

		// Define the objective function.
		{
			IloExpr expr(env);
			for (int i = 0; i < NPerson; ++i) {
				for (int j = 0; j < NDay - 1; ++j) {
					model.add(Y[i][j] >= X[i][j][0] - X[i][j + 1][0]);
					expr += Y[i][j];
				}
			}
			model.add(IloMinimize(env, expr));
			expr.end();
		}

		// Define the constraints.
		
		// Each person works at exactly one shift per day.
		for (int i = 0; i < NPerson; ++i) {
			for (int j = 0; j < NDay; ++j) {
				model.add(IloSum(X[i][j]) == 1);
			}
		}

		// A person cannot work MCWD consecutive days.
		for (int i = 0; i < NPerson; ++i) {
			for (int j = 0; j <= NDay - MCWD; ++j) {
				IloExpr expr(env);

				for (int k = 0; k < MCWD; ++k)
					expr += X[i][j + k][0];

				model.add(expr >= 1);
				expr.end();
			}
		}

		// Each shift must be assigned to the required number of persons per day.
		for (int j = 0; j < NDay; ++j) {
			for (int s = 0; s < NShift; ++s) {
				IloExpr expr(env);

				for (int i = 0; i < NPerson; ++i)
					expr += X[i][j][s];

				model.add(expr == ShiftPerDay[s]);
				expr.end();
			}
		}

		// Each person must be assigned to the required shifts.
		for (int i = 0; i < NPerson; ++i) {
			for (int s = 0; s < NShift; ++s) {
				IloExpr expr(env);

				for (int j = 0; j < NDay; ++j)
					expr += X[i][j][s];

				model.add(expr == PersonShift[i][s]);
				expr.end();
			}
		}

		// Adjacent shifts cannot be assigned to the same person.
		for (int i = 0; i < NPerson; ++i)
			for (int j = 0; j < NDay - 1; ++j)
				for (int s1 = 0; s1 < NShift; ++s1)
					for (int s2 = 0; s2 < NShift; ++s2)
						if (!Adjacent[s1][s2])
							model.add(X[i][j][s1] + X[i][j + 1][s2] <= 1);

		// Presettings.
		for (const auto& elem : PersonDayShift)
			model.add(X[get<0>(elem)][get<1>(elem)][get<2>(elem)] == 1);

		// Solve the model.
		IloCplex cplex(model);
		if (!cplex.solve())
			throw exception();
		cout << "objective = " << cplex.getObjValue() << '\t' << "solution status = " << cplex.getStatus() << endl;

		// Extract the solution.
		solution.resize(NPerson);
		for (int i = 0; i < NPerson; ++i) {
			solution[i].resize(NDay);
			for (int j = 0; j < NDay; ++j) {
				for (int s = 0; s < NShift; ++s) {
					if (equalToReal(cplex.getValue(X[i][j][s]), 1, PPM))
						solution[i][j] = s;
					else if (!equalToReal(cplex.getValue(X[i][j][s]), 0, PPM))
						throw exception();
				}
			}
		}

		cout << "elapsed time (Input::solve): " << runTime(startTime) << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("Input::solve", exc);
	}
	env.end();
}

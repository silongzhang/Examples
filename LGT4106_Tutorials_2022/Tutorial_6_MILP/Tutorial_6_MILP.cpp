#include<ilcplex/ilocplex.h>
#include<iostream>
#include<vector>

typedef IloArray<IloNumVarArray> IloNumVarArray2;

using std::cout;
using std::endl;
using std::vector;


int main(int argc, char** argv) {
	const int NPlant = 2, NWarehouse = 2, NMarket = 3;						// The numbers of plants, warehouses, markets.
	const vector<double> capacity = { 200000,200000 };						// Capacities of plants.
	const vector<double> demand = { 50000,100000,50000 };					// Demand quantities of markets.
	const vector<vector<double>> costTransPltWrh = { {0,5},{4,2} };			// Transportation rates between plants and warehouses.
	const vector<vector<double>> costTransWrhMrk = { {3,4,5},{2,1,2} };		// Transportation rates between warehouses and markets.
	//const vector<double> fixedCost = { 0,0 };								// Fixed costs of plants.
	//const vector<double> fixedCost = { 100000,100000 };						// Fixed costs of plants.
	const vector<double> fixedCost = { 100000,300000 };						// Fixed costs of plants.

	IloEnv env;
	IloModel model(env);

	// Define variables.
	IloNumVarArray2 X(env), Y(env);
	X.setSize(NPlant), Y.setSize(NWarehouse);
	for (int i = 0; i < NPlant; ++i) X[i] = IloNumVarArray(env, NWarehouse, 0, IloInfinity);
	for (int j = 0; j < NWarehouse; ++j) Y[j] = IloNumVarArray(env, NMarket, 0, IloInfinity);
	IloNumVarArray Z(env, NPlant, 0, 1, ILOBOOL);

	// Define the objective.
	IloExpr expr(env);
	for (int i = 0; i < NPlant; ++i)
		for (int j = 0; j < NWarehouse; ++j)
			expr += costTransPltWrh[i][j] * X[i][j];

	for (int j = 0; j < NWarehouse; ++j)
		for (int k = 0; k < NMarket; ++k)
			expr += costTransWrhMrk[j][k] * Y[j][k];

	for (int i = 0; i < NPlant; ++i)
		expr += fixedCost[i] * Z[i];

	model.add(IloMinimize(env, expr));

	// Capacity constraints associated with plants.
	for (int i = 0; i < NPlant; ++i) {
		expr.clear();
		for (int j = 0; j < NWarehouse; ++j)
			expr += X[i][j];

		model.add(expr <= capacity[i] * Z[i]);
	}

	// Flow balance constraints at warehouses.
	for (int j = 0; j < NWarehouse; ++j) {
		expr.clear();
		for (int i = 0; i < NPlant; ++i)
			expr += X[i][j];

		for (int k = 0; k < NMarket; ++k)
			expr -= Y[j][k];

		model.add(expr == 0);
	}

	// Demand constraints associated with markets.
	for (int k = 0; k < NMarket; ++k) {
		expr.clear();
		for (int j = 0; j < NWarehouse; ++j)
			expr += Y[j][k];

		model.add(expr == demand[k]);
	}

	// Solve the model.
	IloCplex cplex(model);
	cplex.solve();

	// Print the solution.
	for (int i = 0; i < NPlant; ++i) {
		cout << "Z[" << i << "] = " << cplex.getValue(Z[i]) << '\t';
		for (int j = 0; j < NWarehouse; ++j)
			cout << "X[" << i << "][" << j << "] = " << cplex.getValue(X[i][j]) << '\t';

		cout << endl;
	}
	for (int j = 0; j < NWarehouse; ++j) {
		for (int k = 0; k < NMarket; ++k)
			cout << "Y[" << j << "][" << k << "] = " << cplex.getValue(Y[j][k]) << '\t';

		cout << endl;
	}
	cout << "The objective function value = " << cplex.getObjValue() << endl;

	env.end();
	return 0;
}


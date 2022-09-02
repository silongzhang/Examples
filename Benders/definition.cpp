#include"definition.h"


void printErrorAndExit(const string &str, const exception &exc) {
	cout << "There is an error in " + str + " ! " << endl;
	cout << "Error information: " << exc.what() << endl;
	exit(1);
}


void terminate(const string& str) {
	cout << "Termination reason: " << str << endl;
	exit(0);
}


double runTime(const clock_t &start)
{
	clock_t finish = clock();
	return (double)(finish - start) / CLOCKS_PER_SEC;
}


string getNowTime()
{
	ostringstream oss;
	oss << time(0);
	return oss.str();
}


bool lessThanReal(const double &lhs, const double &rhs, const double &threshold) {
	return lhs < rhs - threshold;
}


bool greaterThanReal(const double &lhs, const double &rhs, const double &threshold) {
	return lhs > rhs + threshold;
}


bool equalToReal(const double &lhs, const double &rhs, const double &threshold) {
	return abs(lhs - rhs) <= threshold;
}


bool operator<=(const bitset<NMAX>& lhs, const bitset<NMAX>& rhs) {
	return (lhs | rhs) == rhs;
}


bool solveModel(IloCplex& cplex) {
	bool result = false;
	try {
		result = cplex.solve();
		cout << "solution status = " << cplex.getStatus() << endl;
		if (result) cout << "objective = " << cplex.getObjValue() << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("solveModel", exc);
	}
	return result;
}


// Euclidean distance.
double EuclideanDistance(const double x1, const double y1, const double x2, const double y2) {
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}


// Set precision.
double setPrecision(const double num, const bool upDown, const int precision) {
	return upDown ? ceil(num * pow(10, precision)) / pow(10, precision) : floor(num * pow(10, precision)) / pow(10, precision);
}


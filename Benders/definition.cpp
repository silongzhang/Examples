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


bool solveModel(IloCplex cplex) {
	bool feasible = false;
	try {
		feasible = cplex.solve();
		cout << "solution status = " << cplex.getStatus() << endl;
		if (feasible) cout << "objective = " << cplex.getObjValue() << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("solveModel", exc);
	}
	return feasible;
}


// Euclidean distance.
double EuclideanDistance(const double x1, const double y1, const double x2, const double y2) {
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}


// Set precision.
double setPrecision(const double num, const bool upDown, const int precision) {
	return upDown ? ceil(num * pow(10, precision)) / pow(10, precision) : floor(num * pow(10, precision)) / pow(10, precision);
}


bool isInteger(double num, double threshold) {
	return equalToReal(num, ceil(num), threshold) || equalToReal(num, floor(num), threshold);
}


bool isInteger(const vector<double>& nums, double threshold) {
	for (const auto num : nums)
		if (!isInteger(num, threshold))
			return false;

	return true;
}


void print(SolutionStatus status) {
	switch (status)
	{
	case SolutionStatus::Unkown:
		cout << "Unkown";
		break;
	case SolutionStatus::Infeasible:
		cout << "Infeasible";
		break;
	case SolutionStatus::Unbounded:
		cout << "Unbounded";
		break;
	case SolutionStatus::Feasible:
		cout << "Feasible";
		break;
	case SolutionStatus::Optimal:
		cout << "Optimal";
		break;
	default:
		break;
	}
}


vector<double> assign(const IloNumArray nums) {
	vector<double> result(nums.getSize());
	for (int i = 0; i < nums.getSize(); ++i)
		result[i] = nums[i];

	return result;
}


#pragma once

#include"Header.h"

constexpr auto NMAX = 128;

constexpr auto DECI = 1e-1;
constexpr auto CENTI = 1e-2;
constexpr auto MILLI = 1e-3;
constexpr auto PPM = 1e-6;
constexpr auto PPB = 1e-9;
constexpr auto TenTh = 1e4;
constexpr auto InfinityPos = INFINITY;
constexpr auto InfinityNeg = -InfinityPos;

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray3> IloNumVarArray4;
typedef IloArray<IloIntVarArray> IloIntVarArray2;
typedef IloArray<IloBoolVarArray> IloBoolVarArray2;
typedef IloArray<IloBoolVarArray2> IloBoolVarArray3;

enum class SolutionStatus { Unkown, Infeasible, Feasible, Optimal };

void printErrorAndExit(const string &str, const exception &exc);
void terminate(const string& str);
double runTime(const clock_t &start);
string getNowTime();
bool lessThanReal(const double &lhs, const double &rhs, const double &threshold);
bool greaterThanReal(const double &lhs, const double &rhs, const double &threshold);
bool equalToReal(const double &lhs, const double &rhs, const double &threshold);
bool operator<=(const bitset<NMAX>& lhs, const bitset<NMAX>& rhs);
bool solveModel(IloCplex& cplex);
double EuclideanDistance(const double x1, const double y1, const double x2, const double y2);
double setPrecision(const double num, const bool upDown, const int precision);

template<typename T>
void print1(ostream &os, const T &cont, const char character) {
	for (const auto &elem : cont) {
		os << elem << character;
	}
}

template<typename T>
void print2(ostream &os, const T &cont, const char character) {
	for (const auto &first : cont) {
		for (const auto &second : first) {
			os << second << character;
		}
		os << endl;
	}
}

template<typename T>
void read2(istream &ins, T &cont) {
	for (int i = 0; i < cont.size(); ++i) {
		for (int j = 0; j < cont[i].size(); ++j) {
			ins >> cont[i][j];
		}
	}
}

// transfer string to number
template<typename T>
void strToNum(const string &str, T &num) {
	stringstream ss;
	ss << str;
	ss >> num;
}

// transfer number to string
template<typename T>
string numToStr(const T &num) {
	string str;
	stringstream ss;
	ss << num;
	ss >> str;
	return str;
}

class SortReal {
public:
	bool operator()(const double lhs, const double rhs) const {
		return lessThanReal(lhs, rhs, PPM);
	}
};

class SortRealPair {
public:
	bool operator()(const pair<double, double>& lhs, const pair<double, double>& rhs) const {
		if (lessThanReal(lhs.first, rhs.first, PPM)) return true;
		else if (equalToReal(lhs.first, rhs.first, PPM) && lessThanReal(lhs.second, rhs.second, PPM)) return true;
		else return false;
	}
};

class thread_guard
{
private:
	thread &t;
public:
	// constructor
	explicit thread_guard(thread &myT) :t(myT) {}
	// destructor
	~thread_guard()
	{
		if (t.joinable())
		{
			t.join();
		}
	}
	thread_guard(const thread_guard &) = delete;
	thread_guard& operator=(const thread_guard &) = delete;
};


#include"utils.h"

void printErrorAndExit(const string &str, const exception &exc) {
	cerr << "Function: " << str << '\t' << "Exception: " << exc.what() << endl;
	throw;
}

void terminate(const string& str) {
	cerr << "Termination reason: " << str << endl;
	exit(1);
}

double runTime(const time_point<system_clock>& start) {
	std::chrono::duration<double> drt = system_clock::now() - start;
	return drt.count();
}

string getNowTime() {
	ostringstream oss;
	oss << time(0);
	return oss.str();
}

// Euclidean distance.
double EuclideanDistance(const double x1, const double y1, const double x2, const double y2) {
	return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

// Set precision.
double setPrecision(const double num, const bool upDown, const int precision) {
	return upDown ? ceil(num * pow(10, precision)) / pow(10, precision) : floor(num * pow(10, precision)) / pow(10, precision);
}

// Whether it is an elementary sequence without duplicated elements.
bool elementary(const vector<int>& vec) {
	unordered_set<int> ust(vec.begin(), vec.end());
	return ust.size() == vec.size();
}

unordered_set<int> intersection(const vector<int>& lhs, const vector<int>& rhs) {
	unordered_set<int> result;
	unordered_set<int> left(lhs.begin(), lhs.end());
	for (auto elem : rhs)
		if (left.find(elem) != left.end())
			result.insert(elem);

	return result;
}

map<int, int> histogram(const vector<int>& nums) {
	map<int, int> ret;

	for (auto num : nums) {
		auto pos = ret.find(num);
		if (pos == ret.end())
			ret.insert(make_pair(num, 1));
		else
			++(pos->second);
	}

	return ret;
}

map<int, int> histogram(const vector<vector<int>>& nums) {
	map<int, int> ret;

	for (const auto& elem : nums) {
		for (auto num : elem) {
			auto pos = ret.find(num);
			if (pos == ret.end())
				ret.insert(make_pair(num, 1));
			else
				++(pos->second);
		}
	}

	return ret;
}

void printHistogram(ofstream& os, const map<int, int>& mp) {
	try {
		for (const auto& elem : mp)
			os << elem.first << '\t' << elem.second << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("printHistogram", exc);
	}
}

set<int> getKeys(const map<int, set<int>>& mp) {
	set<int> ret;
	for (const auto& elem : mp)
		ret.insert(elem.first);

	return ret;
}

int numPanelDemand(int numPortDemand, int numUnUsedPort, int sizePanel) {
	return int(ceil(max(0, numPortDemand - numUnUsedPort) / double(sizePanel)));
}

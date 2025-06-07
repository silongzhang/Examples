#pragma once

#include"definition.h"

void printErrorAndExit(const string& str, const exception& exc);
void terminate(const string& str);
double runTime(const time_point<system_clock>& start);
string getNowTime();
double EuclideanDistance(const double x1, const double y1, const double x2, const double y2);
double setPrecision(const double num, const bool upDown, const int precision);
bool elementary(const vector<int>& vec);
unordered_set<int> intersection(const vector<int>& lhs, const vector<int>& rhs);
map<int, int> histogram(const vector<int>& nums);
map<int, int> histogram(const vector<vector<int>>& nums);
void printHistogram(ofstream& os, const map<int, int>& mp);
set<int> getKeys(const map<int, set<int>>& mp);
int numPanelDemand(int numPortDemand, int numUnUsedPort, int sizePanel);

inline bool lessThanReal(const double& lhs, const double& rhs, const double& threshold) {
	return lhs < rhs - threshold;
}

inline bool greaterThanReal(const double& lhs, const double& rhs, const double& threshold) {
	return lhs > rhs + threshold;
}

inline bool equalToReal(const double& lhs, const double& rhs, const double& threshold) {
	return abs(lhs - rhs) <= threshold;
}

inline bool operator<=(const bitset<NMAX>& lhs, const bitset<NMAX>& rhs) {
	return (lhs | rhs) == rhs;
}

template<typename T>
void print1(ostream& os, const T& cont, const char character) {
	for (const auto& elem : cont) {
		os << elem << character;
	}
}

template<typename T>
void print2(ostream& os, const T& cont, const char character) {
	for (const auto& first : cont) {
		for (const auto& second : first) {
			os << second << character;
		}
		os << endl;
	}
}

template<typename T>
void read2(istream& ins, T& cont) {
	for (int i = 0; i < cont.size(); ++i) {
		for (int j = 0; j < cont[i].size(); ++j) {
			ins >> cont[i][j];
		}
	}
}

// transfer string to number
template<typename T>
void strToNum(const string& str, T& num) {
	stringstream ss;
	ss << str;
	ss >> num;
}

// transfer number to string
template<typename T>
string numToStr(const T& num) {
	string str;
	stringstream ss;
	ss << num;
	ss >> str;
	return str;
}

class threadGuard {
private:
	std::thread& t;
public:
	explicit threadGuard(std::thread& _t) :t(_t) {}
	~threadGuard() {
		if (t.joinable()) {
			t.join();
		}
	}
	threadGuard(const threadGuard&) = delete;
	threadGuard& operator=(const threadGuard&) = delete;
	threadGuard& operator=(threadGuard&& rhs) noexcept {
		if (this == &rhs)
			return *this;
		t = std::move(rhs.t);
		return *this;
	}
};

#pragma once

#include<ilcplex/ilocplex.h>
#include"utils.h"

class Input {
private:
	vector<vector<int>> solution;						// shift[person][day].

	bool check() const;
public:
	int NDay;
	int NPerson;
	int NShift;
	int MCWD;

	int NSol;

	vector<int> DayName;
	vector<int> PersonName;
	vector<char> ShiftName;
	map<char, int> NameShift;

	vector<int> ShiftPerDay;
	vector<vector<int>> PersonShift;
	vector<vector<bool>> Adjacent;
	vector<tuple<int, int, int>> PersonDayShift;

	Input() : NDay(0), NPerson(0), NShift(0), MCWD(0) {}
	Input(const string& inputFile);
	void solve();
	void printSolution(const string& outputFile) const;
};

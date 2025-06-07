#include"input.h"

Input::Input(const string& inputFile) {
	try {
		ifstream is(inputFile);
		if (!is) {
			cout << "Cannot open file: " << inputFile << endl;
			throw exception();
		}

		int itTemp;
		char chTemp;
		string line;
		getline(is, line);

		is >> NDay >> NPerson >> NShift >> MCWD;

		DayName.resize(NDay);
		for (int i = 0; i < NDay; ++i)
			DayName[i] = i + 1;

		PersonName.resize(NPerson);
		for (int i = 0; i < NPerson; ++i)
			PersonName[i] = i + 1;

		getline(is, line);
		getline(is, line);

		ShiftName.resize(NShift);
		for (int i = 0; i < NShift; ++i)
			is >> ShiftName[i];

		for (int i = 0; i < NShift; ++i)
			NameShift[ShiftName[i]] = i;

		if (NameShift.size() != NShift) {
			cout << "Invalid shift information." << endl;
			throw exception();
		}

		getline(is, line);

		ShiftPerDay.resize(NShift);
		for (int i = 0; i < NShift; ++i)
			is >> ShiftPerDay[i];

		getline(is, line);
		getline(is, line);
		getline(is, line);

		PersonShift.resize(NPerson, vector<int>(NShift, 0));
		for (int i = 0; i < NPerson; ++i) {
			is >> itTemp;
			for (int j = 0; j < NShift; ++j)
				is >> PersonShift[i][j];
		}

		getline(is, line);
		getline(is, line);
		getline(is, line);

		Adjacent.resize(NShift, vector<bool>(NShift, false));
		for (int i = 0; i < NShift; ++i) {
			is >> chTemp;
			for (int j = 0; j < NShift; ++j) {
				is >> itTemp;

				if (itTemp < 0 || itTemp > 1) {
					cout << "Invalid adjacency information." << endl;
					throw exception();
				}

				Adjacent[i][j] = (itTemp == 1) ? true : false;
			}
		}

		getline(is, line);
		getline(is, line);
		getline(is, line);

		while (true) {
			is >> itTemp;
			if (is.eof())
				break;

			int person = itTemp - 1;
			int day, shift;
			is >> day;
			day -= 1;
			is >> chTemp;
			shift = NameShift[chTemp];
			PersonDayShift.push_back(make_tuple(person, day, shift));
		}

		is.close();

		if (!check()) {
			cout << "Input data is invalid." << endl;
			throw exception();
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("Input::Input", exc);
	}
}

bool Input::check() const {
	int sum = 0;
	for (auto num : ShiftPerDay)
		sum += num;

	if (sum != NPerson) {
		cout << "Invalid shift information." << endl;
		return false;
	}

	for (int i = 0; i < NPerson; ++i) {
		int totalShifts = 0;

		for (int j = 0; j < NShift; ++j)
			totalShifts += PersonShift[i][j];

		if (totalShifts != NDay) {
			cout << "Invalid shift assignment for person " << PersonName[i] << "." << endl;
			return false;
		}
	}

	return true;
}

void Input::printSolution(const string& outputFile) const {
	try {
		ofstream os(outputFile);
		if (!os) {
			cout << "Cannot open output file: " << outputFile << endl;
			throw exception();
		}

		os << "Person" << ',';
		for (int i = 0; i < NDay; ++i)
			os << DayName[i] << ',';
		os << endl;

		for (int i = 0; i < NPerson; ++i) {
			os << PersonName[i] << ',';
			for (int j = 0; j < NDay; ++j) {
				int s = solution[i][j];
				os << ShiftName[s] << ',';
			}
			os << endl;
		}

		os.close();
	}
	catch (const exception& exc) {
		printErrorAndExit("Input::printSolution", exc);
	}
}

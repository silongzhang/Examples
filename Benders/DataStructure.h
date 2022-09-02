#pragma once

#include"definition.h"

enum class VariableType { Continuous, Integer };
enum class ObjectiveType { Minimization, Maximization };
enum class ContraintType { Eq, Le, Ge };

class Parameter {
public:

};


class Variable {
public:
	VariableType type;
	int size;
	vector<double> values;

	Variable() :type(VariableType::Continuous), size(0) {}
	Variable(VariableType tp, int sz) :type(tp), size(sz), values(vector<double>(sz, 0)) {}
};


class Objective {
public:
	ObjectiveType type;
	vector<double> coefInt;
	vector<double> coefCont;

	Objective() :type(ObjectiveType::Minimization) {}
	Objective(ObjectiveType tp, const vector<double>& cfIt, const vector<double>& cfCt) :type(tp), coefInt(cfIt), coefCont(cfCt) {}
};


class Contraint {
public:
	ContraintType type;
	vector<double> coefInt;
	vector<double> coefCont;
	double rhs;

	Contraint() :type(ContraintType::Eq), rhs(0) {}
	Contraint(ContraintType tp, const vector<double>& cfIt, const vector<double>& cfCt, double rgh) :type(tp), coefInt(cfIt), coefCont(cfCt), rhs(rgh) {}
};


class Problem {
public:
	Variable varInt;
	Variable varCont;
	Objective obj;
	vector<Contraint> constraints;

	Problem() {}
};



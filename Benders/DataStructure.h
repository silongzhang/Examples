#pragma once

#include"definition.h"

enum class VariableType { Continuous, Integer };
enum class ObjectiveType { Minimization, Maximization };
enum class ContraintType { Eq, Le, Ge };


// Assumption: variables are already transformed to be nonnegative.
// The nonnegativity of continuous variables ensures the existence of extreme points of feasible subprolems.
class Variable {
public:
	VariableType type;
	int size;
	vector<double> upperbounds;

	Variable() :type(VariableType::Continuous), size(0) {}
	Variable(VariableType tp, int sz, const vector<double>& ub) :type(tp), size(sz), upperbounds(ub) {}
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


class Instance {
public:
	Variable varInt;
	Variable varCont;
	Objective obj;
	vector<Contraint> cpCons;				// Coupling constraints.
	vector<Contraint> itCons;				// Constraints involving only integer variables.

	Instance() {}
	Instance(const Variable& vrIt, const Variable& vrCt, const Objective& objective, const vector<Contraint>& cpCs, const vector<Contraint>& itCs) :varInt(vrIt), varCont(vrCt), obj(objective), cpCons(cpCs), itCons(itCs) {
		if (!consistent()) printErrorAndExit("Instance", exception());
	}
	bool consistent() const;
};


// Test
class ParameterTest {
public:

};



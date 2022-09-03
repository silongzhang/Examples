#pragma once

#include"definition.h"

// Assumption: # of integer variables > 0 and # of continuous variables > 0.

constexpr auto INFUB = -1;

enum class VariableType { Continuous, Integer };
enum class ObjectiveType { Minimization, Maximization };
enum class ConstraintType { Eq, Le, Ge };


// Assumption: variables are already transformed to be nonnegative.
// The nonnegativity of continuous variables ensures the existence of extreme points of feasible subprolems.
class Variable {
public:
	VariableType type;
	int size;
	vector<double> upperbounds;					// 0 <= Variable[i] <= + \infty if upperbounds[i] = INFUB.

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


class Constraint {
public:
	ConstraintType type;
	vector<double> coefInt;
	vector<double> coefCont;
	double rhs;

	Constraint() :type(ConstraintType::Eq), rhs(0) {}
	Constraint(ConstraintType tp, const vector<double>& cfIt, const vector<double>& cfCt, double rgh) :type(tp), coefInt(cfIt), coefCont(cfCt), rhs(rgh) {}
};


class Instance {
public:
	Variable varInt;
	Variable varCont;
	Objective obj;
	vector<Constraint> cpCons;				// Coupling constraints.
	vector<Constraint> itCons;				// Constraints involving only integer variables.

	Instance() {}
	Instance(const Variable& vrIt, const Variable& vrCt, const Objective& objective, const vector<Constraint>& cpCs, const vector<Constraint>& itCs) :varInt(vrIt), varCont(vrCt), obj(objective), cpCons(cpCs), itCons(itCs) {
		if (!consistent()) printErrorAndExit("Instance", exception());
	}
	bool consistent() const;
	bool solveSolver() const;
};


// Test
class ParameterTest {
public:
	unsigned seed;
	int NInt, NCont;
	double MaxUBInt, MaxUBCont;
	double ratioInfUBInt, ratioInfUBCont;
	bool minimization;
	pair<double, double> rangeObjCoefInt, rangeObjCoefCont;
	pair<double, double> rangeConsCoefInt, rangeConsCoefCont;
	double sparsity;
	double redundancy;
	double ratioEq, ratioLe, ratioGe;
	int NCpCons, NItCons;

	ParameterTest() :seed(time(0)), NInt(0), NCont(0), MaxUBInt(0), MaxUBCont(0), ratioInfUBInt(0), ratioInfUBCont(0), minimization(true), sparsity(0), redundancy(0), ratioEq(0), ratioLe(0), ratioGe(0), NCpCons(0), NItCons(0) {}
	bool valid() const;
	Instance generateInstance() const;
};


vector<double> generateRandomVector(default_random_engine& engine, int N, pair<double, double> range, double sparsity, bool canEmpty);
void setRandomSignAndRhs(default_random_engine& engine, const vector<double>& midInt, const vector<double>& midCont, Constraint& cons, double ratioEq, double ratioLe, double ratioGe, double redundancy);


// Solver
IloRange genConsSolver(IloEnv env, IloIntVarArray X, IloNumVarArray Y, const Constraint& cons);


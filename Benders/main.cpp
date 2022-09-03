#include"DataStructure.h"

// Command line arguments are required.

int main(int argc, char** argv) {
	try {
		ParameterTest parameter;
		parameter.seed = 11;
		parameter.NInt = 1e3;
		parameter.NCont = 1e4;
		parameter.MaxUBInt = 10;
		parameter.MaxUBCont = 100;
		parameter.ratioInfUBInt = 0.1;
		parameter.ratioInfUBCont = 0.1;
		parameter.minimization = true;
		parameter.rangeObjCoefInt = parameter.rangeObjCoefCont = parameter.rangeConsCoefInt = parameter.rangeConsCoefCont = { 1,100 };
		parameter.sparsity = 0.2;
		parameter.redundancy = 0.05;
		parameter.ratioEq = 0.05;
		parameter.ratioLe = 0.25;
		parameter.ratioGe = 0.7;
		parameter.NCpCons = 1e3;
		parameter.NItCons = 1e2;

		Instance instance = parameter.generateInstance();
	}
	catch (const exception& exc) {
		printErrorAndExit("main", exc);
	}
	return 0;
}


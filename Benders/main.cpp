#include"DataStructure.h"

// Command line arguments are required.

int main(int argc, char** argv) {
	try {
		ParameterAlgorithm prmAlg;
		prmAlg.timeLimit = InfinityPos;
		prmAlg.optThreshold = 0;

		ParameterTest prmTest;
		prmTest.seed = 11;
		prmTest.NInt = 1e2;
		prmTest.NCont = 1e2;
		prmTest.MaxUBInt = 10;
		prmTest.MaxUBCont = 10;
		prmTest.ratioInfUBInt = 0.1;
		prmTest.ratioInfUBCont = 0.1;
		prmTest.minimization = true;
		prmTest.rangeObjCoefInt = prmTest.rangeObjCoefCont = prmTest.rangeConsCoefInt = prmTest.rangeConsCoefCont = { 1,100 };
		prmTest.sparsity = 0.2;
		prmTest.redundancy = 0.25;
		prmTest.ratioEq = 0.05;
		prmTest.ratioLe = 0.25;
		prmTest.ratioGe = 0.7;
		prmTest.NCpCons = 1e2;
		prmTest.NItCons = 1e2;

		Instance instance = prmTest.generateInstance();
		cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		instance.solveSolver();
		instance.standardize();
		cout << endl << endl << endl << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		instance.solveBendersRecursive(prmAlg);
		cout << endl << endl << endl << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		instance.solveBendersLegacyCallback(prmAlg);
		cout << endl << endl << endl << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		prmAlg.iterPrintBC = 1e3;
		instance.solveBendersBC(prmAlg);
	}
	catch (const exception& exc) {
		printErrorAndExit("main", exc);
	}
	return 0;
}


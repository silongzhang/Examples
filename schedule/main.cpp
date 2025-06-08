#include"input.h"

int main(int argc, char** argv) {
	try {
		// Run through command lines.
		if (argc != 2)
			throw exception();
		const string inst(argv[1]);

		//const string inst = "202505";

		const string inputFile = "input_" + inst + ".txt";

		const string outputFile = "output_" + inst + ".csv";

		Input input(inputFile);
		input.solve();
		input.printSolution(outputFile);
	}
	catch (const exception& exc) {
		cerr << "Function: " << "main" << '\t' << "Exception: " << exc.what() << endl;
		exit(1);
	}
	return 0;
}

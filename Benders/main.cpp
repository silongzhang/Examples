#include"DataStructure.h"

// Command line arguments are required.

int main(int argc, char** argv) {
	try {
		if (argc < 3) throw exception();

	}
	catch (const exception& exc) {
		printErrorAndExit("main", exc);
	}
	return 0;
}


#include"DataStructure.h"

// Implement the benders decomposition with branch-and-cut.


TreeNode::TreeNode(const Instance& instance) :depth(1), valueEta(InfinityNeg), objective(InfinityPos), integral(false) {
	branchCons.resize(instance.varInt.size);

	for (int i = 0; i < instance.varInt.size; ++i) {
		branchCons[i].first = 0;
		const double ub(instance.varInt.upperbounds[i]);
		branchCons[i].second = equalToReal(ub, ceil(ub), PPM) ? ceil(ub) : floor(ub);
	}
}


// A naive node selection strategy. Select the node with the smallest lower bound.
TreeNode Tree::selectNode() {
	try {
		if (nodes.empty()) throw exception();
	}
	catch (const exception& exc) {
		printErrorAndExit("Tree::selectNode", exc);
	}
	TreeNode result(nodes.begin()->second);
	nodes.erase(nodes.begin());
	return result;
}


int TreeNode::mostFractional() const {
	int index = -1;
	double smallestDisToMid = 0.5;
	for (int i = 0; i < valueInt.size(); ++i) {
		if (!isInteger(valueInt[i], PPM)) {
			double mid = (ceil(valueInt[i]) + floor(valueInt[i])) / 2;
			double distToMid = abs(valueInt[i] - mid);
			if (greaterThanReal(smallestDisToMid, distToMid, PPM)) {
				smallestDisToMid = distToMid;
				index = i;
			}
		}
	}
	return index;
}


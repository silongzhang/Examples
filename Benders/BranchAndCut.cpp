#include"DataStructure.h"

// Implement the benders decomposition with branch-and-cut.


TreeNode::TreeNode(const Instance& instance) :depth(1), valueEta(InfinityNeg), feasibleLP(false), objective(InfinityPos), integral(false) {
	branchCons.resize(instance.varInt.size);

	for (int i = 0; i < instance.varInt.size; ++i) {
		branchCons[i].first = 0;
		const double ub(instance.varInt.upperbounds[i]);
		if (equalToReal(ub, INFUB, PPM))
			branchCons[i].second = INFUB;
		else
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


// A naive branch strategy. Branch on the most fractional variable.
void Tree::branch(const TreeNode& node) {
	try {
		if (node.integral) throw exception();

		const int index = node.mostFractional();
		const double num = node.valueInt[index];
		const int leftHigh = floor(num), rightLow = ceil(num);

		TreeNode left(node), right(node);
		++left.depth, ++right.depth;
		left.branchCons[index].second = leftHigh;
		right.branchCons[index].first = rightLow;
		nodes.insert(make_pair(node.objective, left));
		nodes.insert(make_pair(node.objective, right));

		cout << "Branch on variable X[" << index << "]: " << '\t' << "[" << node.branchCons[index].first << ", " << leftHigh << "]" << '\t' << "(" << num << ")" << '\t' << "[" << rightLow << ", " << node.branchCons[index].second << "]" << endl;
	}
	catch (const exception& exc) {
		printErrorAndExit("Tree::branch", exc);
	}
}


void TreeNode::setBounds(IloEnv env, IloNumVarArray X) const {
	try {
		if (X.getSize() != branchCons.size()) throw exception();

		for (int i = 0; i < X.getSize(); ++i) {
			X[i].setLB(branchCons[i].first);
			X[i].setUB(equalToReal(branchCons[i].second, INFUB, PPM) ? IloInfinity : branchCons[i].second);
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("TreeNode::setBounds", exc);
	}
}


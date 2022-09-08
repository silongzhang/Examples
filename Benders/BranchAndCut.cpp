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


void TreeNode::solve(IloCplex cplexRMP, IloModel modelRMP, IloNumVarArray X, IloNumVar eta, IloCplex cplexSP, IloNumVarArray Y, IloRangeArray consSP, Solution& incumbent, const Instance& instance) {
	try {
		IloEnv env = cplexRMP.getEnv();
		setBounds(env, X);												// Set branching constraints.

		while (true) {
			solveModel(cplexRMP);										// Solve the LP relaxation of the RMP.
			if (cplexRMP.getStatus() == IloAlgorithm::Status::Infeasible) {
				feasibleLP = false;
				break;
			}
			else if (cplexRMP.getStatus() != IloAlgorithm::Status::Optimal)
				throw exception();
			else {
				valueInt = getValues(cplexRMP, X);
				valueEta = cplexRMP.getValue(eta);
				feasibleLP = true;
				objective = cplexRMP.getObjValue();

				if (!isInteger(valueInt, PPM)) break;
				else {
					// Set the right hand side of constraints in the SP.
					for (int i = 0; i < consSP.getSize(); ++i) {
						double rhs = instance.cpCons[i].rhs - inner_product(instance.cpCons[i].coefInt.begin(), instance.cpCons[i].coefInt.end(), valueInt.begin(), 0.0);
						setRhs(consSP[i], instance.cpCons[i].type, rhs);
					}

					solveModel(cplexSP);								// Solve the subproblem.
					IloNumArray dualSP(getDuals(cplexSP, consSP));		// Get dual values.

					if (cplexSP.getStatus() == IloAlgorithm::Status::Infeasible) {
						modelRMP.add(0 >= instance.exprRhs(env, dualSP, X));			// Add feasibility cut.
						++incumbent.nFeasCutIP;
					}
					else if (cplexSP.getStatus() == IloAlgorithm::Status::Optimal) {
						if (lessThanReal(valueEta, cplexSP.getObjValue(), PPM)) {
							modelRMP.add(eta >= instance.exprRhs(env, dualSP, X));		// Add optimality cut.
							++incumbent.nOptCutIP;
						}
						else if (equalToReal(valueEta, cplexSP.getObjValue(), PPM)) {
							valueCont = getValues(cplexSP, Y);
							integral = true;
							break;
						}
						else throw exception();
					}
					else throw exception();
				}
			}
		}
	}
	catch (const exception& exc) {
		printErrorAndExit("TreeNode::solve", exc);
	}
}


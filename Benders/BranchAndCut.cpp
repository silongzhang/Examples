#include"DataStructure.h"

// Implement the benders decomposition with branch-and-cut.
// A naive implementation, which can be improved by valid inequalities, parallel search, dynamic search, etc.


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
		nNodesGenerated += 2;

		//cout << "Branch on variable X[" << index << "]: " << '\t' << "[" << node.branchCons[index].first << ", " << leftHigh << "]" << '\t' << "(" << num << ")" << '\t' << "[" << rightLow << ", " << node.branchCons[index].second << "]" << endl;
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
			cplexRMP.solve();											// Solve the LP relaxation of the RMP.
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


// Assumption 1. The MILP model has a finite optimal objective value.
// Assumption 2. # of integer variables > 0 and # of continuous variables > 0.
// Assumption 3. Variables are already transformed to be nonnegative.
// Assumption 4. Invoke the method "standardize" in advance if the instance is not "standard".
// Assumption 5. It gives the opposite number of the optimal objective value of the original instance if it is a maximization problem.
// Assumption 6. Relaxed master problems (RMPs) cannot be unbounded.
// Note 1. For a minimization problem, a sufficient condition for Assumption 6 is that integer variables are all nonnegative and their coefficients in the objective function are all nonnegative. A similar sufficient condition can be derived for a maximization problem.
// Note 2. Another sufficient condition for Assumption 6 is that the number of feasible combinations of values of integer variables is finite. For example, each integer variable is bounded from both below and above.
Solution Instance::solveBendersBC(const ParameterAlgorithm& parameter) const {
	Solution incumbent;
	incumbent.clear();
	IloEnv env;
	try {
		const clock_t start = clock();
		cout << "Running Instance::solveBendersBC ..." << endl;
		if (!standard()) throw exception();
		const int NInt = varInt.size, NCont = varCont.size;

		// Declare modeling components.
		IloNumVar eta(env, -IloInfinity, IloInfinity);
		IloNumVarArray X(env), Y(env, NCont, 0, IloInfinity);
		IloModel modelRMP(env), modelSP(env);
		vector<double> currentValInt(NInt, 0);
		double currentValEta = InfinityNeg;
		IloRangeArray consSP(env);
		IloNumArray dualSP(env);

		// Initiate models.
		initiateModels(env, modelRMP, modelSP, X, Y, eta, currentValInt, consSP);
		IloCplex cplexRMP(modelRMP), cplexSP(modelSP);
		cplexSP.setParam(IloCplex::RootAlg, IloCplex::Primal);

		cout << "#############################" << endl;
		cout << "Add Benders cuts based on solutions of the LP relaxation of the RMP." << endl;
		bool proceed = true;
		for (int iter = 1; true; ++iter) {
			cout << "+++++++++++++++++++++++++++++" << endl;
			cout << "Iter = " << iter << '\t' << "LB = " << incumbent.LB << '\t' << "nOptCutLP = " << incumbent.nOptCutLP << '\t' << "nFeasCutLP = " << incumbent.nFeasCutLP << endl;

			// Solve the LP relaxation of the RMP.
			proceed = solveLRMP(cplexRMP, X, eta, incumbent, currentValInt, currentValEta);
			if (!proceed) break;

			// Solve the SP.
			proceed = solveSP(cplexSP, consSP, dualSP, incumbent, cpCons, currentValInt);
			if (!proceed) break;

			// Add Benders cuts.
			if (!addBendersCuts(env, cplexSP, modelRMP, X, eta, dualSP, incumbent, *this, currentValEta, false)) break;
		}

		if (proceed) {
			cplexRMP.setOut(env.getNullStream());
			cplexSP.setOut(env.getNullStream());
			int nInfeas = 0, nBound = 0, nIt = 0, nBranch = 0;
			Tree tree;
			tree.nodes.insert(make_pair(incumbent.LB, TreeNode(*this)));			// The root node.
			tree.nNodesGenerated = 1;
			for (int iter = 1; !tree.nodes.empty() && !parameter.stop(start, incumbent.objective, incumbent.LB); ++iter) {
				incumbent.LB = min(tree.nodes.begin()->first, incumbent.objective);	// Renew LB.

				if (iter % parameter.iterPrintBC == 0) {
					cout << "+++++++++++++++++++++++++++++" << endl;
					cout << "Iter = " << iter << '\t' << "Time = " << runTime(start) << '\t' << "UB = " << incumbent.objective << '\t' << "LB = " << incumbent.LB << '\t' << "# of nodes generated = " << tree.nNodesGenerated << '\t' << "# of nodes remained = " << tree.nodes.size() << '\t' << "nInfeas = " << nInfeas << '\t' << "nBound = " << nBound << '\t' << "nIt = " << nIt << '\t' << "nBranch = " << nBranch << endl;
				}

				TreeNode node(tree.selectNode());									// Node selection.
				if (greaterThanReal(node.objective, incumbent.objective, PPM))
					++nBound;														// Pruned due to boundedness.
				else {
					node.solve(cplexRMP, modelRMP, X, eta, cplexSP, Y, consSP, incumbent, *this);		// Solve the LP relaxation.
					if (!node.feasibleLP)
						++nInfeas;														// Pruned due to infeasibility.
					else if (greaterThanReal(node.objective, incumbent.objective, PPM))
						++nBound;														// Pruned due to boundedness.
					else if (node.integral) {
						++nIt;															// Pruned due to integrality.
						if (greaterThanReal(incumbent.objective, node.objective, PPM)) {
							incumbent.status = SolutionStatus::Feasible;
							incumbent.objective = node.objective;						// Renew UB.
							incumbent.valueInt = node.valueInt;
							incumbent.valueEta = node.valueEta;
							incumbent.valueCont = node.valueCont;
							cout << "Found a better solution with objective = " << incumbent.objective << endl;
						}
					}
					else {
						++nBranch;														// Branched.
						tree.branch(node);
					}
				}
			}

			if (!tree.nodes.empty()) {
				incumbent.LB = min(tree.nodes.begin()->first, incumbent.objective);
				if (equalToReal(incumbent.objective, incumbent.LB, PPM))
					incumbent.status = SolutionStatus::Optimal;
			}
			else if (incumbent.status == SolutionStatus::Feasible) {
				incumbent.LB = incumbent.objective;
				incumbent.status = SolutionStatus::Optimal;
			}
			else incumbent.status = SolutionStatus::Infeasible;

			cout << "UB = " << incumbent.objective << '\t' << "LB = " << incumbent.LB << '\t' << "# of nodes generated = " << tree.nNodesGenerated << '\t' << "# of nodes remained = " << tree.nodes.size() << '\t' << "nInfeas = " << nInfeas << '\t' << "nBound = " << nBound << '\t' << "nIt = " << nIt << '\t' << "nBranch = " << nBranch << endl;
		}
		incumbent.elapsedTime = runTime(start);
		cout << "elapsed time (Instance::solveBendersBC): " << runTime(start) << endl;
		incumbent.print();
	}
	catch (const exception& exc) {
		printErrorAndExit("Instance::solveBendersBC", exc);
	}

	env.end();
	return incumbent;
}


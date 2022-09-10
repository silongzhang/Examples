C++ (Visual studio 2019), Cplex 12.10.

Naive implementations of the classical benders decomposition algorithm for solving general mixed integer linear programming (MILP) models.

Step 1. Add Benders cuts based on solutions of LP relaxations of relaxed master problems (RMPs), until the LP relaxation of the master problem (MP) is solved to optimality.
Step 2. Close the integrality gap, where during the procedure, (i) the current RMP and subproblem (SP) is solved to optimality alternately (Instance::solveBendersRecursive), or (ii) use legacy callback to add lazy constraints (Instance::solveBendersLegacyCallback), or (iii) use generic callback to add lazy constraints (Instance::solveBendersGenericCallback), or (iv) implement the individualized branch-and-cut algorithm (Instance::solveBendersBC).

Assumption 1. The MILP model has a finite optimal objective value.
Assumption 2. # of integer variables > 0 and # of continuous variables > 0.
Assumption 3. Variables are already transformed to be nonnegative.
Assumption 4. Invoke the method "standardize" in advance if the instance is not "standard".
Assumption 5. It gives the opposite number of the optimal objective value of the original instance if it is a maximization problem.
Assumption 6. Relaxed master problems (RMPs) cannot be unbounded.
Note 1. For a minimization problem, a sufficient condition for Assumption 6 is that integer variables are all nonnegative and their coefficients in the objective function are all nonnegative. A similar sufficient condition can be derived for a maximization problem.
Note 2. Another sufficient condition for Assumption 6 is that the number of feasible combinations of values of integer variables is finite. For example, each integer variable is bounded from both below and above.

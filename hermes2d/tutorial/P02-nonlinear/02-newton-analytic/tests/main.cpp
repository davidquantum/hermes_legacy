#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This test makes sure that example 17-newton-elliptic-2 works correctly.

const int P_INIT = 2;                             // Initial polynomial degree.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 9;                    // Maximum allowed number of Newton iterations.
const int INIT_GLOB_REF_NUM = 3;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 4;                   // Number of initial refinements towards boundary.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double alpha = 4.0;
double heat_src = 1.0;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  HermesFunction src(-heat_src);
  WeakFormsH1::DefaultWeakFormPoisson wf(HERMES_ANY, &lambda, &src);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the solution.
  Solution sln;

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)] ;
  CustomInitialCondition init_sln(&mesh);
  OGProjection::project_global(&space, &init_sln, coeff_vec, matrix_solver);

  // Perform Newton's iteration.
  bool verbose = true;
  bool jacobian_changed = true;
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
      NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // Cleanup.
  delete []coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  info("ndof = %d", ndof);
  info("Coordinate (1, 0) value = %lf", sln.get_pt_value(1.0, 0.0));
  info("Coordinate (3, 0) value = %lf", sln.get_pt_value(3.0, 0.0));
  info("Coordinate (5, 0) value = %lf", sln.get_pt_value(5.0, 0.0));
  info("Coordinate (7, 0) value = %lf", sln.get_pt_value(7.0, 0.0));

  double coor_x[4] = {1.0, 3.0, 5.0, 7.0};
  double coor_y = 0.0;
  double t_value[4] = {2.867436, 2.873677, 2.832594, 2.709390};
  bool success = true;

  for (int i = 0; i < 4; i++)
  {
    if (abs(t_value[i] - sln.get_pt_value(coor_x[i], coor_y)) > 1E-6) success = false;
  }

  if (success) {  // should pass with NEWTON_MAX_ITER = 9 and fail with NEWTON_MAX_ITER = 8
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }


}


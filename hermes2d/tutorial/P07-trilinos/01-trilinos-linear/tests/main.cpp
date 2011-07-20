#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace Teuchos;

// This test makes sure that example 40-trilinos-linear works correctly.

const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const bool TRILINOS_JFNK = true;                  // true = Jacobian-free method (for NOX),
                                                  // false = Newton (for NOX).
const bool PRECOND = true;                        // Preconditioning by jacobian in case of JFNK (for NOX),
                                                  // default ML preconditioner in case of Newton.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)
// NOX parameters.
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::Details;
                                                  // Error messages, see NOX_Utils.h.

double ls_tolerance = 1e-5;                       // Tolerance for linear system.
unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-3;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-2;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.

int main(int argc, char **argv)
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set exact solution.
  CustomExactSolution exact(&mesh);
  
  // Initialize the weak formulation.
  CustomWeakFormPoisson wf1;

  // Initialize boundary conditions
  DefaultEssentialBCNonConst bc_essential("Bdy", &exact);
  EssentialBCs bcs(&bc_essential);
 
  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("---- Assembling by DiscreteProblem, solving by %s:", 
       MatrixSolverNames[matrix_solver].c_str());

  // Initialize the linear discrete problem.
  DiscreteProblem dp1(&wf1, &space);
    
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Begin time measurement of assembly.
  cpu_time.tick(HERMES_SKIP);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Perform Newton's iteration.
  if (!hermes2d.solve_newton(coeff_vec, &dp1, solver, matrix, rhs)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution sln1;
  Solution::vector_to_solution(coeff_vec, &space, &sln1);

  // CPU time measurement.
  double time = cpu_time.tick().last();

  // Calculate errors.
  double rel_err_1 = hermes2d.calc_rel_error(&sln1, &exact, HERMES_H1_NORM) * 100;
  info("CPU time: %g s.", time);
  info("Xxact H1 error: %g%%.", rel_err_1);

  delete(matrix);
  delete(rhs);
  delete(solver);
    
  // View the solution and mesh.
  //ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  //sview.show(&sln1);
  //OrderView  oview("Polynomial orders", new WinGeom(450, 0, 400, 350));
  //oview.show(&space);
  
  // TRILINOS PART:

  info("---- Assembling by DiscreteProblem, solving by NOX:");

  // Initialize the weak formulation for Trilinos.
  CustomWeakFormPoisson wf2(TRILINOS_JFNK);
  
  // Initialize DiscreteProblem.
  DiscreteProblem dp2(&wf2, &space);
  
  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Set initial vector for NOX.
  // NOTE: Using zero vector was causing convergence problems.
  info("Projecting to obtain initial vector for the Newton's method.");
  Solution init_sln(&mesh, 0.0);
  OGProjection::project_global(&space, &init_sln, coeff_vec);
  
  // Measure the projection time.
  double proj_time = cpu_time.tick().last();
  
  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  // "" stands for preconditioning that is set later.
  NoxSolver nox_solver(&dp2, message_type, "GMRES", "Newton", ls_tolerance, "", flag_absresid, abs_resid, 
                       flag_relresid, rel_resid, max_iters);
  nox_solver.set_init_sln(coeff_vec);
  
  delete coeff_vec;

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (TRILINOS_JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Assemble and solve using NOX.
  Solution sln2;
  if (nox_solver.solve())
  {
    Solution::vector_to_solution(nox_solver.get_solution(), &space, &sln2);

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else error("NOX failed");

  // CPU time needed by NOX.
  double time2 = cpu_time.tick().last();

  // Show the NOX solution.
  //ScalarView view2("Solution 2", new WinGeom(450, 0, 440, 350));
  //view2.show(&sln2);
  //view2.show(&exact);

  // Calculate errors.
  double rel_err_2 = hermes2d.calc_rel_error(&sln2, &exact, HERMES_H1_NORM) * 100;
  info("Projection time: %g s, NOX assembly/solution time: %g s.", proj_time, time2);
  info("Exact H1 error: %g%%.)", rel_err_2);
 

  /* TESTING */

  info("Coordinate (-0.6, -0.6) hermes value = %lf", sln1.get_pt_value(-0.6, -0.6));
  info("Coordinate ( 0.4, -0.6) hermes value = %lf", sln1.get_pt_value( 0.4, -0.6));
  info("Coordinate ( 0.4,  0.4) hermes value = %lf", sln1.get_pt_value( 0.4,  0.4));
  info("Coordinate (-0.6,  0.0) hermes value = %lf", sln1.get_pt_value(-0.6,  0.0));
  info("Coordinate ( 0.0,  0.0) hermes value = %lf", sln1.get_pt_value( 0.0,  0.0));

  info("Coordinate (-0.6, -0.6) nox value = %lf", sln2.get_pt_value(-0.6, -0.6));
  info("Coordinate ( 0.4, -0.6) nox value = %lf", sln2.get_pt_value( 0.4, -0.6));
  info("Coordinate ( 0.4,  0.4) nox value = %lf", sln2.get_pt_value( 0.4,  0.4));
  info("Coordinate (-0.6,  0.0) nox value = %lf", sln2.get_pt_value(-0.6,  0.0));
  info("Coordinate ( 0.0,  0.0) nox value = %lf", sln2.get_pt_value( 0.0,  0.0));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(sln2.get_pt_value(-0.6, -0.6) - 0.720000) > eps) {
    printf("Coordinate (-0.6, -0.6) nox value is %g\n", sln2.get_pt_value(-0.6, -0.6));
    success = 0;
  }
  if (fabs(sln2.get_pt_value( 0.4, -0.6) - 0.520000) > eps) {
    printf("Coordinate ( 0.4, -0.6) nox value is %g\n", sln2.get_pt_value( 0.4, -0.6));
    success = 0;
  }
  if (fabs(sln2.get_pt_value( 0.4,  0.4) - 0.320000) > eps) {
    printf("Coordinate ( 0.4,  0.4) nox value is %g\n", sln2.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln2.get_pt_value(-0.6,  0.0) - 0.360000) > eps) {
    printf("Coordinate (-0.6,  0.0) nox value is %g\n", sln2.get_pt_value(-0.6,  0.0));
    success = 0;
  }
  if (fabs(sln2.get_pt_value( 0.0,  0.0) - 0.000000) > eps) {
    printf("Coordinate ( 0.0,  0.0) nox value is %g\n", sln2.get_pt_value( 0.0,  0.0));
    success = 0;
  }

  if (fabs(sln1.get_pt_value(-0.6, -0.6) - 0.720000) > eps) {
    printf("Coordinate (-0.6, -0.6) hermes value is %g\n", sln1.get_pt_value(-0.6, -0.6));
    success = 0;
  }
  if (fabs(sln1.get_pt_value( 0.4, -0.6) - 0.520000) > eps) {
    printf("Coordinate ( 0.4, -0.6) hermes value is %g\n", sln1.get_pt_value( 0.4, -0.6));
    success = 0;
  }
  if (fabs(sln1.get_pt_value( 0.4,  0.4) - 0.320000) > eps) {
    printf("Coordinate ( 0.4,  0.4) hermes value is %g\n", sln1.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln1.get_pt_value(-0.6,  0.0) - 0.360000) > eps) {
    printf("Coordinate (-0.6,  0.0) hermes value is %g\n", sln1.get_pt_value(-0.6,  0.0));
    success = 0;
  }
  if (fabs(sln1.get_pt_value( 0.0,  0.0) - 0.000000) > eps) {
    printf("Coordinate ( 0.0,  0.0) hermes value is %g\n", sln1.get_pt_value( 0.0,  0.0));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

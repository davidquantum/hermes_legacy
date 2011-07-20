#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace Teuchos;
using namespace RefinementSelectors;

// This test makes sure that example 42-trilinos-timedep works correctly.

const int INIT_REF_NUM = 4;       // Number of initial uniform mesh refinements.
const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
const double ALPHA = 10.0;        // Coefficient for the Nwwton boundary condition.
const double LAMBDA = 1e5;
const double HEATCAP = 1e6;
const double RHO = 3000.0;
const double TEMP_EXT = 20.0;
const double TEMP_INIT = 10.0;

const double TAU = 50.0;          // Time step.        

const bool JFNK = true;
const bool PRECOND = true;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc("Bdy_bottom", TEMP_INIT);
  EssentialBCs bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  // Define constant initial condition. 
  Solution t_prev_time(&mesh, TEMP_INIT);

  // Initialize the weak formulation.
  CustomWeakForm wf(Hermes::vector<std::string>("Bdy_right", "Bdy_top", "Bdy_left"), HEATCAP, RHO, TAU, LAMBDA, ALPHA, TEMP_EXT, &t_prev_time, JFNK);

  // Initialize the finite element problem.
  DiscreteProblem dp(&wf, &space);

  // Project the function "t_prev_time" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(&space, &t_prev_time, coeff_vec);

  // Initialize NOX solver.
  NoxSolver solver(&dp);

  // Select preconditioner.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) solver.set_precond(pc);
    else solver.set_precond("ML");
  }

  // Time stepping loop:
  double total_time = 0.0;
  for (int ts = 1; total_time <= 2000.0; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time += TAU);

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
      Solution::vector_to_solution(solver.get_solution(), &space, &t_prev_time);
    else
      error("NOX failed.");

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }

  info("Coordinate ( 0.6,  0.6) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.6,  0.6));
  info("Coordinate ( 0.4,  0.6) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.4,  0.6));
  info("Coordinate ( 0.4,  0.4) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.4,  0.4));
  info("Coordinate ( 0.6,  0.0) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.6,  0.0));
  info("Coordinate ( 0.5,  0.5) t_prev_time value = %lf", t_prev_time.get_pt_value( 0.5,  0.5));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(t_prev_time.get_pt_value(0.6, 0.6) - 13.826133) > eps) {
    printf("Coordinate (0.6, 0.6) t_prev_time value is %g\n", t_prev_time.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.4, 0.6) - 13.826133) > eps) {
    printf("Coordinate ( 0.4, 0.6) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.4,  0.4) - 12.599711) > eps) {
    printf("Coordinate ( 0.4,  0.4) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value(0.6,  0.0) - 10.000000) > eps) {
    printf("Coordinate (0.6,  0.0) t_prev_time value is %g\n", t_prev_time.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if (fabs(t_prev_time.get_pt_value( 0.5,  0.5) - 12.916592) > eps) {
    printf("Coordinate ( 0.5,  0.5) t_prev_time value is %g\n", t_prev_time.get_pt_value( 0.5,  0.5));
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

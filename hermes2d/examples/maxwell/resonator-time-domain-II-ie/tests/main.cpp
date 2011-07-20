#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

// This test makes sure that the example "maxwell/resonator-time-domain-II-ie" works correctly.

// The following parameters can be changed:
const int P_INIT = 6;                              // Initial polynomial degree of all elements.
const int INIT_REF_NUM = 1;                        // Number of initial uniform mesh refinements.
const double time_step = 0.05;                     // Time step.
const double T_FINAL = time_step*5;                // Final time.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,

// Problem parameters.
const double C_SQUARED = 1;                        // Square of wave speed.                     

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  CustomInitialConditionWave E_sln(&mesh);
  Solution F_sln(&mesh, 0.0, 0.0);
  Hermes::vector<Solution*> slns(&E_sln, &F_sln);

  // Initialize the weak formulation.
  CustomWeakFormWave wf(time_step, C_SQUARED, &E_sln, &F_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential("Perfect conductor", 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace E_space(&mesh, &bcs, P_INIT);
  HcurlSpace F_space(&mesh, &bcs, P_INIT);
  Hermes::vector<Space *> spaces = Hermes::vector<Space *>(&E_space, &F_space);

  info("ndof = %d.", Space::get_num_dofs(spaces));

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, spaces);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

  // Initialize views.
  // ScalarView E1_view("Solution E1", new WinGeom(0, 0, 400, 350));
  // E1_view.fix_scale_width(50);
  // ScalarView E2_view("Solution E2", new WinGeom(410, 0, 400, 350));
  // E2_view.fix_scale_width(50);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Perform one implicit Euler time step.
    info("Implicit Euler time step (t = %g s, time_step = %g s).", current_time, time_step);

    // First time assemble both the stiffness matrix and right-hand side vector,
    // then just the right-hand side vector.
    if (ts == 1) {
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);
    }
    else {
      info("Assembling the right-hand side vector (only).");
      dp.assemble(NULL, rhs);
    }

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), spaces, slns);
    else error ("Matrix solver failed.\n");

    /*
    // Visualize the solutions.
    char title[100];
    sprintf(title, "E1, t = %g", current_time);
    E1_view.set_title(title);
    E1_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    sprintf(title, "E2, t = %g", current_time);
    E2_view.set_title(title);
    E2_view.show(&E_sln, HERMES_EPS_NORMAL, H2D_FN_VAL_1);
    */

    // Update time.
    current_time += time_step;
  } while (current_time < T_FINAL);

  double coord_x[4] = {0.3, 0.6, 0.9, 1.4};
  double coord_y[4] = {0, 0.3, 0.5, 0.7};

  info("Coordinate (0.3, 0.0) value = %lf", F_sln.get_pt_value(coord_x[0], coord_y[0]));
  info("Coordinate (0.6, 0.3) value = %lf", F_sln.get_pt_value(coord_x[1], coord_y[1]));
  info("Coordinate (0.9, 0.5) value = %lf", F_sln.get_pt_value(coord_x[2], coord_y[2]));
  info("Coordinate (1.4, 0.7) value = %lf", F_sln.get_pt_value(coord_x[3], coord_y[3]));

  double t_value[4] = {-0.142680, -0.260440, -0.331901, -0.363901};
  bool success = true;

  for (int i = 0; i < 4; i++)
  {
    if (fabs(t_value[i] - F_sln.get_pt_value(coord_x[i], coord_y[i])) > 1E-6) success = false;
  }

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

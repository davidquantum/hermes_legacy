#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace RefinementSelectors;

// This test makes sure that example 19-newton-timedep-heat-rk works correctly.

const int INIT_GLOB_REF_NUM = 3;                   // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                    // Number of initial refinements towards boundary.
const int P_INIT = 2;                              // Initial polynomial degree.
const double time_step = 0.001;                    // Time step.
const double T_FINAL = 2*time_step;                // Time interval length.
const double NEWTON_TOL = 1e-5;                    // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2,
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4,
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded,
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded,
//   Implicit_DIRK_ISMAIL_7_45_embedded.
ButcherTableType butcher_table_type = Implicit_RK_1;

// Problem parameters.
const double alpha = 4.0;                         // For the nonlinear thermal conductivity.
const double heat_src = 1.0;

// Main function.
int main(int argc, char* argv[])
{
  // Check number of command-line parameters.
  if (argc < 2)
    error("Not enough parameters: Provide a Butcher's table type.");

  int b_type = atoi(argv[1]);
  info ("%d", b_type);

  switch (b_type)
  {
    case 1: butcher_table_type = Explicit_RK_1; break;
    case 2: butcher_table_type = Implicit_RK_1; break;
    case 3: butcher_table_type = Explicit_RK_2; break;
    case 4: butcher_table_type = Implicit_Crank_Nicolson_2_2; break;
    case 5: butcher_table_type = Implicit_SDIRK_2_2; break;
    case 6: butcher_table_type = Implicit_Lobatto_IIIA_2_2; break;
    case 7: butcher_table_type = Implicit_Lobatto_IIIB_2_2; break;
    case 8: butcher_table_type = Implicit_Lobatto_IIIC_2_2; break;
    case 9: butcher_table_type = Explicit_RK_3; break;
    case 10: butcher_table_type = Explicit_RK_4; break;
    case 11: butcher_table_type = Implicit_Lobatto_IIIA_3_4; break;
    case 12: butcher_table_type = Implicit_Lobatto_IIIB_3_4; break;
    case 13: butcher_table_type = Implicit_Lobatto_IIIC_3_4; break;
    case 14: butcher_table_type = Implicit_Radau_IIA_3_5; break;
    case 15: butcher_table_type = Implicit_SIRK_2_2; break;
    case 16: butcher_table_type = Implicit_ESIRK_2_2; break;
    case 17: butcher_table_type = Implicit_SDIRK_5_4; break;

    default: error("Admissible command-line options are from 1 to 17.");
  }

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d.", ndof);

  // Previous time level solution (initialized by the initial condition).
  CustomInitialCondition sln_time_prev(&mesh);

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  HermesFunction f(heat_src);
  WeakFormsH1::DefaultWeakFormPoisson wf(HERMES_ANY, &lambda, &f);

  // Previous and next time level solutions.
  Solution sln_time_new(&mesh);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g, tau = %g, stages: %d).",
         current_time, time_step, bt.get_size());
    bool verbose = true;
    Hermes::vector<Solution*> slns_time_prev;
    slns_time_prev.push_back(&sln_time_prev);
    Hermes::vector<Solution*> slns_time_new;
    slns_time_new.push_back(&sln_time_new);

    if (!runge_kutta.rk_time_step(current_time, time_step, slns_time_prev, slns_time_new, true, verbose, NEWTON_TOL, NEWTON_MAX_ITER))
      error("Runge-Kutta time step failed, try to decrease time step size.");

    // Update time.
    current_time += time_step;

    // Copy solution for the new time step.
    sln_time_prev.copy(&sln_time_new);
    // Increase counter of time steps.
    ts++;
  }
  while (current_time < T_FINAL);

  info("Coordinate (-8.0, -8.0) value = %lf", sln_time_new.get_pt_value(-8.0, -8.0));
  info("Coordinate (-5.0, -5.0) value = %lf", sln_time_new.get_pt_value(-5.0, -5.0));
  info("Coordinate (-3.0, -3.0) value = %lf", sln_time_new.get_pt_value(-3.0, -3.0));
  info("Coordinate ( 0.0,  0.0) value = %lf", sln_time_new.get_pt_value(0.0,  0.0));
  info("Coordinate ( 3.0,  3.0) value = %lf", sln_time_new.get_pt_value(3.0,  3.0));
  info("Coordinate ( 5.0,  5.0) value = %lf", sln_time_new.get_pt_value(5.0,  5.0));
  info("Coordinate ( 8.0,  8.0) value = %lf", sln_time_new.get_pt_value(8.0,  8.0));

  double coor_x_y[7] = {-8.0, -5.0, -3.0, 0.0, 3.0, 5.0, 8.0};
  bool success = true;

  switch (b_type)
  {
    case 1: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042560) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492025) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002150) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693363) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255564) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.266989) > 1E-6) success = false;
            break;

    case 2: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251887) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002152) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693351) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255921) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.264214) > 1E-6) success = false;
            break;

    case 3: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251887) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693354) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255822) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265515) > 1E-6) success = false;
            break;

    case 4: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693356) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255784) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265491) > 1E-6) success = false;
            break;

    case 5: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255791) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265487) > 1E-6) success = false;
            break;

    case 6: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693356) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255784) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265491) > 1E-6) success = false;
            break;

    case 7: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693356) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255784) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265491) > 1E-6) success = false;
            break;

    case 8: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251887) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693354) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255819) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265435) > 1E-6) success = false;
            break;

    case 9: if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255795) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265458) > 1E-6) success = false;
            break;

    case 10:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255797) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265487) > 1E-6) success = false;
            break;

    case 11:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
            break;

    case 12:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
            break;

    case 13:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
            break;

    case 14:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
            break;

    case 15:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255791) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265486) > 1E-6) success = false;
            break;

    case 16:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255791) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265484) > 1E-6) success = false;
            break;

    case 17:if (fabs(sln_time_new.get_pt_value(coor_x_y[0], coor_x_y[0]) - 0.042559) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[1], coor_x_y[1]) - 0.251886) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[2], coor_x_y[2]) - 0.492024) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[3], coor_x_y[3]) - 1.002151) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[4], coor_x_y[4]) - 1.693355) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[5], coor_x_y[5]) - 2.255798) > 1E-6) success = false;
            if (fabs(sln_time_new.get_pt_value(coor_x_y[6], coor_x_y[6]) - 3.265482) > 1E-6) success = false;
            break;

    default: error("Admissible command-line options are from 1 to 17.");
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

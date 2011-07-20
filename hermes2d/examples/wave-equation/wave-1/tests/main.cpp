#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

// This test makes sure that the example "wave-equation/wave-1" works correctly.

// The following parameters can be changed:
const int P_INIT = 6;                              // Initial polynomial degree of all elements.
const double time_step = 0.01;                     // Time step.
const double T_FINAL = time_step*10;               // Final time.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,

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
const double C_SQUARED = 100;                      // Square of wave speed.                     

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);
  if (bt.is_explicit()) info("Using a %d-stage explicit R-K method.", bt.get_size());
  if (bt.is_diagonally_implicit()) info("Using a %d-stage diagonally implicit R-K method.", bt.get_size());
  if (bt.is_fully_implicit()) info("Using a %d-stage fully implicit R-K method.", bt.get_size());

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Convert to quadrilaterals.
  mesh.convert_triangles_to_quads();

  // Refine towards boundary.
  mesh.refine_towards_boundary("Bdy", 1, true);

  // Refine once towards vertex #4.
  mesh.refine_towards_vertex(4, 1);

  // Initialize solutions.
  CustomInitialConditionWave u_sln(&mesh);
  Solution v_sln(&mesh, 0.0);
  Hermes::vector<Solution*> slns(&u_sln, &v_sln);

  // Initialize the weak formulation.
  CustomWeakFormWave wf(time_step, C_SQUARED, &u_sln, &v_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential("Bdy", 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u_space(&mesh, &bcs, P_INIT);
  H1Space v_space(&mesh, &bcs, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)));

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u_space, &v_space));

  // Initialize Runge-Kutta time stepping.
  RungeKutta runge_kutta(&dp, &bt, matrix_solver);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    info("Runge-Kutta time step (t = %g s, time_step = %g s, stages: %d).", 
         current_time, time_step, bt.get_size());
    bool jacobian_changed = false;
    bool verbose = true;

    if (!runge_kutta.rk_time_step(current_time, time_step, slns, 
                                  slns, jacobian_changed, verbose))
      error("Runge-Kutta time step failed, try to decrease time step size.");

    // Update time.
    current_time += time_step;

  } while (current_time < T_FINAL);

  double coord_x[4] = {1, 3, 5, 7};
  double coord_y[4] = {1, 3, 5, 7};

  info("Coordinate (1.0, 0.0) value = %lf", u_sln.get_pt_value(coord_x[0], coord_y[0]));
  info("Coordinate (3.0, 3.0) value = %lf", u_sln.get_pt_value(coord_x[1], coord_y[1]));
  info("Coordinate (5.0, 5.0) value = %lf", u_sln.get_pt_value(coord_x[2], coord_y[2]));
  info("Coordinate (7.0, 7.0) value = %lf", u_sln.get_pt_value(coord_x[3], coord_y[3]));

  info("Coordinate (1.0, 0.0) value = %lf", v_sln.get_pt_value(coord_x[0], coord_y[0]));
  info("Coordinate (3.0, 3.0) value = %lf", v_sln.get_pt_value(coord_x[1], coord_y[1]));
  info("Coordinate (5.0, 5.0) value = %lf", v_sln.get_pt_value(coord_x[2], coord_y[2]));
  info("Coordinate (7.0, 7.0) value = %lf", v_sln.get_pt_value(coord_x[3], coord_y[3]));

  double t_value[8] = {0.212655, 0.000163, 0.000000, 0.000000, -0.793316, 0.007255, 0.000001, 0.000000};
  bool success = true;

  for (int i = 0; i < 4; i++)
  {
    if (fabs(t_value[i] - u_sln.get_pt_value(coord_x[i], coord_y[i])) > 1E-6) success = false;
  }

  for (int i = 4; i < 8; i++)
  {
    if (fabs(t_value[i] - v_sln.get_pt_value(coord_x[i-4], coord_y[i-4])) > 1E-6) success = false;
  }

  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }

  // Wait for the view to be closed.
  View::wait();

  return 0;
}

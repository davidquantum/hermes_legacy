#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves the Rayleigh-Benard convection problem
// http://en.wikipedia.org/wiki/Rayleigh%E2%80%93B%C3%A9nard_convection.
// In this problem, a steady fluid is heated from the bottom and 
// it starts to move. The time-dependent laminar incompressible Navier-Stokes 
// equations are coupled with a heat transfer equation and discretized 
// in time via the implicit Euler method. The Newton's method is used 
// to solve the nonlinear problem at each time step. Flow pressure can be 
// approximated using either continuous (H1) elements or discontinuous (L2)
// elements. The L2 elements for pressure make the velocity dicreetely
// divergence-free. 
//
// PDE: incompressible Navier-Stokes equations in the form
//      \partial v / \partial t = \Delta v / Pr - (v \cdot \nabla) v - \nabla p - Ra(T)Pr(0, -T),
//      div v = 0,
//      \partial T / \partial t = -v \cdot \nabla T + \Delta T.
//
// BC: velocity... zero on the entire boundary,
//     temperature... constant on the bottom,
//                    zero Neumann on vetrical edges,
//                    Newton (heat loss) (1 / Pr) * du/dn = ALPHA_AIR * (TEMP_EXT - u) on the top edge.
//
// Geometry: Rectangle (0, Lx) x (0, Ly)... see the file domain.mesh.
//
// The following parameters can be changed:

#define PRESSURE_IN_L2                            // If this is defined, the pressure is approximated using
                                                  // discontinuous L2 elements (making the velocity discreetely
                                                  // divergence-free, more accurate than using a continuous
                                                  // pressure approximation). Otherwise the standard continuous
                                                  // elements are used. The results are striking - check the
                                                  // tutorial for comparisons.
const int P_INIT_VEL = 2;                         // Initial polynomial degree for velocity components.
const int P_INIT_PRESSURE = 1;                    // Initial polynomial degree for pressure.
                                                  // Note: P_INIT_VEL should always be greater than
                                                  // P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_TEMP = 1;                        // Initial polynomial degree for temperature
const double time_step = 0.1;                     // Time step.
const double T_FINAL = 3600.0;                    // Time interval length.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double Pr = 7.0;                            // Prandtl number (water around 20 degrees Celsius).
const double Ra = 50;                             // Rayleigh number.
const double TEMP_INIT = 20;
const double TEMP_BOTTOM = 25;
const double TEMP_EXT = 20;                       // External temperature above the surface of the water.
const double ALPHA_AIR = 5.0;                     // Heat transfer coefficient between water and air on top edge.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  Hermes2D hermes_2D;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < 3; i++) mesh.refine_all_elements(2);
  for (int i=0; i < 3; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(HERMES_ANY, 2);

  // Initialize boundary conditions.
  DefaultEssentialBCConst zero_vel_bc(Hermes::vector<std::string>("Bottom", "Right", "Top", "Left"), 0.0);
  EssentialBCs bcs_vel_x(&zero_vel_bc);
  EssentialBCs bcs_vel_y(&zero_vel_bc);
  DefaultEssentialBCConst bc_temp_bottom("Bottom", TEMP_BOTTOM);
  EssentialBCs bcs_temp(&bc_temp_bottom);
  EssentialBCs bcs_pressure;

  // Spaces for velocity components and pressure.
  H1Space xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#endif
  H1Space t_space(&mesh, &bcs_temp, P_INIT_TEMP);
  Hermes::vector<Space *> spaces = Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space, &t_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space, &t_space));
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif
  ProjNormType t_proj_norm = HERMES_H1_NORM;

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  Solution xvel_prev_time, yvel_prev_time, p_prev_time, t_prev_time; 
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);
  t_prev_time.set_const(&mesh, TEMP_INIT);
  Hermes::vector<Solution*> slns = Hermes::vector<Solution*>(&xvel_prev_time, &yvel_prev_time, 
                                                             &p_prev_time, &t_prev_time);

  // Initialize weak formulation.
  WeakForm* wf = new WeakFormRayleighBenard(Pr, Ra, "Top", TEMP_EXT, ALPHA_AIR, time_step, 
                                            &xvel_prev_time, &yvel_prev_time, &t_prev_time);

  // Initialize the FE problem.
  DiscreteProblem dp(wf, spaces);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize views.
  VectorView vview("velocity", new WinGeom(0, 0, 1000, 200));
  ScalarView tview("temperature", new WinGeom(0, 255, 1000, 200));
  ScalarView pview("pressure", new WinGeom(0, 485, 1000, 200));
  //vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  pview.fix_scale_width(80);
  tview.fix_scale_width(80);
  pview.show_mesh(true);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  OGProjection::project_global(spaces, slns, coeff_vec, matrix_solver, 
                               Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm, t_proj_norm));

  // Time-stepping loop:
  char title[100];
  double current_time = 0;
  int num_time_steps = T_FINAL / time_step;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d, time = %g:", ts, current_time);

    // Perform Newton's iteration.
    info("Solving nonlinear problem:");
    bool verbose = true;
    if (!hermes_2D.solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
        NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

    // Update previous time level solutions.
    Solution::vector_to_solutions(coeff_vec, spaces, slns);
    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", current_time);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time);
    sprintf(title, "Pressure, time %g", current_time);
    pview.set_title(title);
    pview.show(&p_prev_time);
    tview.show(&t_prev_time);

    // Update current time.
    current_time += time_step;
  }

  // Clean up.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

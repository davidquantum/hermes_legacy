#define HERMES_REPORT_ALL
#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_nist-04 benchmarks/nist-04
 *  \{
 *  \brief This test makes sure that the benchmark "nist-04" works correctly.
 *
 *  \section s_params Parameters
 *   - INIT_REF_NUM=2
 *   - P_INIT=1
 *   - THRESHOLD=0.3
 *   - STRATEGY=0
 *   - CAND_LIST=H2D_HP_ANISO;
 *   - MESH_REGULARITY=-1
 *   - CONV_EXP=0.5
 *   - ERR_STOP=0.01
 *   - NDOF_STOP=60000
 *   - matrix_solver = SOLVER_UMFPACK
 *   - ALPHA = 1000
 *   - X_LOC = 0.5
 *   - Y_LOC = 0.5
 *
 *
 *  \section s_res Results
 *   - DOFs: 941
 *   - Adaptivity steps: 11
 */

const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 0.5;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.                      
double alpha = 1000;        // This problem has and exponential peak in the interior of the domain.
double x_loc = 0.5;         // (x_loc, y_loc) is the location of the peak, and alpha determines the strenghth of the peak.
double y_loc = 0.5;

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../square_quad.mesh", &mesh);     // quadrilaterals
  // mloader.load("square_tri.mesh", &mesh);   // triangles

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Set exact solution.
  CustomExactSolution exact(&mesh, alpha, x_loc, y_loc);

  // Define right-hand side.
  CustomFunction f(alpha, x_loc, y_loc);

  // Initialize the weak formulation.
  HermesFunction lambda(1.0);
  WeakFormsH1::DefaultWeakFormPoisson wf(HERMES_ANY, &lambda, &f);
  
  // Initialize boundary conditions
  DefaultEssentialBCNonConst bc("Bdy", &exact);
  EssentialBCs bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = Space::construct_refined_space(&space);
    int ndof_ref = Space::get_num_dofs(ref_space);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem dp(&wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(scalar));

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution ref_sln;
    Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    Solution sln;
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    Adapt* adaptivity = new Adapt(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Calculate exact error.
    double err_exact_rel = hermes2d.calc_rel_error(&sln, &exact, HERMES_H1_NORM) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d", Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space::get_num_dofs(&space), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete [] coeff_vec;
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = Space::get_num_dofs(&space);

  int n_dof_allowed = 950;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}


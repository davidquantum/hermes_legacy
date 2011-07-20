#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_nist-03 Benchmarks/nist-03
 *  \{
 *  \brief This test makes sure that the benchmark "nist-03" works correctly.
 *
 *  \section s_params Parameters
 *   - INIT_REF_NUM = 1
 *   - P_INIT_U = 2
 *   - P_INIT_V = 2
 *   - THRESHOLD = 0.3
 *   - STRATEGY = 0
 *   - CAND_LIST = H2D_HP_ANISO;
 *   - MESH_REGULARITY = -1
 *   - CONV_EXP = 1.0
 *   - ERR_STOP = 0.1
 *   - NDOF_STOP = 60000
 *   - matrix_solver = SOLVER_UMFPACK
 *
 *  \section s_res Results
 *   - DOFs: 1309
 *   - Adaptivity steps: 35
 */

const int P_INIT_U = 2;                           // Initial polynomial degree for u.
const int P_INIT_V = 2;                           // Initial polynomial degree for v.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
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
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1;                        // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double E = 1.0;                             // Young modulus.
const double nu = 0.3;                            // Poisson ratio.
#ifdef MODE_1
  const double lambda = 0.5444837367825;          // lambda for mode-1 solution.
  const double mu = E / (2 * (1 + nu));           // mu for mode-1 solution (mu is the same as G)
  const double Q = 0.5430755788367;               // Q for mode-1 solution.
#else
  const double lambda = 0.9085291898461;          // lambda for mode-2 solution.
  const double mu = E / (2 * (1 + nu));           // mu for mode-2 solution (mu is the same as G).
  const double Q = -0.2189232362488;              // Q for mode-2 solution.
#endif

// Exact solution, weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u_mesh, v_mesh;
  H2DReader mloader;
  mloader.load("../elasticity.mesh", &u_mesh);

  // Create initial mesh (master mesh).
  v_mesh.copy(&u_mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) u_mesh.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM; i++) v_mesh.refine_all_elements();

  // Set exact solution for each displacement component.
  CustomExactSolutionU exact_u(&u_mesh, E, nu, lambda, Q);
  CustomExactSolutionV exact_v(&v_mesh, E, nu, lambda, Q);

  // Initialize boundary conditions
  DefaultEssentialBCNonConst bc_u("Bdy", &exact_u);
  EssentialBCs bcs_u(&bc_u);
  DefaultEssentialBCNonConst bc_v("Bdy", &exact_v);
  EssentialBCs bcs_v(&bc_v);

  // Initialize the weak formulation.
  // NOTE: We pass all four parameters since in Mitchell's 
  // papers they are mutually inconsistent.
  CustomWeakFormElasticityNIST wf(E, nu, mu, lambda);
  
  // Create H1 spaces with default shapeset for both displacement components.
  H1Space u_space(&u_mesh, &bcs_u, P_INIT_U);
  H1Space v_space(&v_mesh, &bcs_v, P_INIT_V);

  Solution u_sln, v_sln, u_ref_sln, v_ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&u_space, &v_space));
    int ndof_ref = Space::get_num_dofs(*ref_spaces);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
    solver->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING);

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem dp(&wf, *ref_spaces);

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(scalar));

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solutions(coeff_vec, *ref_spaces, Hermes::vector<Solution *>(&u_ref_sln, &v_ref_sln));

    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space *>(&u_space, &v_space), 
                                 Hermes::vector<Solution *>(&u_ref_sln, &v_ref_sln), 
                                 Hermes::vector<Solution *>(&u_sln, &v_sln), matrix_solver); 

    // Calculate element errors.
    info("Calculating error estimate and exact error."); 
    Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&u_space, &v_space));
    
    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&u_sln, &v_sln), 
                               Hermes::vector<Solution *>(&u_ref_sln, &v_ref_sln), &err_est_rel) * 100;

    // Calculate exact error for each solution component and the total exact error.
    Hermes::vector<double> err_exact_rel;
    bool solutions_for_adapt = false;
    double err_exact_rel_total = adaptivity->calc_err_exact(Hermes::vector<Solution *>(&u_sln, &v_sln), 
                                 Hermes::vector<Solution *>(&exact_u, &exact_v), 
                                 &err_exact_rel, solutions_for_adapt) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse[u]: %d, ndof_fine[u]: %d",
         u_space.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[0]));

    info("err_est_rel[u]: %g%%, err_exact_rel[u]: %g%%", err_est_rel[0]*100, err_exact_rel[0]*100);

    info("ndof_coarse[v]: %d, ndof_fine[v]: %d",
         v_space.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[1]));

    info("err_est_rel[v]: %g%%, err_exact_rel[v]: %g%%", err_est_rel[1]*100, err_exact_rel[1]*100);

    info("ndof_coarse_total: %d, ndof_fine_total: %d",
         Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)), Space::get_num_dofs(*ref_spaces));

    info("err_est_rel_total: %g%%, err_est_exact_total: %g%%", err_est_rel_total, err_exact_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)), err_exact_rel_total);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector), 
                               THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)) >= NDOF_STOP) done = true;

    // Clean up.
    delete [] coeff_vec;
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
    delete ref_spaces;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space));

  int n_dof_allowed = 135;
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

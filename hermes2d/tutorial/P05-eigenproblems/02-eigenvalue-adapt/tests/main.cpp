#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include <stdio.h>

using namespace RefinementSelectors;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

// This test makes sure that example 51-eigenvalue-adapt works correctly.

const int NUMBER_OF_EIGENVALUES = 6;              // Desired number of eigenvalues.
int P_INIT = 2;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial mesh refinements.
double TARGET_VALUE = 2.0;                        // PySparse parameter: Eigenvalues in the vicinity of 
                                                  // this number will be computed. 
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
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
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

int main(int argc, char* argv[])
{
  if (NUMBER_OF_EIGENVALUES > 6) error("Maximum number of eigenvalues is 6.");
  info("Desired number of eigenvalues: %d.", NUMBER_OF_EIGENVALUES);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions. 
  DefaultEssentialBCConst bc_essential("Bdy", 0.0);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);

  // Initialize the weak formulation.
  WeakFormEigenLeft wf_left;
  WeakFormEigenRight wf_right;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  Solution sln[NUMBER_OF_EIGENVALUES], ref_sln[NUMBER_OF_EIGENVALUES];

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = Space::construct_refined_space(&space);
    int ref_ndof = Space::get_num_dofs(ref_space);
    info("ref_ndof: %d.", ref_ndof);

    // Initialize matrices and matrix solver on referenc emesh.
    RCP<SparseMatrix> matrix_left = rcp(new CSCMatrix());
    RCP<SparseMatrix> matrix_right = rcp(new CSCMatrix());
    Solver* solver = create_linear_solver(matrix_solver, matrix_left.get());

    // Assemble the matrices on reference mesh.
    DiscreteProblem* dp_left = new DiscreteProblem(&wf_left, ref_space);
    dp_left->assemble(matrix_left.get());
    DiscreteProblem* dp_right = new DiscreteProblem(&wf_right, ref_space);
    dp_right->assemble(matrix_right.get());

    // Time measurement.
    cpu_time.tick();

    // Time measurement.
    cpu_time.tick(HERMES_SKIP);

    EigenSolver es(matrix_left, matrix_right);
    info("Calling Pysparse...");
    es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
    info("Pysparse finished.");
    es.print_eigenvalues();

    // Initializing solution vector, solution and ScalarView.
    double* ref_coeff_vec;

    // Reading solution vectors from file and visualizing.
    double eigenval[NUMBER_OF_EIGENVALUES];
    int neig = es.get_n_eigs();
    if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
    for (int ieig = 0; ieig < NUMBER_OF_EIGENVALUES; ieig++) {
      eigenval[ieig] = es.get_eigenvalue(ieig);
      int n;
      es.get_eigenvector(ieig, &ref_coeff_vec, &n);
      // Convert coefficient vector into a Solution.
      Solution::vector_to_solution(ref_coeff_vec, ref_space, &(ref_sln[ieig]));

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution %d on coarse mesh.", ieig);
      OGProjection::project_global(&space, &(ref_sln[ieig]), &(sln[ieig]), matrix_solver);
    }  

     // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    Hermes::vector<Space *> spaces;
    for(int i = 0; i < NUMBER_OF_EIGENVALUES; i++) spaces.push_back(&space);
    Adapt* adaptivity = new Adapt(spaces);
 
    Hermes::vector<Solution *> slns;
    if (NUMBER_OF_EIGENVALUES > 0) slns.push_back(&sln[0]);
    if (NUMBER_OF_EIGENVALUES > 1) slns.push_back(&sln[1]);
    if (NUMBER_OF_EIGENVALUES > 2) slns.push_back(&sln[2]);
    if (NUMBER_OF_EIGENVALUES > 3) slns.push_back(&sln[3]);
    if (NUMBER_OF_EIGENVALUES > 4) slns.push_back(&sln[4]);
    if (NUMBER_OF_EIGENVALUES > 5) slns.push_back(&sln[5]);
    
    Hermes::vector<Solution *> ref_slns;
    if (NUMBER_OF_EIGENVALUES > 0) ref_slns.push_back(&ref_sln[0]);
    if (NUMBER_OF_EIGENVALUES > 1) ref_slns.push_back(&ref_sln[1]);
    if (NUMBER_OF_EIGENVALUES > 2) ref_slns.push_back(&ref_sln[2]);
    if (NUMBER_OF_EIGENVALUES > 3) ref_slns.push_back(&ref_sln[3]);
    if (NUMBER_OF_EIGENVALUES > 4) ref_slns.push_back(&ref_sln[4]);
    if (NUMBER_OF_EIGENVALUES > 5) ref_slns.push_back(&ref_sln[5]);
    Hermes::vector<double> component_errors;
    double err_est_rel = adaptivity->calc_err_est(slns, ref_slns, &component_errors) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d.", Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    if (NUMBER_OF_EIGENVALUES > 0) info("err_est_rel[0]: %g%%", component_errors[0] * 100);
    if (NUMBER_OF_EIGENVALUES > 1) info("err_est_rel[1]: %g%%", component_errors[1] * 100);
    if (NUMBER_OF_EIGENVALUES > 2) info("err_est_rel[2]: %g%%", component_errors[2] * 100);
    if (NUMBER_OF_EIGENVALUES > 3) info("err_est_rel[3]: %g%%", component_errors[3] * 100);
    if (NUMBER_OF_EIGENVALUES > 4) info("err_est_rel[4]: %g%%", component_errors[4] * 100);
    if (NUMBER_OF_EIGENVALUES > 5) info("err_est_rel[5]: %g%%", component_errors[5] * 100);
   
    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      info("Adapting coarse mesh.");
      Hermes::vector<RefinementSelectors::Selector *> selectors;
      for(int i = 0; i < NUMBER_OF_EIGENVALUES; i++)
        selectors.push_back(&selector);
      done = adaptivity->adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    delete dp_left;
    delete dp_right;
  }
  while (done == false);

  int ndof = Space::get_num_dofs(&space);

  if (ndof < 450) {          // Was 401 when this test was created.
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
 
}

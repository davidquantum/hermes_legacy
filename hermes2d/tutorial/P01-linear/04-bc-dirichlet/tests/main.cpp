#define HERMES_REPORT_ALL
#include "../definitions.h"

// This example is a continuation of the previous one and it shows how to 
// use non-constant Dirichlet boundary conditions.
//
// PDE: Poisson equation -div(LAMBDA grad u) = VOLUME_HEAT_SRC 
//
// Boundary conditions: Dirichlet u = bdy_function(x, y).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 5;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 6e3;        // Volume heat sources generated by electric current.        
const double BDY_A_PARAM = 1.0;
const double BDY_B_PARAM = 2.0;
const double BDY_C_PARAM = 20.0;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoissonDirichlet wf("Aluminum", LAMBDA_AL, "Copper", 
                                    LAMBDA_CU, VOLUME_HEAT_SRC);
  
  // Initialize boundary conditions.
  CustomDirichletCondition bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
                                        BDY_A_PARAM, BDY_B_PARAM, BDY_C_PARAM);
  EssentialBCs bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Perform Newton's iteration.
  if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution sln;
  Solution::vector_to_solution(coeff_vec, &space, &sln);

  // VTK output.
  if (VTK_VISUALIZATION) {
    // Output solution in VTK format.
    Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  ndof = Space::get_num_dofs(&space);
  printf("ndof = %d\n", ndof);
  double sum = 0;
  for (int i=0; i < ndof; i++) sum += coeff_vec[i];
  printf("coefficient sum = %g\n", sum);

  bool success = true;
  if (fabs(sum + 4.2471) > 1e-3) success = 0;

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }

}

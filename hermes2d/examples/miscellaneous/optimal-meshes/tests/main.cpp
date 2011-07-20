#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace RefinementSelectors;

// This test makes sure that example "optimal-meshes" works correctly.

const int P_INIT = 1;                             // Initial polynomial degree.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double K = 10.0;

int main(int argc, char* argv[])
{
  // Initialize the library's global functions.
  Hermes2D hermes2D;

  // Read the command-line arguments.
  if (argc != 10) error("You must provide 5 real numbers (mesh vertices) and 4 integers (poly degrees).");
  double x0 = atof(argv[1]);
  double x1 = atof(argv[2]);
  double x2 = atof(argv[3]);
  double x3 = atof(argv[4]);
  double x4 = atof(argv[5]);
  int o0 = atoi(argv[6]);
  int o1 = atoi(argv[7]);
  int o2 = atoi(argv[8]);
  int o3 = atoi(argv[9]);

  // Prepare mesh geometry.
  int nv = 10;
  double2 verts[10];
  verts[0][0] = x0; verts[0][1] = 0;
  verts[1][0] = x1; verts[1][1] = 0;
  verts[2][0] = x2; verts[2][1] = 0;
  verts[3][0] = x3; verts[3][1] = 0;
  verts[4][0] = x4; verts[4][1] = 0;
  verts[5][0] = x0; verts[5][1] = 1;
  verts[6][0] = x1; verts[6][1] = 1;
  verts[7][0] = x2; verts[7][1] = 1;
  verts[8][0] = x3; verts[8][1] = 1;
  verts[9][0] = x4; verts[9][1] = 1;
  int nt = 0;
  int4* tris = NULL;
  int nq = 4;
  int5 quads[4];
  quads[0][0] = 0; quads[0][1] = 1; quads[0][2] = 6; quads[0][3] = 5; quads[0][4] = 0;
  quads[1][0] = 1; quads[1][1] = 2; quads[1][2] = 7; quads[1][3] = 6; quads[1][4] = 0;
  quads[2][0] = 2; quads[2][1] = 3; quads[2][2] = 8; quads[2][3] = 7; quads[2][4] = 0;
  quads[3][0] = 3; quads[3][1] = 4; quads[3][2] = 9; quads[3][3] = 8; quads[3][4] = 0;
  int nm = 10;
  int3 mark[10];
  mark[0][0] = 0; mark[0][1] = 1; mark[0][2] = 1;
  mark[1][0] = 1; mark[1][1] = 2; mark[1][2] = 1;
  mark[2][0] = 2; mark[2][1] = 3; mark[2][2] = 1;
  mark[3][0] = 3; mark[3][1] = 4; mark[3][2] = 1;
  mark[4][0] = 4; mark[4][1] = 9; mark[4][2] = 1;
  mark[5][0] = 9; mark[5][1] = 8; mark[5][2] = 1;
  mark[6][0] = 8; mark[6][1] = 7; mark[6][2] = 1;
  mark[7][0] = 7; mark[7][1] = 6; mark[7][2] = 1;
  mark[8][0] = 6; mark[8][1] = 5; mark[8][2] = 1;
  mark[9][0] = 5; mark[9][1] = 0; mark[9][2] = 1;

  // Create a mesh with 10 vertices, 4 elements and 10 boundary 
  // edges from the above data.
  Mesh mesh;
  mesh.create(nv, verts, nt, tris, nq, quads, nm, mark);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, P_INIT);

  // Set element poly orders.
  space.set_element_order(0, H2D_MAKE_QUAD_ORDER(o0, 1));
  space.set_element_order(1, H2D_MAKE_QUAD_ORDER(o1, 1));
  space.set_element_order(2, H2D_MAKE_QUAD_ORDER(o2, 1));
  space.set_element_order(3, H2D_MAKE_QUAD_ORDER(o3, 1));

  // Perform orthogonal projection in the H1 norm.
  Solution sln_approx;
  CustomExactSolution sln_exact(&mesh, K);
  OGProjection::project_global(&space, &sln_exact, &sln_approx);

  // Calculate the error.
  double err = hermes2D.calc_abs_error(&sln_approx, &sln_exact, HERMES_H1_NORM);
  printf("\nMesh: %g, %g, %g, %g, %g\n", x0, x1, x2, x3, x4);
  printf("Poly degrees: %d, %d, %d, %d\n", o0, o1, o2, o3);
  printf("err = %g, err_squared = %g\n\n", err, err*err);

  // Mesh: 0, 1, 2, 3, 4
  // Poly degrees: 10, 10, 10, 10
  if ((err - 0.04381394) < 1E-6) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}


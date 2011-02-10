#define HERMES_REPORT_INFO

// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "hermes2d.h"

// TODO LIST: 
//
// (1) Incorporate spatial adaptivity into the time stepping.
//
// (2) Enable more equations than one. Right now rk_time_step() does not 
//     work for systems.
//
// (3) Enable all other matrix solvers, so far UMFPack is hardwired here.
//
// (4) We do not take advantage of the fact that all blocks in the 
//     Jacobian matrix have the same structure. Thus it is enough to 
//     assemble the matrix M (one block) and copy the sparsity structure
//     into all remaining nonzero blocks (and diagonal blocks). Right 
//     now, the sparsity structure is created expensively in each block 
//     again.
//
// (5) If space does not change, the sparsity does not change. Right now 
//     we discard everything at the end of every time step, we should not 
//     do it.  
//
// (6) If the problem does not depend explicitly on time, then all the blocks 
//     in the Jacobian matrix of the stationary residual are the same up 
//     to a multiplicative constant. Thus they do not have to be aassembled 
//     from scratch.
// 
// (7) If the problem is linear, then the Jacobian is constant. If Space 
//     does not change between time steps, we should keep it. 

void create_stage_wf(double current_time, double time_step, ButcherTable* bt, 
                     DiscreteProblem* dp, WeakForm* stage_wf_left, 
                     WeakForm* stage_wf_right) 
{
  // First let's do the mass matrix (only one block ndof times ndof).
  WeakForm::MatrixFormVol mfv_00;
  mfv_00.i = 0;
  mfv_00.j = 0;
  mfv_00.sym = HERMES_SYM;
  mfv_00.area = HERMES_ANY;
  mfv_00.fn = l2_form<double, scalar>;
  mfv_00.ord = l2_form<Ord, Ord>;
  mfv_00.ext = Hermes::vector<MeshFunction*> ();
  mfv_00.scaling_factor = 1.0;
  mfv_00.u_ext_offset = 0;
  stage_wf_left->add_matrix_form(&mfv_00);

  // In the rest we will take the stationary jacobian and residual forms 
  // (right-hand side) and use them to create a block Jacobian matrix of
  // size (num_stages*ndof times num_stages*ndof) and a block residual 
  // vector of length num_stages*ndof.

  // Number of stages.
  int num_stages = bt->get_size();

  // Original weak formulation.
  WeakForm* wf = dp->get_weak_formulation();

  // Extract mesh from (the first space of) the original discrete problem.
  Mesh* mesh = dp->get_space(0)->get_mesh();

  // Create constant Solutions to represent the stage times,
  // stage_time = current_time + c_i*time_step.
  // WARNING - THIS IS A TEMPORARY HACK. THE STAGE TIME SHOULD BE ENTERED
  // AS A NUMBER, NOT IN THIS WAY. IT WILL BE ADDED AFTER EXISTING EXTERNAL
  // SOLUTIONS IN ExtData.
  Solution** stage_time_sol = new Solution*[num_stages];
  for (int i = 0; i < num_stages; i++) {
    stage_time_sol[i] = new Solution(mesh);
    stage_time_sol[i]->set_const(mesh, current_time + bt->get_C(i)*time_step);
  }

  // Extracting volume and surface matrix and vector forms from the
  // original weak formulation.
  if (wf->get_neq() != 1) error("wf->neq != 1 not implemented yet.");
  Hermes::vector<WeakForm::MatrixFormVol> mfvol_base = wf->get_mfvol();
  Hermes::vector<WeakForm::MatrixFormSurf> mfsurf_base = wf->get_mfsurf();
  Hermes::vector<WeakForm::VectorFormVol> vfvol_base = wf->get_vfvol();
  Hermes::vector<WeakForm::VectorFormSurf> vfsurf_base = wf->get_vfsurf();

  // Duplicate matrix volume forms, scale them according
  // to the Butcher's table, enhance them with additional
  // external solutions, and anter them as blocks to the
  // new stage Jacobian.
  for (unsigned int m = 0; m < mfvol_base.size(); m++) {
    WeakForm::MatrixFormVol mfv_base = mfvol_base[m];
    for (int i = 0; i < num_stages; i++) {
      for (int j = 0; j < num_stages; j++) {
        WeakForm::MatrixFormVol mfv_ij;
        mfv_ij.i = i;
        mfv_ij.j = j;
        mfv_ij.sym = mfv_base.sym;
        mfv_ij.area = mfv_base.area;
        mfv_ij.fn = mfv_base.fn;
        mfv_ij.ord = mfv_base.ord;
        mfv_ij.scaling_factor = -time_step * bt->get_A(i, j);

        // Duplicate external solutions.
        // THIS WAS WRONG, ONE CANNOT COPY INTO AN EMPTY STD::VECTOR:
        // std::copy(mfv_base.ext.begin(), mfv_base.ext.end(), mfv_ij.ext.begin()); 
        for (unsigned int f_idx = 0; f_idx < mfv_base.ext.size(); f_idx++) 
          mfv_ij.ext.push_back(mfv_base.ext[f_idx]);

        // Add stage_time_sol[i] as an external function to the form.
        mfv_ij.ext.push_back(stage_time_sol[i]);

        // Set offset for u_ext[] external solutions.
        mfv_ij.u_ext_offset = i;

        // Add the matrix form to the corresponding block of the
        // stage Jacobian matrix.
        stage_wf_right->add_matrix_form(&mfv_ij);
      }
    }
  }

  // Duplicate matrix surface forms, enhance them with
  // additional external solutions, and anter them as
  // blocks of the stage Jacobian.
  for (unsigned int m = 0; m < mfsurf_base.size(); m++) {
    WeakForm::MatrixFormSurf mfs_base = mfsurf_base[m];
    for (int i = 0; i < num_stages; i++) {
      for (int j = 0; j < num_stages; j++) {
        WeakForm::MatrixFormSurf mfs_ij;
        mfs_ij.i = i;
        mfs_ij.j = j;
        mfs_ij.area = mfs_base.area;
        mfs_ij.fn = mfs_base.fn;
        mfs_ij.ord = mfs_base.ord;
        mfs_ij.scaling_factor = -time_step * bt->get_A(i, j);

        // Duplicate external solutions.
        // THIS WAS WRONG, ONE CANNOT COPY INTO AN EMPTY STD::VECTOR:
        // std::copy(mfs_base.ext.begin(), mfs_base.ext.end(), mfs_ij.ext.begin()); 
        for (unsigned int f_idx = 0; f_idx < mfs_base.ext.size(); f_idx++) 
          mfs_ij.ext.push_back(mfs_base.ext[f_idx]);

        // Add stage_time_sol[i] as an external function to the form.
        mfs_ij.ext.push_back(stage_time_sol[i]);

        // Set offset for u_ext[] external solutions.
        mfs_ij.u_ext_offset = i;

        // Add the matrix form to the corresponding block of the
        // stage Jacobian matrix.
        stage_wf_right->add_matrix_form_surf(&mfs_ij);
      }
    }
  }

  // Duplicate vector volume forms, enhance them with
  // additional external solutions, and anter them as
  // blocks of the stage residual.
  for (unsigned int m = 0; m < vfvol_base.size(); m++) {
    WeakForm::VectorFormVol vfv_base = vfvol_base[m];
    for (int i = 0; i < num_stages; i++) {
      WeakForm::VectorFormVol vfv_i;
      vfv_i.i = i;
      vfv_i.area = vfv_base.area;
      vfv_i.fn = vfv_base.fn;
      vfv_i.ord = vfv_base.ord;
      vfv_i.scaling_factor = -1.0;

      // Duplicate external solutions.
      // THIS WAS WRONG, ONE CANNOT COPY INTO AN EMPTY STD::VECTOR:
      // std::copy(vfv_base.ext.begin(), vfv_base.ext.end(), vfv_ij.ext.begin()); 
      for (unsigned int f_idx = 0; f_idx < vfv_base.ext.size(); f_idx++) 
        vfv_i.ext.push_back(vfv_base.ext[f_idx]);

      // Add stage_time_sol[i] as an external function to the form.
      vfv_i.ext.push_back(stage_time_sol[i]);

      // Set offset for u_ext[] external solutions.
      vfv_i.u_ext_offset = i;

      // Add the matrix form to the corresponding block of the
      // stage Jacobian matrix.
      stage_wf_right->add_vector_form(&vfv_i);
    }
  }

  // Duplicate vector surface forms, enhance them with
  // additional external solutions, and anter them as
  // blocks of the stage residual.
  for (unsigned int m = 0; m < vfsurf_base.size(); m++) {
    WeakForm::VectorFormSurf vfs_base = vfsurf_base[m];
    for (int i = 0; i < num_stages; i++) {
      WeakForm::VectorFormSurf vfs_i;
      vfs_i.i = i;
      vfs_i.area = vfs_base.area;
      vfs_i.fn = vfs_base.fn;
      vfs_i.ord = vfs_base.ord;
      vfs_i.scaling_factor = -1.0;

      // Duplicate external solutions.
      // THIS WAS WRONG, ONE CANNOT COPY INTO AN EMPTY STD::VECTOR:
      // std::copy(vfv_base.ext.begin(), vfv_base.ext.end(), vfv_ij.ext.begin()); 
      for (unsigned int f_idx = 0; f_idx < vfs_base.ext.size(); f_idx++) 
        vfs_i.ext.push_back(vfs_base.ext[f_idx]);
  
      // Add stage_time_sol[i] as an external function to the form.
      vfs_i.ext.push_back(stage_time_sol[i]);

      // Set offset for u_ext[] external solutions.
      vfs_i.u_ext_offset = i;

      // Add the matrix form to the corresponding block of the
      // stage Jacobian matrix.
      stage_wf_right->add_vector_form_surf(&vfs_i);
    }
  }
}

// This takes a matrix, and uses it to formally construct a block-diagonal
// matrix. There are num_blocks blocks on the diagonal. The block diagonal
// matrix is then multiplied with the vector source_vec.
void multiply_as_diagonal_block_matrix(UMFPackMatrix* matrix, int num_blocks,
                                       scalar* source_vec, scalar* target_vec)
{
  int size = matrix->get_size();
  for (int i = 0; i < num_blocks; i++) {
    matrix->multiply_with_vector(source_vec + i*size, target_vec + i*size);
  }
}

bool HERMES_RESIDUAL_AS_VECTOR_RK = true;
bool rk_time_step(double current_time, double time_step, ButcherTable* const bt,
                  Solution* sln, Space* sln_space, Solution* error_fn, DiscreteProblem* dp, 
                  MatrixSolverType matrix_solver,
                  bool verbose, bool is_linear, double newton_tol, int newton_max_iter,
                  double newton_damping_coeff, double newton_max_allowed_residual_norm)
{
  // Check for not implemented features.
  if (matrix_solver != SOLVER_UMFPACK)
    error("Sorry, rk_time_step() still only works with UMFpack.");
  if (dp->get_weak_formulation()->get_neq() > 1)
    error("Sorry, rk_time_step() does not work with systems yet.");

  // Get number of stages from the Butcher's table.
  int num_stages = bt->get_size();

  // Check whether the user provided a nonzero B2-row if he wants temporal error estimation.
  if(error_fn != NULL) {
    double b2_coeff_sum = 0;
    for (int i=0; i < num_stages; i++) b2_coeff_sum += fabs(bt->get_B2(i)); 
    if (b2_coeff_sum < 1e-10) 
      error("error_fn != NULL but the B2 row in the Butcher's table is zero in rk_time_step().");
  }

  // Matrix for the time derivative part of the equation (left-hand side).
  UMFPackMatrix* matrix_left = new UMFPackMatrix();

  // Matrix and vector for the rest (right-hand side).
  UMFPackMatrix* matrix_right = new UMFPackMatrix();
  UMFPackVector* vector_right = new UMFPackVector();

  // Create matrix solver.
  Solver* solver = create_linear_solver(matrix_solver, matrix_right, vector_right);

  // Get space, mesh, and ndof for the stage solutions in the R-K method (K_i vectors).
  Space* K_space = dp->get_space(0);
  Mesh* K_mesh = K_space->get_mesh();
  int ndof = K_space->get_num_dofs();

  // Create spaces for stage solutions K_i. This is necessary
  // to define a num_stages x num_stages block weak formulation.
  Hermes::vector<Space*> stage_spaces;
  stage_spaces.push_back(K_space);
  for (int i = 1; i < num_stages; i++) {
    stage_spaces.push_back(K_space->dup(K_mesh));
  }
  Space::assign_dofs(stage_spaces);

  // Create a multistage weak formulation.
  WeakForm stage_wf_left;                   // For the matrix M (size ndof times ndof).
  WeakForm stage_wf_right(num_stages);      // For the rest of equation (written on the right),
                                            // size num_stages*ndof times num_stages*ndof.
  create_stage_wf(current_time, time_step, bt, dp, &stage_wf_left, &stage_wf_right); 

  // Initialize discrete problems for the assembling of the
  // matrix M and the stage Jacobian matrix and residual.
  DiscreteProblem stage_dp_left(&stage_wf_left, K_space);
  DiscreteProblem stage_dp_right(&stage_wf_right, stage_spaces);

  // Vector K_vector of length num_stages * ndof. will represent
  // the 'K_i' vectors in the usual R-K notation.
  scalar* K_vector = new scalar[num_stages*ndof];
  memset(K_vector, 0, num_stages * ndof * sizeof(scalar));

  // Vector u_ext_vec will represent h \sum_{j=1}^s a_{ij} K_i.
  scalar* u_ext_vec = new scalar[num_stages*ndof];

  // Vector for the left part of the residual.
  scalar* vector_left = new scalar[num_stages*ndof];

  // Prepare residuals of stage solutions.
  Hermes::vector<Solution*> residuals;
  Hermes::vector<bool> add_dir_lift;
  for (int i = 0; i < num_stages; i++) {
    residuals.push_back(new Solution(K_mesh));
    add_dir_lift.push_back(false);
  }

  // Assemble the block-diagonal mass matrix M of size ndof times ndof.
  // The corresponding part of the global residual vector is obtained 
  // just by multiplication.
  stage_dp_left.assemble(matrix_left);

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (true)
  {
    // Prepare vector h\sum_{j=1}^s a_{ij} K_j.
    for (int i = 0; i < num_stages; i++) {                // block row
      for (int idx = 0; idx < ndof; idx++) {
        scalar increment = 0;
        for (int j = 0; j < num_stages; j++) {
          increment += bt->get_A(i, j) * K_vector[j*ndof + idx];
        }
        u_ext_vec[i*ndof + idx] = time_step * increment;
      }
    }

    multiply_as_diagonal_block_matrix(matrix_left, num_stages, K_vector, vector_left);

    // Assemble the block Jacobian matrix of the stationary residual F
    // Diagonal blocks are created even if empty, so that matrix_left
    // can be added later.
    bool rhs_only = false;
    bool force_diagonal_blocks = true;
    stage_dp_right.assemble(u_ext_vec, matrix_right, vector_right,
                            rhs_only, force_diagonal_blocks, false); // false = do not add Dirichlet lift while
                                                                     // converting u_ext_vec into Solutions.

    matrix_right->add_to_diagonal_blocks(num_stages, matrix_left);

    vector_right->add_vector(vector_left);

    // Multiply the residual vector with -1 since the matrix
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    vector_right->change_sign();

    // Measure the residual norm.
    if (HERMES_RESIDUAL_AS_VECTOR_RK) {
      // Calculate the l2-norm of residual vector.
      residual_norm = get_l2_norm(vector_right);
    }
    else {
      // Translate residual vector into residual functions.
      Solution::vector_to_solutions(vector_right, stage_dp_right.get_spaces(),
                                    residuals, add_dir_lift);
      residual_norm = calc_norms(residuals);
    }

    // Info for the user.
    if (verbose) info("---- Newton iter %d, ndof %d, residual norm %g",
                      it, ndof, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > newton_max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", newton_max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number
    // of iteration has been reached, or the problem is linear, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaK^{n+1} to K^n.
    for (int i = 0; i < num_stages*ndof; i++) {
      K_vector[i] += newton_damping_coeff * solver->get_solution()[i];
    }

    // If the problem is linear, quit.
    if (is_linear) {
      if (verbose) {
        info("Terminating Newton's loop as problem is linear.");
      }
      break;
    }

    // Increase iteration counter.
    it++;
  }

  // If max number of iterations was exceeded, fail.
  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }

  // Project previous time level solution on the stage space,
  // to be able to add them together. The result of the projection 
  // will be stored in the vector coeff_vec.
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(K_space, sln, coeff_vec, matrix_solver);

  // Calculate new time level solution in the stage space (u_{n+1} = u_n + h \sum_{j=1}^s b_j k_j).
  for (int i = 0; i < ndof; i++) {
    for (int j = 0; j < num_stages; j++) {
      coeff_vec[i] += time_step * bt->get_B(j) * K_vector[j*ndof + i];
    }
  }
  Solution::vector_to_solution(coeff_vec, K_space, sln);

  // If error_fn is not NULL, use the B2-row in the Butcher's
  // table to calculate the temporal error estimate.
  if (error_fn != NULL) {
    for (int i = 0; i < ndof; i++) {
      coeff_vec[i] = 0;
      for (int j = 0; j < num_stages; j++) {
        coeff_vec[i] += (bt->get_B(j) - bt->get_B2(j)) * K_vector[j*ndof + i];
      }
      coeff_vec[i] *= time_step;
    }
    Solution::vector_to_solution(coeff_vec, K_space, error_fn, false);
  }

  // Clean up.
  delete matrix_left;
  delete matrix_right;
  delete vector_right;
  delete solver;

  // Delete stage spaces, but not the first (original) one.
  for (int i = 1; i < num_stages; i++) delete stage_spaces[i];

  // Delete all residuals.
  for (int i = 0; i < num_stages; i++) delete residuals[i];

  // TODO: Delete stage_wf, in particular its external solutions
  // stage_time_sol[i], i = 0, 1, ..., num_stages-1.

  // Clean up.
  delete [] K_vector;
  delete [] u_ext_vec;
  delete [] coeff_vec;
  delete [] vector_left;

  return true;
}

// This is the same as the rk_time_step() function above but it does not have the error_fn parameter.
bool rk_time_step(double current_time, double time_step, ButcherTable* const bt,
                  Solution* sln, Space* sln_space, DiscreteProblem* dp, MatrixSolverType matrix_solver,
                  bool verbose, bool is_linear, double newton_tol, int newton_max_iter,
                  double newton_damping_coeff, double newton_max_allowed_residual_norm) 
{
  return rk_time_step(current_time, time_step, bt,
	              sln, sln_space, NULL, dp, matrix_solver,
	              verbose, is_linear, newton_tol, newton_max_iter,
                      newton_damping_coeff, newton_max_allowed_residual_norm);
}

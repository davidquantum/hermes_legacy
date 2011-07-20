#include "hermes2d.h"

void OGProjection::project_internal(Hermes::vector<Space *> spaces, WeakForm* wf,
                                    scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  unsigned int n = spaces.size();

  // sanity checks
  if (n <= 0 || n > 10) error("Wrong number of projected functions in project_internal().");
  for (unsigned int i = 0; i < n; i++) if(spaces[i] == NULL) error("this->spaces[%d] == NULL in project_internal().", i);
  if (spaces.size() != n) error("Number of spaces must match number of projected functions in project_internal().");

  // this is needed since spaces may have their DOFs enumerated only locally.
  int ndof = Space::assign_dofs(spaces);

  // Initialize DiscreteProblem.
  DiscreteProblem* dp = new DiscreteProblem(wf, spaces);

  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Perform Newton's iteration.
  if (!hermes2d.solve_newton(coeff_vec, dp, solver, matrix, rhs)) error("Newton's iteration failed.");

  if (target_vec != NULL)
    for (int i=0; i < ndof; i++) target_vec[i] = coeff_vec[i];

  delete [] coeff_vec;
  delete solver;
  delete matrix;
  delete rhs;
  delete dp;
  //delete wf;
}

void OGProjection::project_global(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction*> source_meshfns,
                   scalar* target_vec, MatrixSolverType matrix_solver, Hermes::vector<ProjNormType> proj_norms)
{
  _F_
  int n = spaces.size();

  // define temporary projection weak form
  WeakForm* proj_wf = new WeakForm(n);
  int found[100];
  for (int i = 0; i < 100; i++) found[i] = 0;
  for (int i = 0; i < n; i++)
  {
    ProjNormType norm = HERMES_UNSET_NORM;
    if (proj_norms == Hermes::vector<ProjNormType>()) {
      ESpaceType space_type = spaces[i]->get_type();
      switch (space_type) {
        case HERMES_H1_SPACE: norm = HERMES_H1_NORM; break;
        case HERMES_HCURL_SPACE: norm = HERMES_HCURL_NORM; break;
        case HERMES_HDIV_SPACE: norm = HERMES_HDIV_NORM; break;
        case HERMES_L2_SPACE: norm = HERMES_L2_NORM; break;
        default: error("Unknown space type in OGProjection::project_global().");
      }
    }
    else norm = proj_norms[i];

    // FIXME - memory leak - create Projection class and encapsulate this function project_global(...)
    // maybe in more general form
    found[i] = 1;
    // Jacobian.
    proj_wf->add_matrix_form(new ProjectionMatrixFormVol(i, i, norm));
    // Residual.
    proj_wf->add_vector_form(new ProjectionVectorFormVol(i, source_meshfns[i], norm));
  }
  for (int i=0; i < n; i++)
  {
    if (found[i] == 0)
    {
      warn("Index of component: %d\n", i);
      error("Wrong projection norm in project_global().");
    }
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
  
  delete proj_wf;
}

void OGProjection::project_global(Hermes::vector<Space *> spaces, Hermes::vector<Solution*> source_sols,
                   scalar* target_vec, MatrixSolverType matrix_solver, Hermes::vector<ProjNormType> proj_norms)
{
  Hermes::vector<MeshFunction *> mesh_fns;
  for(unsigned int i = 0; i < source_sols.size(); i++)
    mesh_fns.push_back(source_sols[i]);
  project_global(spaces, mesh_fns, target_vec, matrix_solver, proj_norms);
}

void OGProjection::project_global(Space* space, MeshFunction* source_meshfn,
                             scalar* target_vec, MatrixSolverType matrix_solver,
                             ProjNormType proj_norm)
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<MeshFunction *> source_meshfns;
  source_meshfns.push_back(source_meshfn);
  Hermes::vector<ProjNormType> proj_norms;
  proj_norms.push_back(proj_norm);
  project_global(spaces, source_meshfns, target_vec, matrix_solver, proj_norms);
}


void OGProjection::project_global(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> sols_src,
                                  Hermes::vector<Solution *> sols_dest, MatrixSolverType matrix_solver,
                                  Hermes::vector<ProjNormType> proj_norms, bool delete_old_meshes)
{
  _F_

  scalar* target_vec = new scalar[Space::get_num_dofs(spaces)];
  Hermes::vector<MeshFunction *> ref_slns_mf;
  for (unsigned int i = 0; i < sols_src.size(); i++)
    ref_slns_mf.push_back(static_cast<MeshFunction*>(sols_src[i]));

  OGProjection::project_global(spaces, ref_slns_mf, target_vec, matrix_solver, proj_norms);

  if(delete_old_meshes)
    for(unsigned int i = 0; i < sols_src.size(); i++) {
      delete sols_src[i]->get_mesh();
      sols_src[i]->own_mesh = false;
    }

  Solution::vector_to_solutions(target_vec, spaces, sols_dest);

  delete [] target_vec;
}

void OGProjection::project_global(Space * space,
                             Solution* sol_src, Solution* sol_dest,
                             MatrixSolverType matrix_solver,
                             ProjNormType proj_norm)
{
  Hermes::vector<Space *> spaces;
  spaces.push_back(space);
  Hermes::vector<Solution *> sols_src;
  sols_src.push_back(sol_src);
  Hermes::vector<Solution *> sols_dest;
  sols_dest.push_back(sol_dest);
  Hermes::vector<ProjNormType> proj_norms;
  if(proj_norm != HERMES_UNSET_NORM)
    proj_norms.push_back(proj_norm);
  
  project_global(spaces, sols_src, sols_dest, matrix_solver, proj_norms);
}

void OGProjection::project_global(Hermes::vector<Space *> spaces,
                                  Hermes::vector<WeakForm::MatrixFormVol *> custom_projection_jacobian,
                                  Hermes::vector<WeakForm::VectorFormVol *> custom_projection_residual,
                                  scalar* target_vec, MatrixSolverType matrix_solver)
{
  _F_
  unsigned int n = spaces.size();
  unsigned int n_biforms = custom_projection_jacobian.size();
  if (n_biforms == 0)
    error("Please use the simpler version of project_global with the argument Hermes::vector<ProjNormType> proj_norms if you do not provide your own projection norm.");
  if (n_biforms != custom_projection_residual.size())
    error("Mismatched numbers of projection forms in project_global().");
  if (n != n_biforms)
    error("Mismatched numbers of projected functions and projection forms in project_global().");

  // Define local projection weak form.
  WeakForm* proj_wf = new WeakForm(n);
  for (unsigned int i = 0; i < n; i++) {
    proj_wf->add_matrix_form(custom_projection_jacobian[i]);
    proj_wf->add_vector_form(custom_projection_residual[i]);
  }

  project_internal(spaces, proj_wf, target_vec, matrix_solver);
  
  delete proj_wf;
}

void OGProjection::project_global(Hermes::vector<Space *> spaces,
                                  Hermes::vector<WeakForm::MatrixFormVol *> custom_projection_jacobian,
                                  Hermes::vector<WeakForm::VectorFormVol *> custom_projection_residual,
                                  Hermes::vector<Solution *> sols_dest, MatrixSolverType matrix_solver)
{
  _F_
  scalar* target_vec = new scalar[Space::get_num_dofs(spaces)];
  OGProjection::project_global(spaces, custom_projection_jacobian, custom_projection_residual, target_vec, matrix_solver);
  Solution::vector_to_solutions(target_vec, spaces, sols_dest);
  delete [] target_vec;
}

void OGProjection::project_global(Space* space,
                                  WeakForm::MatrixFormVol* custom_projection_jacobian,
                                  WeakForm::VectorFormVol* custom_projection_residual,
                                  Solution* sol_dest, MatrixSolverType matrix_solver)
{
  _F_
  project_global(Hermes::vector<Space*>(space),
                 Hermes::vector<WeakForm::MatrixFormVol*>(custom_projection_jacobian),
                 Hermes::vector<WeakForm::VectorFormVol*>(custom_projection_residual),
                 Hermes::vector<Solution*>(sol_dest), matrix_solver);
}
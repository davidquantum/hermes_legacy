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


#ifndef __H2D_OGPROJECTION_H
#define __H2D_OGPROJECTION_H

#include "../hermes_common/common.h"
#include "function/solution.h"
#include "function/forms.h"

class HERMES_API OGProjection
{
public:
  static void project_global(Hermes::vector<Space *> spaces, Hermes::vector<MeshFunction *> source_meshfns,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK,
                             Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

  static void project_global(Hermes::vector<Space *> spaces, Hermes::vector<Solution *> source_sols,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK,
                             Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());

  static void project_global(Space* space, MeshFunction* source_meshfn,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK,
                             ProjNormType proj_norm = HERMES_H1_NORM);

  static void project_global(Hermes::vector<Space *> spaces,
                             Hermes::vector<Solution*> sols_src, Hermes::vector<Solution*> sols_dest,
                             MatrixSolverType matrix_solver = SOLVER_UMFPACK,
                             Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>(), bool delete_old_mesh = false);

  static void project_global(Space * space,
                             Solution* sol_src, Solution* sol_dest,
                             MatrixSolverType matrix_solver = SOLVER_UMFPACK,
                             ProjNormType proj_norm = HERMES_UNSET_NORM);

  static void project_global(Hermes::vector<Space *> spaces,
                             Hermes::vector<WeakForm::MatrixFormVol *> custom_projection_jacobian,
                             Hermes::vector<WeakForm::VectorFormVol *> custom_projection_residual,
                             scalar* target_vec, MatrixSolverType matrix_solver = SOLVER_UMFPACK);
  
  static void project_global(Hermes::vector<Space *> spaces,
                             Hermes::vector<WeakForm::MatrixFormVol *> custom_projection_jacobian,
                             Hermes::vector<WeakForm::VectorFormVol *> custom_projection_residual,
                             Hermes::vector<Solution *> sols_dest, MatrixSolverType matrix_solver = SOLVER_UMFPACK);
                             
  static void project_global(Space* space,
                             WeakForm::MatrixFormVol* custom_projection_jacobian,
                             WeakForm::VectorFormVol* custom_projection_residual,
                             Solution* sol_dest, MatrixSolverType matrix_solver = SOLVER_UMFPACK);      

  // Underlying function for global orthogonal projection.
  // Not intended for the user. NOTE: the weak form here must be
  // a special projection weak form, which is different from
  // the weak form of the PDE. If you supply a weak form of the
  // PDE, the PDE will just be solved.
protected:
  static void project_internal(Hermes::vector<Space *> spaces, WeakForm *proj_wf, scalar* target_vec,
                               MatrixSolverType matrix_solver = SOLVER_UMFPACK);

  // Jacobian matrix (same as stiffness matrix since projections are linear).
  class ProjectionMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    ProjectionMatrixFormVol(int i, int j, ProjNormType projNormType) : WeakForm::MatrixFormVol(i, j)
    {
      this->adapt_eval = false;
      this->projNormType = projNormType;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                 Geom<double> *e, ExtData<scalar> *ext) const
    {
      switch (projNormType)
      {
      case HERMES_L2_NORM:
        return l2_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_NORM:
        return h1_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_SEMINORM:
        return h1_semi_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HCURL_NORM:
        return hcurl_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HDIV_NORM:
        return hdiv_projection_biform<double, scalar>(n, wt, u_ext, u, v, e, ext);
      default:
        error("Unknown projection type");
        return 0.0;
      }
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      switch (projNormType)
      {
      case HERMES_L2_NORM:
        return l2_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_NORM:
        return h1_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_H1_SEMINORM:
        return h1_semi_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HCURL_NORM:
        return hcurl_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      case HERMES_HDIV_NORM:
        return hdiv_projection_biform<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
      default:
        error("Unknown projection type");
        return Ord();
      }
    }

  private:
    ProjNormType projNormType;

    template<typename Real, typename Scalar>
    static Scalar h1_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i] + u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar h1_semi_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                            Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar l2_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (u->val[i] * v->val[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar hcurl_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                          Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * (u->curl[i] * conj(v->curl[i]));
        result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      }
      return result;
    }

    template<typename Real, typename Scalar>
    static Scalar hdiv_projection_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * (u->div[i] * conj(v->div[i]));
        result += wt[i] * (u->val0[i] * conj(v->val0[i]) + u->val1[i] * conj(v->val1[i]));
      }
      return result;
    }
  };

  // Residual.
  class ProjectionVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    ProjectionVectorFormVol(int i, MeshFunction* ext, ProjNormType projNormType) : WeakForm::VectorFormVol(i)
    {
      this->adapt_eval = false;
      this->projNormType = projNormType;
      this->ext = Hermes::vector<MeshFunction *>();
      this->ext.push_back(ext);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                 Geom<double> *e, ExtData<scalar> *ext) const
    {
      switch (projNormType)
      {
      case HERMES_L2_NORM:
        return l2_projection_residual<double, scalar>(n, wt, u_ext, v, e, ext);
      case HERMES_H1_NORM:
        return h1_projection_residual<double, scalar>(n, wt, u_ext, v, e, ext);
      case HERMES_H1_SEMINORM:
        return h1_semi_projection_residual<double, scalar>(n, wt, u_ext, v, e, ext);
      case HERMES_HCURL_NORM:
        return hcurl_projection_residual<double, scalar>(n, wt, u_ext, v, e, ext);
      case HERMES_HDIV_NORM:
        return hdiv_projection_residual<double, scalar>(n, wt, u_ext, v, e, ext);
      default:
        error("Unknown projection type");
        return 0.0;
      }
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
            Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      switch (projNormType)
      {
      case HERMES_L2_NORM:
        return l2_projection_residual<Ord, Ord>(n, wt, u_ext, v, e, ext);
      case HERMES_H1_NORM:
        return h1_projection_residual<Ord, Ord>(n, wt, u_ext, v, e, ext);
      case HERMES_H1_SEMINORM:
        return h1_semi_projection_residual<Ord, Ord>(n, wt, u_ext, v, e, ext);
      case HERMES_HCURL_NORM:
        return hcurl_projection_residual<Ord, Ord>(n, wt, u_ext, v, e, ext);
      case HERMES_HDIV_NORM:
        return hdiv_projection_residual<Ord, Ord>(n, wt, u_ext, v, e, ext);
      default:
        error("Unknown projection type");
        return Ord();
      }
    }

  private:
    ProjNormType projNormType;

    template<typename Real, typename Scalar>
    Scalar h1_projection_residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                         Geom<Real> *e, ExtData<Scalar> *ext) const
     {
       _F_
       Scalar result = 0;
       for (int i = 0; i < n; i++)
         result += wt[i] * ((u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i] 
                            + (u_ext[this->i]->dx[i] - ext->fn[0]->dx[i]) * v->dx[i] 
                            + (u_ext[this->i]->dy[i] - ext->fn[0]->dy[i]) * v->dy[i]);
       return result;
    }

    template<typename Real, typename Scalar>
    Scalar h1_semi_projection_residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                              Geom<Real> *e, ExtData<Scalar> *ext) const
    {
      _F_
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * ((u_ext[this->i]->dx[i] - ext->fn[0]->dx[i]) * v->dx[i] 
                           + (u_ext[this->i]->dy[i] - ext->fn[0]->dy[i]) * v->dy[i]);
      return result;
    }

    template<typename Real, typename Scalar>
    Scalar l2_projection_residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                         Geom<Real> *e, ExtData<Scalar> *ext) const
     {
       _F_
       Scalar result = 0;
       for (int i = 0; i < n; i++)
         result += wt[i] * (u_ext[this->i]->val[i] - ext->fn[0]->val[i]) * v->val[i];
       return result;
    }

    template<typename Real, typename Scalar>
    Scalar hcurl_projection_residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                            Geom<Real> *e, ExtData<Scalar> *ext) const
     {
       _F_
       Scalar result = 0;
       for (int i = 0; i < n; i++) {
         result += wt[i] * (u_ext[this->i]->curl[i] - ext->fn[0]->curl[i]) * conj(v->curl[i]);
         result += wt[i] * ((u_ext[this->i]->val0[i] - ext->fn[0]->val0[i]) * conj(v->val0[i])
			                      + (u_ext[this->i]->val1[i] - ext->fn[0]->val1[i]) * conj(v->val1[i]));
       }

       return result;
    }

    template<typename Real, typename Scalar>
    Scalar hdiv_projection_residual(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                                           Geom<Real> *e, ExtData<Scalar> *ext) const
     {
       _F_
       Scalar result = 0;
       for (int i = 0; i < n; i++) {
         result += wt[i] * (u_ext[this->i]->div[i] - ext->fn[0]->div[i]) * conj(v->div[i]);
         result += wt[i] * ((u_ext[this->i]->val0[i] - ext->fn[0]->val0[i]) * conj(v->val0[i])
			                      + (u_ext[this->i]->val1[i] - ext->fn[0]->val1[i]) * conj(v->val1[i]));
       }

       return result;
    }
  };
};

#endif

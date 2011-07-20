// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
// This file was taken from hermes2d and adjusted for hermes3d
//

#include "../h3d_common.h"
#include "weakform.h"
#include "../../../hermes_common/matrix.h"
#include "../../../hermes_common/error.h"
#include "../../../hermes_common/trace.h"
#include "../../../hermes_common/callstack.h"

WeakForm::WeakForm(int neq, bool mat_free)
{
  _F_
  this->neq = neq;
  seq = 0;
  this->is_matfree = mat_free;
}

WeakForm::~WeakForm()
{
	_F_
}

void WeakForm::add_matrix_form(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, SymFlag sym, int area,
                               Hermes::vector<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) error("Invalid equation number.");
	if (sym != HERMES_ANTISYM && sym != HERMES_NONSYM && sym != HERMES_SYM) error("\"sym\" must be HERMES_ANTISYM, HERMES_NONSYM or HERMES_SYM.");
	if (sym < 0 && i == j) error("Only off-diagonal forms can be antisymmetric.");
	if (area != HERMES_ANY_INT && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");
	if (mfvol.size() > 100) warning("Large number of forms (> 100). Is this the intent?");

        MatrixFormVol form = { i, j, sym, area, fn, ord, ext };
	mfvol.push_back(form);
}

void WeakForm::add_matrix_form_surf(int i, int j, matrix_form_val_t fn, matrix_form_ord_t ord, int area, 
                                    Hermes::vector<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq || j < 0 || j >= neq) error("Invalid equation number.");
	if (area != HERMES_ANY_INT && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");

        MatrixFormSurf form = { i, j, area, fn, ord, ext };
	mfsurf.push_back(form);
}

void WeakForm::add_vector_form(int i, vector_form_val_t fn, vector_form_ord_t ord, int area, 
                               Hermes::vector<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq) error("Invalid equation number.");
	if (area != HERMES_ANY_INT && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");

        VectorFormVol form = { i, area, fn, ord, ext };
	vfvol.push_back(form);
}

void WeakForm::add_vector_form_surf(int i, vector_form_val_t fn, vector_form_ord_t ord, int area, 
                                    Hermes::vector<MeshFunction*> ext)
{
	_F_
	if (i < 0 || i >= neq) error("Invalid equation number.");
	if (area != HERMES_ANY_INT && area < 0 && -area > (signed) areas.size()) error("Invalid area number.");

        VectorFormSurf form = { i, area, fn, ord, ext };
	vfsurf.push_back(form);
}

void WeakForm::set_ext_fns(void *fn, Hermes::vector<MeshFunction*> ext)
{
	EXIT(HERMES_ERR_NOT_IMPLEMENTED);
}


//// stages ////////////////////////////////////////////////////////////////////////////////////////

/// Constructs a list of assembling stages. Each stage contains a list of forms
/// that share the same meshes. Each stage is then assembled separately. This
/// improves the performance of multi-mesh assembling.
/// This function is identical in H2D and H3D.
///
void WeakForm::get_stages(Hermes::vector<Space *> spaces, Hermes::vector<Solution *>& u_ext, 
               std::vector<WeakForm::Stage>& stages, bool want_matrix, bool want_vector)
{
  _F_

  if (!want_matrix && !want_vector) return;

  unsigned i;
  stages.clear();

  if (want_matrix || want_vector) {    // This is because of linear problems where 
                                       // matrix terms with the Dirichlet lift go to rhs.
    // process volume matrix forms
    for (i = 0; i < mfvol.size(); i++) 
    {
      int ii = mfvol[i].i, jj = mfvol[i].j;
      Mesh *m1 = spaces[ii]->get_mesh();
      Mesh *m2 = spaces[jj]->get_mesh();
      Stage *s = find_stage(stages, ii, jj, m1, m2, 
                            mfvol[i].ext, u_ext);
      s->mfvol.push_back(&mfvol[i]);
    }

    // process surface matrix forms
    for (i = 0; i < mfsurf.size(); i++) 
    {
      int ii = mfsurf[i].i, jj = mfsurf[i].j;
      Mesh *m1 = spaces[ii]->get_mesh();
      Mesh *m2 = spaces[jj]->get_mesh();
      Stage *s = find_stage(stages, ii, jj, m1, m2, 
                            mfsurf[i].ext, u_ext);
      s->mfsurf.push_back(&mfsurf[i]);
    }
  }	
		
  if (want_vector) {
    // process volume vector forms
    for (unsigned i = 0; i < vfvol.size(); i++) {
      int ii = vfvol[i].i;
      Mesh *m = spaces[ii]->get_mesh();
      Stage *s = find_stage(stages, ii, ii, m, m, 
                            vfvol[i].ext, u_ext);
      s->vfvol.push_back(&vfvol[i]);
    }

    // process surface vector forms
    for (unsigned i = 0; i < vfsurf.size(); i++) {
      int ii = vfsurf[i].i;
      Mesh *m = spaces[ii]->get_mesh();
      Stage *s = find_stage(stages, ii, ii, m, m, 
                            vfsurf[i].ext, u_ext);
      s->vfsurf.push_back(&vfsurf[i]);
    }
  }

  // helper macro for iterating in a set
  #define set_for_each(myset, type) \
    for (std::set<type>::iterator it = (myset).begin(); it != (myset).end(); it++)

  // initialize the arrays meshes and fns needed by Traverse for each stage
  for (i = 0; i < stages.size(); i++)
  {
    Stage* s = &stages[i];
    
    // First, initialize arrays for the test functions. A pointer to the PrecalcShapeset 
    // corresponding to each space will be assigned to s->fns later during assembling.
    set_for_each(s->idx_set, int)
    {
      s->idx.push_back(*it);
      s->meshes.push_back(spaces[*it]->get_mesh());
      s->fns.push_back(NULL);
    }
    
    // Next, append to the existing arrays the external functions (including the solutions
    // from previous Newton iteration) and their meshes. Also fill in a special array with
    // these external functions only.
    set_for_each(s->ext_set, MeshFunction*)
    {
      s->ext.push_back(*it);
      s->meshes.push_back((*it)->get_mesh());
      s->fns.push_back(*it);
    }
    
    s->idx_set.clear();
    s->seq_set.clear();
    s->ext_set.clear();
  }
}


/// Finds an assembling stage with the same set of meshes as [m1, m2, ext, u_ext]. If no such
/// stage can be found, a new one is created and returned.
/// This function is the same in H2D and H3D.
///
WeakForm::Stage* WeakForm::find_stage(std::vector<WeakForm::Stage>& stages, int ii, int jj,
                                      Mesh* m1, Mesh* m2, 
                                      std::vector<MeshFunction*>& ext, Hermes::vector<Solution*>& u_ext)
{
  _F_
  // first create a list of meshes the form uses
  std::set<unsigned> seq;
  seq.insert(m1->get_seq());
  seq.insert(m2->get_seq());
  Mesh *mmm;
  for (unsigned i = 0; i < ext.size(); i++) {
    mmm = ext[i]->get_mesh();
    if (mmm == NULL) error("NULL Mesh pointer detected in ExtData during assembling.\n  Have you initialized all external functions?");
    seq.insert(mmm->get_seq());
  }
  for (unsigned i = 0; i < u_ext.size(); i++) {
    if (u_ext[i] != NULL) {
      mmm = u_ext[i]->get_mesh();
      if (mmm == NULL) error("NULL Mesh pointer detected in u_ext during assembling.");
      seq.insert(mmm->get_seq());
    }
  }
  
  // find a suitable existing stage for the form
  Stage* s = NULL;
  for (unsigned i = 0; i < stages.size(); i++)
    if (seq.size() == stages[i].seq_set.size() &&
        equal(seq.begin(), seq.end(), stages[i].seq_set.begin()))
      { s = &stages[i]; break; }

  // create a new stage if not found
  if (s == NULL)
  {
    Stage newstage;
    stages.push_back(newstage);
    s = &stages.back();
    s->seq_set = seq;
  }

  // update and return the stage
  for (unsigned int i = 0; i < ext.size(); i++)
    s->ext_set.insert(ext[i]);
  for (unsigned int i = 0; i < u_ext.size(); i++)
    if (u_ext[i] != NULL)
      s->ext_set.insert(u_ext[i]);
  
  s->idx_set.insert(ii);
  s->idx_set.insert(jj);
  return s;
}

/// Returns a (neq x neq) array containing true in each element, if the corresponding
/// block of weak forms is used, and false otherwise.
/// This function is the same in H2D and H3D.
///
bool **WeakForm::get_blocks()
{
  _F_
  bool** blocks = new_matrix<bool>(neq, neq);
  for (int i = 0; i < neq; i++)
    for (int j = 0; j < neq; j++)
      blocks[i][j] = false;

  for (unsigned i = 0; i < mfvol.size(); i++) {
    blocks[mfvol[i].i][mfvol[i].j] = true;
    if (mfvol[i].sym)
      blocks[mfvol[i].j][mfvol[i].i] = true;
  }

  for (unsigned i = 0; i < mfsurf.size(); i++)
    blocks[mfsurf[i].i][mfsurf[i].j] = true;

  return blocks;
}

//// areas /////////////////////////////////////////////////////////////////////////////////////////

int WeakForm::def_area(Hermes::vector<int> area_markers)
{
	_F_
	Area newarea;
        int n = area_markers.size();
	for (int i = 0; i < n; i++) newarea.markers.push_back(area_markers[i]);

	areas.push_back(newarea);
	return -areas.size();
}


bool WeakForm::is_in_area_2(int marker, int area) const
{
	_F_
	if (-area > (signed) areas.size()) error("Invalid area number.");
	const Area *a = &areas[-area - 1];

	for (unsigned i = 0; i < a->markers.size(); i++)
		if (a->markers[i] == marker)
			return true;

	return false;
}

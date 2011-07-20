// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef _Space_H_
#define _Space_H_

#include "../../hermes_common/common.h"
#include "../../hermes_common/matrix.h"
#include "legendre.h"
#include "lobatto.h"

class HERMES_API Space;

class HERMES_API Element {
public:
    Element();
    Element(double x_left, double x_right, int level, int deg, 
            int n_eq, int n_sln, int marker);
    void free_element();

    ~Element();

    void init(double x1, double x2, int p_init, 
	      int id, int active, int level, int n_eq, int n_sln, int marker);
    void copy_into(Element *e_trg);
    void copy_recursively_into(Element *e_trg);
    double get_x_phys(double x_ref);             // Gets physical coordinate of a reference point.
    double calc_elem_norm_squared(int norm);
    void get_coeffs_from_vector(double *y, int sln=0);
    void copy_coeffs_to_vector(double *y, int sln=0);
    void copy_dofs(int sln_src, int sln_trg);
    void get_solution_quad(int flag, int quad_order, 
                           double val_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
			   double der_phys[MAX_EQN_NUM][MAX_QUAD_PTS_NUM], int sln=0);
    void get_solution_plot(double x_phys[MAX_PLOT_PTS_NUM], int pts_num,
                           double val_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], 
			   double der_phys[MAX_EQN_NUM][MAX_PLOT_PTS_NUM], int sln=0);
    double get_solution_value(double x_phys, int c);
    double get_solution_deriv(double x_phys, int c);
    void get_coeffs(int sln, int comp, double coeffs[]);
    void get_solution_point(double x_phys, 
			    double val[MAX_EQN_NUM], double der[MAX_EQN_NUM], int sln=0);
    int create_cand_list(int adapt_type, int p_ref_left, int p_ref_right, int3 *cand_list);
    void print_cand_list(int num_cand, int3 *cand_list);
    void refine(int3 cand);
    void refine(int type, int p_left, int p_right);
    unsigned is_active();
    unsigned active;   // flag used by assembling algorithm
    double x1, x2;     // endpoints
    int p;             // poly degrees
    int marker;        // can be used to distinguish between material parameters
    int n_eq;          // number of equations (= number of solution components)
    int n_sln;         // number of solution copies
    int dof[MAX_EQN_NUM][MAX_P + 1];   // connectivity array of length p+1 
                                       // for every solution component
    double coeffs[MAX_SLN_NUM][MAX_EQN_NUM][MAX_P + 1];   // solution coefficient array of length p+1 
                                                          // for every component and every solution 
    int id;
    unsigned level;    // refinement level (zero for initial space elements) 
    Element *sons[2];  // for refinement
};

typedef Element* ElemPtr2[2];

void HERMES_API get_coeff_vector(Space *space, scalar *y, int sln = 0);
void HERMES_API set_coeff_vector(scalar *y, Space *space, int sln = 0);
void HERMES_API get_coeff_vector(Space *space, Vector *y, int sln = 0);
void HERMES_API set_coeff_vector(Vector *y, Space *space, int sln =  0);

class HERMES_API Space {
    public:
        Space();
        // Creates equidistant space with uniform polynomial degree of elements.
        // All elements will have the same (zero) marker.
        Space(double a, double b, int n_base_elem, int p_init = 1, int n_eq = 1, 
              int n_sln = 1, bool print_banner = true);
        Space(double a, double b, int n_base_elem, 
              Hermes::vector<std::pair<int, double> *> left_boundary_conditions = 
                ( Hermes::vector<std::pair<int, double> *>() ), 
              Hermes::vector<std::pair<int, double> *> right_boundary_conditions = 
                ( Hermes::vector<std::pair<int, double> *>() ), 
              int p_init = 1, 
              int n_eq = 1, 
              int n_sln = 1, 
              bool print_banner = true);

        Space(double a, double b, int n_base_elem, std::pair<int, double> left_boundary_condition, 
              std::pair<int, double> right_boundary_condition, int p_init = 1, 
              int n_eq = 1, 
              int n_sln = 1, 
              bool print_banner = true);
        // Creates a general space (used, e.g., in example "neutronics").
        // n_macro_elem... number of macro elements
        // pts_array[]...  array of macroelement grid points
        // p_array[]...    array of macroelement poly degrees
        // m_array[]...    array of macroelement material markers
        // div_array[]...  array of macroelement equidistant divisions
        Space(int n_macro_elem, 
              double *pts_array, int *p_array, int *m_array, int *div_array, 
              Hermes::vector<std::pair<int, double> *> left_boundary_conditions = 
                ( Hermes::vector<std::pair<int, double> *>() ), 
              Hermes::vector<std::pair<int, double> *> right_boundary_conditions = 
                ( Hermes::vector<std::pair<int, double> *>() ), 
              int n_eq=1, 
              int n_sln=1, 
              bool print_banner=true);
        
        Space(int n_macro_elem, 
              double *pts_array, int *p_array, int *m_array, int *div_array, 
              std::pair<int, double> left_boundary_condition, 
              std::pair<int, double> right_boundary_condition, 
              int n_eq=1, 
              int n_sln=1, 
              bool print_banner=true);
        //FIXME: We currently do not have a Python wrapper for Hermes::vector, so we
        // cannot wrap the above constructor. This constructor is required to ensure
        // compatibility with the wrappers, but should not be used.
        Space(int n_macro_elem, 
              double *pts_array, int *p_array, int *m_array, int *div_array, 
              int n_eq=1, 
              int n_sln=1, 
              bool print_banner=true);
        
        void init(double a, double b, int n_base_elem, int p_init, int n_eq, 
                  int n_sln, bool print_banner);
        ~Space();

        void free_elements();

        int assign_dofs();

        Element *get_base_elems();

        int get_n_base_elem();

        void set_n_base_elem(int n_base_elem);

        int get_n_active_elem();

        void set_n_active_elem(int n);

        static int get_num_dofs(Space* space);

        int get_num_dofs() { return Space::get_num_dofs(this); };

        void set_n_dof(int n);

        int get_n_eq();

        void set_n_eq(int n_eq);

        int get_n_sln();

        void set_n_sln(int n_sln);

        double get_left_endpoint();

        void set_left_endpoint(double a);
        
        double get_right_endpoint();

        void set_right_endpoint(double b);

        Element* first_active_element();
        Element* last_active_element();
        void set_bc_left_dirichlet(int eqn, double val);
        void set_bc_right_dirichlet(int eqn, double val);
        void refine_single_elem(int id, int3 cand);
        void refine_elems(int elem_num, int *id_array, int3 *cand_array);
        void reference_refinement(int start_elem_id, int elem_num);
        Space *replicate(); 
        void plot(const char* filename); // plots the space and polynomial degrees of elements
        void plot_element_error_p(int norm, FILE *f, Element *p, Element *e_ref,  
                                  int subdivision = 20); // plots error wrt. reference solution
        void plot_element_error_hp(int norm, FILE *f, Element *p, 
                                   Element *e_ref_left, Element *e_ref_right, 
                                   int subdivision = 20); // plots error wrt. reference solution
                                                          // if ref. refinement was hp-refinement
        void plot_element_error_exact(int norm, FILE *f, Element *p, 
				   exact_sol_type exact_sol,
                                   int subdivision = 20); // plots error wrt. reference solution
                                                          // if ref. refinement was hp-refinement
        void plot_error_estimate(int norm, Space* space_ref, const char *filename, 
                        int subdivision = 500);  // plots error wrt. reference solution
        void plot_error_estimate(int norm, ElemPtr2* elem_ref_pairs, const char *filename,  
                        int subdivision = 500);  // plots error wrt. reference solution
        void plot_error_exact(int norm, exact_sol_type exact_sol, const char *filename,  
                        int subdivision = 500); // plots error wrt. exact solution
        void assign_elem_ids();
        int n_active_elem;
        void set_coeff_vector(scalar *y, int sln=0);

        void get_coeff_vector(scalar *y, int sln=0);

    private:
        double left_endpoint, right_endpoint;
        int n_eq;            // number of equations in the system
        int n_sln;           // number of solution copies
        int n_base_elem;     // number of elements in the base space
        int n_dof;           // number of DOF (in each solution copy)
        Element *base_elems; // base space

};

// Returns updated coarse and reference spacees, with the last 
// coarse and reference space solutions on them, respectively. 
// The coefficient vectors and numbers of degrees of freedom 
// on both spacees are also updated. 
// Refine coarse mesh elements whose id_array >= 0, and 
// adjust the reference space accordingly.  
// Returns updated coarse and reference spacees, with the last 
// coarse and reference space solutions on them, respectively. 
// The coefficient vectors and numbers of degrees of freedom 
// on both spacees are also updated. 

void HERMES_API adapt(int norm, int adapt_type, double threshold, 
           double *err_squared_array,
           Space* & space, Space* & space_ref);

// Returns updated coarse mesh, with the last 
// coarse solution on it. 
// The coefficient vector and number of degrees of freedom 
// also is updated. 

void HERMES_API adapt(int norm, int adapt_type, double threshold, 
           double *err_array, 
           Space* &space, ElemPtr2 *ref_elem_pairs);

void HERMES_API adapt_plotting(Space *space, Space *space_ref,
                    int norm, int exact_sol_provided, 
                    exact_sol_type exact_sol); 

void HERMES_API adapt_plotting(Space *space, ElemPtr2* ref_elem_pairs,
                    int norm, int exact_sol_provided, 
                    exact_sol_type exact_sol); 
#endif

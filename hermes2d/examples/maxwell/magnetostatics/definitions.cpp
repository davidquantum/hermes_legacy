#include "weakform/weakform.h"
#include "weakform_library/h1.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

using namespace WeakFormsH1;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

/* Weak forms */

class CustomWeakFormMagnetostatics : public WeakForm
{
public:
  CustomWeakFormMagnetostatics(std::string material_iron_1, std::string material_iron_2, 
                               CubicSpline* mu_inv_iron, std::string material_air,
                               std::string material_copper, double mu_vacuum, double current_density) 
    : WeakForm(1) {

    // Jacobian.
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_air, 1.0, HERMES_SYM, HERMES_AXISYM_Y));
    add_matrix_form(new DefaultLinearMagnetostatics(0, 0, material_copper, 1.0, HERMES_SYM, HERMES_AXISYM_Y));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, material_iron_1, mu_inv_iron, HERMES_SYM, HERMES_AXISYM_Y));
    add_matrix_form(new DefaultJacobianNonlinearMagnetostatics(0, 0, material_iron_2, mu_inv_iron, HERMES_SYM, HERMES_AXISYM_Y));

    // Residual.
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_air, 1.0, HERMES_AXISYM_Y));
    add_vector_form(new DefaultResidualLinearMagnetostatics(0, material_copper, 1.0, HERMES_AXISYM_Y));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, material_iron_1, mu_inv_iron, HERMES_AXISYM_Y));
    add_vector_form(new DefaultResidualNonlinearMagnetostatics(0, material_iron_2, mu_inv_iron, HERMES_AXISYM_Y));
    add_vector_form(new DefaultVectorFormConst(0, material_copper, -current_density * mu_vacuum, HERMES_AXISYM_Y));
  };
};




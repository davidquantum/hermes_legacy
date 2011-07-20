#include "definitions.h"

CustomWeakFormPoissonNewton::CustomWeakFormPoissonNewton(double lambda, double alpha, double T0, 
                                                         std::string bdy_heat_flux) : WeakForm(1)
{
  // Jacobian form - volumetric.
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0, HERMES_ANY, new HermesFunction(lambda),
                                                            HERMES_SYM, HERMES_AXISYM_Y));

  // Jacobian form - surface.
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf(0, 0, bdy_heat_flux, new HermesFunction(alpha),
                                                              HERMES_AXISYM_Y));

  // Residual forms - volumetric.
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion(0, HERMES_ANY, new HermesFunction(lambda),
                                                            HERMES_AXISYM_Y));

  // Residual form - surface.
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf(0, bdy_heat_flux, new HermesFunction(alpha),
                                                            HERMES_AXISYM_Y));
  add_vector_form_surf(new WeakFormsH1::DefaultVectorFormSurf(0, bdy_heat_flux, new HermesFunction(-alpha * T0),
                                                              HERMES_AXISYM_Y));
};

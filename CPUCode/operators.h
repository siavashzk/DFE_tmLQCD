#ifndef DFE_OPERATORS_H
#define DFE_OPERATORS_H

#include "su3.h"

void dfe_Qtm_pm_psi(spinor * const l, spinor * const k);
void dfe_H_eo_tm_inv_psi(spinor * const l, spinor * const k,
                                const int ieo, const double _sign);
void dfe_Qtm_minus_psi(spinor * const l, spinor * const k);


#endif

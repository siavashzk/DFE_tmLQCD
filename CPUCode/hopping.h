#ifndef HOPPING_H
#define HOPPING_H

void tm_times_Hopping_Matrix(spinor * const l, spinor * const k,
                                 su3 * U0, su3 * U1, complex float const cfactor,
                                 int ieo);
void tm_sub_Hopping_Matrix(spinor * const l, spinor * const p, spinor * const k,
                                 su3 * U0, su3 * U1, complex float const cfactor,
                                 int ieo);
void Qtm_pm_psi(spinor * const l, spinor * const k, su3 * UOdd, su3 * UEven,
                complex float cfactor0, complex float cfactor1);

#endif

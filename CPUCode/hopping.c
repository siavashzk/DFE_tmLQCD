/*
 * hopping.c
 *
 *  Created on: 21 May 2015
 *      Author: skamali
 */
#include <stdio.h>

#include "su3.h"
#include "halfspinor_hopping.h"
#include "global.h"

void tm_times_Hopping_Matrix(spinor * const l, spinor * const k, su3 * U0,
		su3 * U1, complex float const cfactor, int ieo) {

#define _MUL_G5_CMPLX
#  include "halfspinor_body.c"
#undef _MUL_G5_CMPLX

	return;
}

void tm_sub_Hopping_Matrix(spinor * const l, spinor * const p, spinor * const k, su3 * U0,
		su3 * U1, complex float const cfactor, int ieo) {

	spinor * pn;
#define _TM_SUB_HOP
#  include "halfspinor_body.c"
#undef _TM_SUB_HOP

	return;
}

void Qtm_pm_psi(spinor * const l, spinor * const k, su3 * UOdd, su3 * UEven,
                complex float cfactor0, complex float cfactor1) {
  spinor temp[VOLUME/2];
  spinor temp2[VOLUME/2];

  tm_times_Hopping_Matrix(l, k, UOdd, UEven, cfactor0, 0);
  /*tm_sub_Hopping_Matrix(temp2, k, temp, UEven, UOdd, cfactor1, 1);
  tm_times_Hopping_Matrix(temp, temp2, UOdd, UEven, -cfactor0, 0);
  tm_sub_Hopping_Matrix(l, temp2, temp, UEven, UOdd, -cfactor1, 1);*/

}

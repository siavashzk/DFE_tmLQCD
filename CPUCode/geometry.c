/*
 * geometry.c
 *
 *  Created on: 21 May 2015
 *      Author: skamali
 */

#include <printf.h>
#include <stdlib.h>

#include "su3.h"
#include "global.h"

static int *iup = NULL, *idn = NULL;

static halfspinor ** NBPointer_;
static halfspinor * HalfSpinor_;
static halfspinor * HalfSpinor;

int init_geometry_indices(const int V);
int Index(const int x0, const int x1, const int x2, const int x3);
void geometry();
int init_dirac_halfspinor();

/**************** Public Function Definitions ****************/
void init_geometry(void) {
	init_geometry_indices(VOLUME);
	geometry();
	init_dirac_halfspinor();
}

/**************** Private Function Definitions ****************/
int init_geometry_indices(const int V) {
	int i = 0;

	g_idn = (int**) calloc(V, sizeof(int*));
	if ((void*) g_idn == NULL)
		return (1);
	g_iup = (int**) calloc(V, sizeof(int*));
	if ((void*) g_iup == NULL)
		return (2);

	idn = (int*) calloc(4 * V, sizeof(int));
	if ((void*) idn == NULL)
		return (6);
	iup = (int*) calloc(4 * V, sizeof(int));
	if ((void*) iup == NULL)
		return (7);

	g_lexic2eo = (int*) calloc(V, sizeof(int));
	if ((void*) g_lexic2eo == NULL)
		return (9);
	/* this +2 is for sanity reasons */
	g_lexic2eosub = (int*) calloc(V + 2, sizeof(int));
	if ((void*) g_lexic2eosub == NULL)
		return (10);
	g_eo2lexic = (int*) calloc(V, sizeof(int));
	if ((void*) g_eo2lexic == NULL)
		return (11);

	g_idn[0] = idn;
	g_iup[0] = iup;

	for (i = 1; i < V; i++) {
		g_idn[i] = g_idn[i - 1] + 4;
		g_iup[i] = g_iup[i - 1] + 4;
	}

	return (0);
}

int Index(const int x0, const int x1, const int x2, const int x3) {
	int y0, y1, y2, y3, ix;

	y0 = (x0 + T) % T;

	y1 = (x1 + LX) % LX;
	y2 = (x2 + LY) % LY;
	y3 = (x3 + LZ) % LZ;
	ix = ((y0 * LX + y1) * LY + y2) * LZ + y3;

	return (ix);
}

void geometry() {

	int x0, x1, x2, x3, ix;
	int i_even, i_odd;
	int * xeven;

	xeven = malloc(VOLUME * sizeof(int));

	/* extended for boundary slices */
	for (x0 = 0; x0 < T; x0++) {
		for (x1 = 0; x1 < LX; x1++) {
			for (x2 = 0; x2 < LY; x2++) {
				for (x3 = 0; x3 < LZ; x3++) {
					ix = Index(x0, x1, x2, x3);

					if ((x0 + x1 + x2 + x3) % 2 == 0) {
						xeven[ix] = 1;
					} else {
						xeven[ix] = 0;
					}
					g_iup[ix][0] = Index(x0 + 1, x1, x2, x3);
					g_idn[ix][0] = Index(x0 - 1, x1, x2, x3);

					g_iup[ix][1] = Index(x0, x1 + 1, x2, x3);
					g_idn[ix][1] = Index(x0, x1 - 1, x2, x3);

					g_iup[ix][2] = Index(x0, x1, x2 + 1, x3);
					g_idn[ix][2] = Index(x0, x1, x2 - 1, x3);

					g_iup[ix][3] = Index(x0, x1, x2, x3 + 1);
					g_idn[ix][3] = Index(x0, x1, x2, x3 - 1);
				}
			}
		}
	}

	i_even = 0;
	i_odd = 0;
	/*For the spinor fields we need only till VOLUME+RAND */
	for (ix = 0; ix < (VOLUME); ix++) {
		if (xeven[ix] == 1) {
			g_lexic2eo[ix] = i_even;
			g_lexic2eosub[ix] = i_even;
			g_eo2lexic[i_even] = ix;
			i_even++;
		} else {
			g_lexic2eo[ix] = (VOLUME) / 2 + i_odd;
			g_lexic2eosub[ix] = i_odd;
			g_eo2lexic[(VOLUME) / 2 + i_odd] = ix;
			i_odd++;
		}
	}

	free(xeven);
}

int init_dirac_halfspinor() {
	int j = 0;

	NBPointer = (halfspinor***) calloc(4, sizeof(halfspinor**));
	NBPointer_ = (halfspinor**) calloc(16, (VOLUME) * sizeof(halfspinor*));
	NBPointer[0] = NBPointer_;
	NBPointer[1] = NBPointer_ + (8 * (VOLUME) / 2);
	NBPointer[2] = NBPointer_ + (16 * (VOLUME) / 2);
	NBPointer[3] = NBPointer_ + (24 * (VOLUME) / 2);

	if ((void*) (HalfSpinor_ = (halfspinor*) calloc(4 * (VOLUME) + 1,
			sizeof(halfspinor))) == NULL) {
		return (1);
	}

	HalfSpinor = (halfspinor*) (((unsigned long int) (HalfSpinor_) + 0xf + 1)
			& ~0xf);

	for (int ieo = 0; ieo < 2; ieo++) {
		for (int i = 0; i < VOLUME / 2; i++) {
			j = g_eo2lexic[i + ((ieo + 1) % 2) * (VOLUME) / 2];
			for (int mu = 0; mu < 4; mu++) {
				NBPointer[ieo][8 * i + 2 * mu + 0] = &HalfSpinor[8 * i + 2 * mu
						+ 0];
				NBPointer[ieo][8 * i + 2 * mu + 1] = &HalfSpinor[8 * i + 2 * mu
						+ 1];
			}
		}
	}
	for (int ieo = 2; ieo < 4; ieo++) {
		for (int i = 0; i < VOLUME / 2; i++) {
			for (int mu = 0; mu < 4; mu++) {
				j = g_eo2lexic[i + ((ieo) % 2) * (VOLUME) / 2];
				NBPointer[ieo][8 * i + 2 * mu + 0] = &HalfSpinor[8
						* (g_iup[j][mu] / 2) + 2 * mu + 0];
				NBPointer[ieo][8 * i + 2 * mu + 1] = &HalfSpinor[8
						* (g_idn[j][mu] / 2) + 2 * mu + 1];

			}
		}
	}
	return (0);

}


#ifndef GLOBAL_H
#define GLOBAL_H

#ifdef TEST_MAIN
#  define EXTERN
#else
#  define EXTERN extern
#endif

EXTERN int T, VOLUME;
EXTERN int LX, LY, LZ;

EXTERN double ka0, ka1, ka2, ka3;

EXTERN int **g_iup, **g_idn;

EXTERN int * g_eo2lexic;
EXTERN int * g_lexic2eo;
EXTERN int * g_lexic2eosub;


EXTERN halfspinor *** NBPointer;

#endif

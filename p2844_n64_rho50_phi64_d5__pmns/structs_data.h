#ifndef STRUCTS_DATA
#define STRUCTS_DATA


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>

#define PHI_LOG2 64
#define POLY_DEG 63
#define NB_COEFF 64
#define NB_ADD_MAX 5

#define RHO_LOG2 50  // rho = 1 << RHO_LOG2.

#define CONV_MASK 1125899906842623UL  // CONV_MASK = (1 << RHO_LOG2) - 1, for conversion

typedef __int128 int128;



#endif


#ifndef STRUCTS_DATA
#define STRUCTS_DATA

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <gmp.h>


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 4
#define NB_COEFF 5
#define NB_ADD_MAX 0

#define RHO_LOG2 55  // rho = 1 << RHO_LOG2.

#define BETA_LOG2 64  // beta = 1 << BETA_LOG2, for conversion

typedef __int128 int128;


#endif


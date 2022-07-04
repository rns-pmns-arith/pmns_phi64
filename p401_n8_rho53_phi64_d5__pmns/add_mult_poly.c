#include "p__conv_data.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

//~ computes : pa + 2.pb
void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + 2*pb[j];
}

//~ computes : pa - 2.pb
void double_sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - 2*pb[j];
}

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void double_poly_coeffs(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << 1;
}

//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.
void lshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << nb_pos;
}

//~ Computes pa(X)*pb(X) mod(E)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];
	int128 c8, c9, c10, c11, c12, c13, c14;

	c8 = (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1];
	c9 = (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2];
	c10 = (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3];
	c11 = (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4];
	c12 = (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5];
	c13 = (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6];
	c14 = (int128)pa[7] * pb[7];

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + c8;
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - 2*c8 + c9;
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - 2*c9 + c10;
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - 2*c10 + c11;
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - 2*c11 + c12;
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - 2*c12 + c13;
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - 2*c13 + c14;
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - 2*c14;

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];
	int128 c8, c9, c10, c11, c12, c13, c14;

	c8 = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4];
	c9 = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 1);
	c10 = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5];
	c11 = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 1);
	c12 = (((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6];
	c13 = (((int128)pa[7] * pa[6]) << 1);
	c14 = (int128)pa[7] * pa[7];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + c8;
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - 2*c8 + c9;
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - 2*c9 + c10;
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - 2*c10 + c11;
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - 2*c11 + c12;
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - 2*c12 + c13;
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - 2*c13 + c14;
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - 2*c14;

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmpQ[8];
	int128 tmpZero[8];

	//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)
	tmpQ[0] = (uint64_t)op[0] * 0x26eee722adf556f4UL + (uint64_t)op[1] * 0x241aee20b221d85aUL + (uint64_t)op[2] * 0x4e8c51c94e887851UL + (uint64_t)op[3] * 0x7464ce96e327dc8eUL + (uint64_t)op[4] * 0x72e924f9e5115d54UL + (uint64_t)op[5] * 0x668f31604a8407b6UL + (uint64_t)op[6] * 0x81088b477497f1eeUL + (uint64_t)op[7] * 0xc4f6a21b8a9058e2UL;
	tmpQ[1] = (uint64_t)op[0] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[1] * 0xdeb90ae149b1a640UL + (uint64_t)op[2] * 0x87024a8e1510e7b8UL + (uint64_t)op[3] * 0x65c2b49b8838bf35UL + (uint64_t)op[4] * 0x8e9284a3190521e6UL + (uint64_t)op[5] * 0xa5cac23950094de8UL + (uint64_t)op[6] * 0x647e1ad1615423daUL + (uint64_t)op[7] * 0xf71b47105f77402aUL;
	tmpQ[2] = (uint64_t)op[0] * 0x81088b477497f1eeUL + (uint64_t)op[1] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[2] * 0xdeb90ae149b1a640UL + (uint64_t)op[3] * 0x87024a8e1510e7b8UL + (uint64_t)op[4] * 0x65c2b49b8838bf35UL + (uint64_t)op[5] * 0x8e9284a3190521e6UL + (uint64_t)op[6] * 0xa5cac23950094de8UL + (uint64_t)op[7] * 0x647e1ad1615423daUL;
	tmpQ[3] = (uint64_t)op[0] * 0x668f31604a8407b6UL + (uint64_t)op[1] * 0x81088b477497f1eeUL + (uint64_t)op[2] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[3] * 0xdeb90ae149b1a640UL + (uint64_t)op[4] * 0x87024a8e1510e7b8UL + (uint64_t)op[5] * 0x65c2b49b8838bf35UL + (uint64_t)op[6] * 0x8e9284a3190521e6UL + (uint64_t)op[7] * 0xa5cac23950094de8UL;
	tmpQ[4] = (uint64_t)op[0] * 0x72e924f9e5115d54UL + (uint64_t)op[1] * 0x668f31604a8407b6UL + (uint64_t)op[2] * 0x81088b477497f1eeUL + (uint64_t)op[3] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[4] * 0xdeb90ae149b1a640UL + (uint64_t)op[5] * 0x87024a8e1510e7b8UL + (uint64_t)op[6] * 0x65c2b49b8838bf35UL + (uint64_t)op[7] * 0x8e9284a3190521e6UL;
	tmpQ[5] = (uint64_t)op[0] * 0x7464ce96e327dc8eUL + (uint64_t)op[1] * 0x72e924f9e5115d54UL + (uint64_t)op[2] * 0x668f31604a8407b6UL + (uint64_t)op[3] * 0x81088b477497f1eeUL + (uint64_t)op[4] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[5] * 0xdeb90ae149b1a640UL + (uint64_t)op[6] * 0x87024a8e1510e7b8UL + (uint64_t)op[7] * 0x65c2b49b8838bf35UL;
	tmpQ[6] = (uint64_t)op[0] * 0x4e8c51c94e887851UL + (uint64_t)op[1] * 0x7464ce96e327dc8eUL + (uint64_t)op[2] * 0x72e924f9e5115d54UL + (uint64_t)op[3] * 0x668f31604a8407b6UL + (uint64_t)op[4] * 0x81088b477497f1eeUL + (uint64_t)op[5] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[6] * 0xdeb90ae149b1a640UL + (uint64_t)op[7] * 0x87024a8e1510e7b8UL;
	tmpQ[7] = (uint64_t)op[0] * 0x241aee20b221d85aUL + (uint64_t)op[1] * 0x4e8c51c94e887851UL + (uint64_t)op[2] * 0x7464ce96e327dc8eUL + (uint64_t)op[3] * 0x72e924f9e5115d54UL + (uint64_t)op[4] * 0x668f31604a8407b6UL + (uint64_t)op[5] * 0x81088b477497f1eeUL + (uint64_t)op[6] * 0xc4f6a21b8a9058e2UL + (uint64_t)op[7] * 0xdeb90ae149b1a640UL;

	//~ computation of : tmp_q*red_int_coeff mod(E)
	tmpZero[0] = (int128)tmpQ[0] * 0x2bcf3242a2d44L + (int128)tmpQ[1] * 0x9b49a8cf2486L + (int128)tmpQ[2] * 0x27c9efdb0cfdeL + (int128)tmpQ[3] * 0x1a8a150e0b6faL + (int128)tmpQ[4] * 0x1fda11e0fc6e4L + (int128)tmpQ[5] * 0x1dd73b7ceba00L + (int128)tmpQ[6] * 0x11a52f86e9387L + (int128)tmpQ[7] * 0x46c583c95a56L;
	tmpZero[1] = (int128)tmpQ[0] * 0x46c583c95a56L + (int128)tmpQ[1] * 0x1865fd28be438L + (int128)tmpQ[2] * -0x45df452927b36L + (int128)tmpQ[3] * -0xd4a3a4109e16L + (int128)tmpQ[4] * -0x252a0eb3ed6ceL + (int128)tmpQ[5] * -0x1bd46518dad1cL + (int128)tmpQ[6] * -0x5732390e6d0eL + (int128)tmpQ[7] * 0x8cc7f0dbdedbL;
	tmpZero[2] = (int128)tmpQ[0] * 0x11a52f86e9387L + (int128)tmpQ[1] * 0x46c583c95a56L + (int128)tmpQ[2] * 0x1865fd28be438L + (int128)tmpQ[3] * -0x45df452927b36L + (int128)tmpQ[4] * -0xd4a3a4109e16L + (int128)tmpQ[5] * -0x252a0eb3ed6ceL + (int128)tmpQ[6] * -0x1bd46518dad1cL + (int128)tmpQ[7] * -0x5732390e6d0eL;
	tmpZero[3] = (int128)tmpQ[0] * 0x1dd73b7ceba00L + (int128)tmpQ[1] * 0x11a52f86e9387L + (int128)tmpQ[2] * 0x46c583c95a56L + (int128)tmpQ[3] * 0x1865fd28be438L + (int128)tmpQ[4] * -0x45df452927b36L + (int128)tmpQ[5] * -0xd4a3a4109e16L + (int128)tmpQ[6] * -0x252a0eb3ed6ceL + (int128)tmpQ[7] * -0x1bd46518dad1cL;
	tmpZero[4] = (int128)tmpQ[0] * 0x1fda11e0fc6e4L + (int128)tmpQ[1] * 0x1dd73b7ceba00L + (int128)tmpQ[2] * 0x11a52f86e9387L + (int128)tmpQ[3] * 0x46c583c95a56L + (int128)tmpQ[4] * 0x1865fd28be438L + (int128)tmpQ[5] * -0x45df452927b36L + (int128)tmpQ[6] * -0xd4a3a4109e16L + (int128)tmpQ[7] * -0x252a0eb3ed6ceL;
	tmpZero[5] = (int128)tmpQ[0] * 0x1a8a150e0b6faL + (int128)tmpQ[1] * 0x1fda11e0fc6e4L + (int128)tmpQ[2] * 0x1dd73b7ceba00L + (int128)tmpQ[3] * 0x11a52f86e9387L + (int128)tmpQ[4] * 0x46c583c95a56L + (int128)tmpQ[5] * 0x1865fd28be438L + (int128)tmpQ[6] * -0x45df452927b36L + (int128)tmpQ[7] * -0xd4a3a4109e16L;
	tmpZero[6] = (int128)tmpQ[0] * 0x27c9efdb0cfdeL + (int128)tmpQ[1] * 0x1a8a150e0b6faL + (int128)tmpQ[2] * 0x1fda11e0fc6e4L + (int128)tmpQ[3] * 0x1dd73b7ceba00L + (int128)tmpQ[4] * 0x11a52f86e9387L + (int128)tmpQ[5] * 0x46c583c95a56L + (int128)tmpQ[6] * 0x1865fd28be438L + (int128)tmpQ[7] * -0x45df452927b36L;
	tmpZero[7] = (int128)tmpQ[0] * 0x9b49a8cf2486L + (int128)tmpQ[1] * 0x27c9efdb0cfdeL + (int128)tmpQ[2] * 0x1a8a150e0b6faL + (int128)tmpQ[3] * 0x1fda11e0fc6e4L + (int128)tmpQ[4] * 0x1dd73b7ceba00L + (int128)tmpQ[5] * 0x11a52f86e9387L + (int128)tmpQ[6] * 0x46c583c95a56L + (int128)tmpQ[7] * 0x1865fd28be438L;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> PHI_LOG2;
	rop[1] = (op[1] + tmpZero[1]) >> PHI_LOG2;
	rop[2] = (op[2] + tmpZero[2]) >> PHI_LOG2;
	rop[3] = (op[3] + tmpZero[3]) >> PHI_LOG2;
	rop[4] = (op[4] + tmpZero[4]) >> PHI_LOG2;
	rop[5] = (op[5] + tmpZero[5]) >> PHI_LOG2;
	rop[6] = (op[6] + tmpZero[6]) >> PHI_LOG2;
	rop[7] = (op[7] + tmpZero[7]) >> PHI_LOG2;
}

void exact_coeffs_reduction(int64_t *rop, int64_t *op){

	int i;
	int128 tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp[i] = (int128) op[i];

	internal_reduction(rop, tmp);
	int128 c8, c9, c10, c11, c12, c13, c14;

	c8 = (int128)rop[1] * polys_P[0][7] + (int128)rop[2] * polys_P[0][6] + (int128)rop[3] * polys_P[0][5] + (int128)rop[4] * polys_P[0][4] + (int128)rop[5] * polys_P[0][3] + (int128)rop[6] * polys_P[0][2] + (int128)rop[7] * polys_P[0][1];
	c9 = (int128)rop[2] * polys_P[0][7] + (int128)rop[3] * polys_P[0][6] + (int128)rop[4] * polys_P[0][5] + (int128)rop[5] * polys_P[0][4] + (int128)rop[6] * polys_P[0][3] + (int128)rop[7] * polys_P[0][2];
	c10 = (int128)rop[3] * polys_P[0][7] + (int128)rop[4] * polys_P[0][6] + (int128)rop[5] * polys_P[0][5] + (int128)rop[6] * polys_P[0][4] + (int128)rop[7] * polys_P[0][3];
	c11 = (int128)rop[4] * polys_P[0][7] + (int128)rop[5] * polys_P[0][6] + (int128)rop[6] * polys_P[0][5] + (int128)rop[7] * polys_P[0][4];
	c12 = (int128)rop[5] * polys_P[0][7] + (int128)rop[6] * polys_P[0][6] + (int128)rop[7] * polys_P[0][5];
	c13 = (int128)rop[6] * polys_P[0][7] + (int128)rop[7] * polys_P[0][6];
	c14 = (int128)rop[7] * polys_P[0][7];

	tmp[0] = (int128)rop[0] * polys_P[0][0] + c8;
	tmp[1] = (int128)rop[0] * polys_P[0][1] + (int128)rop[1] * polys_P[0][0] - 2*c8 + c9;
	tmp[2] = (int128)rop[0] * polys_P[0][2] + (int128)rop[1] * polys_P[0][1] + (int128)rop[2] * polys_P[0][0] - 2*c9 + c10;
	tmp[3] = (int128)rop[0] * polys_P[0][3] + (int128)rop[1] * polys_P[0][2] + (int128)rop[2] * polys_P[0][1] + (int128)rop[3] * polys_P[0][0] - 2*c10 + c11;
	tmp[4] = (int128)rop[0] * polys_P[0][4] + (int128)rop[1] * polys_P[0][3] + (int128)rop[2] * polys_P[0][2] + (int128)rop[3] * polys_P[0][1] + (int128)rop[4] * polys_P[0][0] - 2*c11 + c12;
	tmp[5] = (int128)rop[0] * polys_P[0][5] + (int128)rop[1] * polys_P[0][4] + (int128)rop[2] * polys_P[0][3] + (int128)rop[3] * polys_P[0][2] + (int128)rop[4] * polys_P[0][1] + (int128)rop[5] * polys_P[0][0] - 2*c12 + c13;
	tmp[6] = (int128)rop[0] * polys_P[0][6] + (int128)rop[1] * polys_P[0][5] + (int128)rop[2] * polys_P[0][4] + (int128)rop[3] * polys_P[0][3] + (int128)rop[4] * polys_P[0][2] + (int128)rop[5] * polys_P[0][1] + (int128)rop[6] * polys_P[0][0] - 2*c13 + c14;
	tmp[7] = (int128)rop[0] * polys_P[0][7] + (int128)rop[1] * polys_P[0][6] + (int128)rop[2] * polys_P[0][5] + (int128)rop[3] * polys_P[0][4] + (int128)rop[4] * polys_P[0][3] + (int128)rop[5] * polys_P[0][2] + (int128)rop[6] * polys_P[0][1] + (int128)rop[7] * polys_P[0][0] - 2*c14;

	internal_reduction(rop, tmp);
}


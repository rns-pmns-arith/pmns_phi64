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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmpQ[5];
	int128 tmpZero[5];

	//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)
	tmpQ[0] = (uint64_t)op[0] * 4619465585820490359UL + (uint64_t)op[1] * 10663724694026160790UL + (uint64_t)op[2] * 703491987568749296UL + (uint64_t)op[3] * 9932607448416294812UL + (uint64_t)op[4] * 12568926044450165254UL;
	tmpQ[1] = (uint64_t)op[0] * 6284463022225082627UL + (uint64_t)op[1] * 4619465585820490359UL + (uint64_t)op[2] * 10663724694026160790UL + (uint64_t)op[3] * 703491987568749296UL + (uint64_t)op[4] * 9932607448416294812UL;
	tmpQ[2] = (uint64_t)op[0] * 14189675761062923214UL + (uint64_t)op[1] * 6284463022225082627UL + (uint64_t)op[2] * 4619465585820490359UL + (uint64_t)op[3] * 10663724694026160790UL + (uint64_t)op[4] * 703491987568749296UL;
	tmpQ[3] = (uint64_t)op[0] * 9575118030639150456UL + (uint64_t)op[1] * 14189675761062923214UL + (uint64_t)op[2] * 6284463022225082627UL + (uint64_t)op[3] * 4619465585820490359UL + (uint64_t)op[4] * 10663724694026160790UL;
	tmpQ[4] = (uint64_t)op[0] * 14555234383867856203UL + (uint64_t)op[1] * 9575118030639150456UL + (uint64_t)op[2] * 14189675761062923214UL + (uint64_t)op[3] * 6284463022225082627UL + (uint64_t)op[4] * 4619465585820490359UL;

	//~ computation of : tmp_q*red_int_coeff mod(E)
	tmpZero[0] = (int128)tmpQ[0] * 377098328410739L + (int128)tmpQ[1] * 3230851069417532L + (int128)tmpQ[2] * 2084744545767434L + (int128)tmpQ[3] * 2767020584642314L + (int128)tmpQ[4] * 2795886552619278L;
	tmpZero[1] = (int128)tmpQ[0] * 1397943276309639L + (int128)tmpQ[1] * 377098328410739L + (int128)tmpQ[2] * 3230851069417532L + (int128)tmpQ[3] * 2084744545767434L + (int128)tmpQ[4] * 2767020584642314L;
	tmpZero[2] = (int128)tmpQ[0] * 1383510292321157L + (int128)tmpQ[1] * 1397943276309639L + (int128)tmpQ[2] * 377098328410739L + (int128)tmpQ[3] * 3230851069417532L + (int128)tmpQ[4] * 2084744545767434L;
	tmpZero[3] = (int128)tmpQ[0] * 1042372272883717L + (int128)tmpQ[1] * 1383510292321157L + (int128)tmpQ[2] * 1397943276309639L + (int128)tmpQ[3] * 377098328410739L + (int128)tmpQ[4] * 3230851069417532L;
	tmpZero[4] = (int128)tmpQ[0] * 1615425534708766L + (int128)tmpQ[1] * 1042372272883717L + (int128)tmpQ[2] * 1383510292321157L + (int128)tmpQ[3] * 1397943276309639L + (int128)tmpQ[4] * 377098328410739L;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmpZero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmpZero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmpZero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmpZero[4]) >> WORD_SIZE;
}

void exact_coeffs_reduction(int64_t *rop, int64_t *op){

	int i;
	int128 tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp[i] = (int128) op[i];

	internal_reduction(rop, tmp);

	tmp[0] = (int128)rop[0] * polys_P[0][0] + (((int128)rop[1] * polys_P[0][4] + (int128)rop[2] * polys_P[0][3] + (int128)rop[3] * polys_P[0][2] + (int128)rop[4] * polys_P[0][1]) << 1);
	tmp[1] = (int128)rop[0] * polys_P[0][1] + (int128)rop[1] * polys_P[0][0] + (((int128)rop[2] * polys_P[0][4] + (int128)rop[3] * polys_P[0][3] + (int128)rop[4] * polys_P[0][2]) << 1);
	tmp[2] = (int128)rop[0] * polys_P[0][2] + (int128)rop[1] * polys_P[0][1] + (int128)rop[2] * polys_P[0][0] + (((int128)rop[3] * polys_P[0][4] + (int128)rop[4] * polys_P[0][3]) << 1);
	tmp[3] = (int128)rop[0] * polys_P[0][3] + (int128)rop[1] * polys_P[0][2] + (int128)rop[2] * polys_P[0][1] + (int128)rop[3] * polys_P[0][0] + (((int128)rop[4] * polys_P[0][4]) << 1);
	tmp[4] = (int128)rop[0] * polys_P[0][4] + (int128)rop[1] * polys_P[0][3] + (int128)rop[2] * polys_P[0][2] + (int128)rop[3] * polys_P[0][1] + (int128)rop[4] * polys_P[0][0];

	internal_reduction(rop, tmp);
}


#include "useful_functs.h"
#include "add_mult_poly.c"


//~ Assumes allocation already done for 'rop'.
//~ IMPORTANT : convertion to montgomery domain will be done here
void from_int_to_amns(int64_t *rop, mpz_t op){
	int i;
	mpz_t tmp;
	int128 tmp_poly[NB_COEFF];
	int128 tmp_sum[NB_COEFF];

	mpz_init_set(tmp, op);

	for(i=0; i<NB_COEFF; i++){
		rop[i] = 0;
		tmp_sum[i] = 0;
	}

	if(tmp->_mp_size == 0)
		return;

	i = 0;
	while(tmp->_mp_size && (i < NB_COEFF)){
		scalar_mult_lpoly(tmp_poly, polys_P[i++], (tmp->_mp_d[0]) & CONV_MASK);
		add_lpoly(tmp_sum, tmp_sum, tmp_poly);

		mpz_tdiv_q_2exp (tmp, tmp, RHO_LOG2);
	}

	internal_reduction(rop, tmp_sum);

	mpz_clear(tmp);
}

//~ Assumes "rop" already initialized.
//~ IMPORTANT : convertion from montgomery domain will be done here.
void from_amns_to_int(mpz_t rop, int64_t *op){
	int i;
	mpz_t tmp_sum;
	int64_t tmp_conv[NB_COEFF];

	mpz_init(tmp_sum);

	//~ convertion out of mont domain
	from_mont_domain(tmp_conv, op);

	mpz_set_si(rop, tmp_conv[0]);
	for(i=0; i<POLY_DEG; i++){
		mpz_mul_si(tmp_sum, gama_pow[i], tmp_conv[i+1]);
		mpz_add(rop, rop, tmp_sum);
	}
	mpz_mod (rop, rop, modul_p);

	mpz_clear(tmp_sum);
}

//~ computes : op/phi
void from_mont_domain(int64_t *rop, int64_t *op){

	int i;
	int128 tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp[i] = (int128) op[i];

	internal_reduction(rop, tmp);
}

void init_data(){

	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_init (gama_pow[i]);

	mpz_init (modul_p);


	mpz_set_str (modul_p, "15926785601034017652597688965642382288067825175930698716358344026572368065986672663012414611984553289419183160005921730973833172922104722418619405421123439799786864450243239730812998583804723116335621502019505976390375029830451440114798215201069277857994638135635659393811802681357471501026460455234964770195044337163284826543047878226247513417879926777450626391710833885101073846548044443214873631989320807409429937976557972213764087705454773110694164608530550543060968501684715090365738317974156965140615237650491169382331346715007267520044707061064695349759840758203818248732191051277899830638189202761019641789630623949450619549067324117844609086999534727777539125675140524476569271034435623443713764525024256006072071247340801939", 10);

	mpz_set_str (gama_pow[0], "7554275921064473242862655071777361516719321372748330292438712731756126492952954089682848577896760643510735755509316613357662048983909670905755465459878245426992948974679640937568286257909513053044433943318235555364609627820173790976735639539788864309141566177789888075031431758051224031467313831664947635843792753973028844815248966477281483197578977789511493090909663099123025212045608698045769437602801193089236132239734027118065388991296739707310023860942171128388091935486360257714275620702341146831531346343088240262468063550265471616320298190559140478945509789359376131900555046315647427068911134955340376600343892409871394245233233270638502688447430979131113495272978112703311393856612694010391342044702921689818305350227240467", 10);
	for(i=1; i<POLY_DEG; i++){
		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);
		mpz_mod (gama_pow[i], gama_pow[i], modul_p);
	}
}

void free_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_clear (gama_pow[i]);

	mpz_clear (modul_p);
}

//~ ----------------------------------------------------------------------------------------

//~ return a positive value if pa > pb, zero if pa = pb, or a negative value if pa < pb.
//~ Important : evaluation is done using the corresponding integers modulo 'p'.
int cmp_poly_evals(int64_t *pa, int64_t *pb){
	int rep;
	mpz_t a, b;
	mpz_inits (a, b, NULL);
	from_amns_to_int(a, pa);
	from_amns_to_int(b, pb);
	rep = mpz_cmp (a, b);
	mpz_clears (a, b, NULL);
	return rep;
}

void copy_poly(int64_t *rop, int64_t *op){
	int i;
	for(i=0; i<NB_COEFF; i++)
		rop[i] = op[i];
}

void add_lpoly(int128 *rop, int128 *pa, int128 *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_lpoly(int128 *rop, int64_t *op, uint64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = (int128)op[j] * scalar;
}

void print_element(int64_t *poly){
	int i;
	printf("[");
	for (i=0; i<POLY_DEG; i++)
		printf("%2ld, ", poly[i]);
	printf("%2ld]", poly[POLY_DEG]);
}


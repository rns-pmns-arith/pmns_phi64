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


	mpz_set_str (modul_p, "3488943444255984364423638560731936428715107471723434054443952878868966919559379757428364349975781134968509334011980095527219462812461399771310244963248735999426700992090962902533883577229606875094117258741240688354631659733595050827192335631065187160955376529971403733009933931908648669900859281899372924380630271154780195968191696533001063230839809281455678505250494866870781190074054358018048148778014881789500760707801931505502405218013479877169969055294695492106267316438134502313779569945022254877476581257902983722757598039734354320138759717939447050291943920789957487741792306326465351936129660504659379073456789393144147455519373745409507050627214056259245664156620348835904974622867848099392592473209749829084065122386071954467135944925632782097627505111743756462706029263492488732968950814159737211755062471407851531914473431428918003810660327243409693357152864209017449181337883131283254172508954989673076520796121761852136826119438420446446568445641649839317681986361", 10);

	mpz_set_str (gama_pow[0], "3210860690608991668675022854880194121595353138615824408468484282421497457462099795001258984331744731311117717805064898994543358334558274648768106208835334061851290585561927817924748122708461688185451249592558395452014167999957525380317214183423439121634583466113067228928344387086726514798677962804860859894753126391223554099055912501508990765441287398996758837937592217400632847996120545598908490878465303606971167906341418257925998631333236198966014261491537301435114413432703312009486690283420211651586313979962451317327877023090956664255666453719824271036758823849235456815853152808841999482466324090871519376447929456619708754798272193129930963475572777034380074214084485170195774616103772184082586934359228894064375996301862208426639079799440492807611615313593390649331764968641239345951670178055599163940764320583865159732659434393607641184977952666447193554319243537336499473249068410091475723935813443605290932980824789620628393107421150358205318874833004226585836990100", 10);
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


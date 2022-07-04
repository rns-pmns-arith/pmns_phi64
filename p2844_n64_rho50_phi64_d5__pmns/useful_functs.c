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


	mpz_set_str (modul_p, "7822126819527056218802607089799723688277078923429881893774431984780683979050944261621918589336396407952502557270927247590973487015040028923590853476103749692938431661010438523369196953471629677631432885259364894860912554815667925599435059109549093446109655405477787683166270363660638570848643117803512279036561922794756949776849740051447341326277035272537590858582927158208347236769078631867973411626657148525743870686053389163912608774915888111233282369715873434752039635716449640739656593831649265944256993481497293342063978844373452286067972893588920526862959538504448909415504347870095870767902454984722647912788177994367632615555380874928105521103488524838463769790894781508916382678310574318484278971645145998208894686678157739959198290600286793663408917427697455276762106872308434241571910556553001625981763454185596072831183063410527103715208043159", 10);

	mpz_set_str (gama_pow[0], "3555783638016291246895304712065878081624418371634637869811792584019762624501900703041379425744921408031518577601860403729289747969843931096392530973811640043079795929864157067952712815076245539245847562620222342827558409193519577785222500194695243038227032966135822191262555423901286978615164126583617200598896225569956937471531839476472030586504045583886892057345400798132904443345099191994048467653334663099591027732289999583546994935589781539347535961221255598256132113741234611763326915843735919814725558576875052726982293049219697378029228248404251403635170720395666028461450353457142127511070167109805112934159311858905718323436104927105291895960463415356371878888105171367725503459114486462974780026881835782197975806917298881749999036483695787245479408570345491770340259219180454423828584004909799835976412047992480433403053802932074885525805595200", 10);
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


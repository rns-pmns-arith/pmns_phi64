#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <gmp.h>
#include <openssl/bn.h>

#include "ccount.h"
#include "gmp_stuff.c"
#include "useful_functs.h"

extern mpz_t modul_p;

#define NSAMPLES 50
#define NTEST 1000

unsigned long long int START, STOP, START1,STOP1;


int main(int argc, char* argv[]){

	uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;

	int64_t pa[NB_COEFF];
	int64_t pb[NB_COEFF];
	int64_t pc[NB_COEFF];
	
	srand(time(NULL));
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	int nb_limbs;
	mpz_t A, B, C, E, F, R2, R3, G, H;
	mpz_inits (A, B, C, E, F, R2, R3, G, H, NULL);
	
	mp_limb_t mip0;
	
	const mp_limb_t *p_limbs;
	mp_limb_t *r2_limbs, *scratch_limbs, *mip_limbs;
	mp_limb_t *a_limbs, *b_limbs, *am1_limbs, *am2_limbs, *bm1_limbs, *bm2_limbs;
	
	init_data();
	
	mpz_urandomm(A, r, modul_p);
	mpz_urandomm(B, r, modul_p);
	mpz_set(E, A);

	
	nb_limbs = mpz_size (modul_p);
	
	p_limbs = mpz_limbs_read (modul_p);   
	
	binvert_limb (mip0, p_limbs[0]);
	mip0 = -mip0;
	
	mpz_setbit (R2, 2*nb_limbs*8*sizeof(mp_limb_t)); 
	mpz_mod(R2, R2, modul_p);
	r2_limbs = mpz_limbs_modify (R2, nb_limbs); 
	
	a_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	am2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm1_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	bm2_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	mip_limbs = (mp_limb_t*) calloc (nb_limbs, sizeof(mp_limb_t));
	scratch_limbs = (mp_limb_t*) calloc (2*nb_limbs, sizeof(mp_limb_t));

	b_limbs = mpz_limbs_modify (B, nb_limbs);
	copy_limbs(a_limbs, A, nb_limbs);
	copy_limbs(am1_limbs, A, nb_limbs);
	copy_limbs(am2_limbs, A, nb_limbs);
	copy_limbs(bm1_limbs, B, nb_limbs);
	copy_limbs(bm2_limbs, B, nb_limbs);
	
	mpn_binvert (mip_limbs, p_limbs, nb_limbs, scratch_limbs);
	

	/***********************************************/
	printf("\nGMP mult_mod :\n----------\n");
	gmp_printf("p       : %Zd\n\n", modul_p);
	gmp_printf("A       : %Zd\n", E);
	gmp_printf("B       : %Zd\n\n", B);

	//~ conversion to Montgomery domain (block mont)
	mpn_mont_mul_red_1(am1_limbs, am1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	mpn_mont_mul_red_1(bm1_limbs, bm1_limbs, r2_limbs, p_limbs, mip0, nb_limbs);
	
	//~ Montgomery modular multiplication (block mont)
	mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(int i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am1_limbs[i];
	mpn_redc_1 (am1_limbs, scratch_limbs, p_limbs, nb_limbs, mip0);
	from_limbs_to_mpz_t(G, am1_limbs, nb_limbs);
	gmp_printf("r_mbgmp : %Zd\n", G);
	
	//~ conversion to Montgomery domain (classic mont)
	mpn_mont_mul_red_n(am2_limbs, am2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	mpn_mont_mul_red_n(bm2_limbs, bm2_limbs, r2_limbs, p_limbs, mip_limbs, nb_limbs);
	
	//~ Montgomery modular multiplication (block mont)
	mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
	
	clean_limbs(scratch_limbs, 2*nb_limbs);
	for(int i=0; i<nb_limbs; i++)
		scratch_limbs[i] = am2_limbs[i];
	mpn_redc_n (am2_limbs, scratch_limbs, p_limbs, nb_limbs, mip_limbs);
	from_limbs_to_mpz_t(H, am2_limbs, nb_limbs);
	gmp_printf("r_mcgmp : %Zd\n", H);
	
	
	//~ Montgomery low level modular multiplication
	mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
	from_limbs_to_mpz_t(F, a_limbs, nb_limbs);
	gmp_printf("r_lgmp  : %Zd\n", F);
	
#ifdef TEST
	goto chrono;
#endif
	
	mpz_mul (E, A, B);
	mpz_mod (E, E, modul_p);
	gmp_printf("r_gmp   : %Zd\n", E);


	/***********************************************/
	printf("\nmult_mod_poly n=%d, phi=2^64 :\n------------------------------------\n",NB_COEFF);
	
	from_int_to_amns(pa, A);
	from_int_to_amns(pb, B);
	printf("pA	:");print_element(pa); printf("\n");
	printf("pB	:");print_element(pb); printf("\n");
	//goto fin;

	mult_mod_poly(pc, pa, pb);
	printf("pC	:");print_element(pc); printf("\n");
	
	from_amns_to_int(C, pc);
	gmp_printf("r_pmns  : %Zd\n", C);
	
	printf("\n\n");
	
	//~ printf("\nAll results :\n-------------\n");
	
	//~ gmp_printf("p       : %Zd\n\n", modul_p);
	//~ gmp_printf("A       : %Zd\n", A);
	//~ gmp_printf("B       : %Zd\n\n", B);
	//~ gmp_printf("r_pmns  : %Zd\n", C);
	
	//~ printf("\nGMP mult_mod :\n----------\n");
	//~ gmp_printf("r_mbgmp             : %Zd\n", G);
	//~ gmp_printf("r_mcgmp             : %Zd\n", H);
	//~ gmp_printf("r_lgmp              : %Zd\n", F);
	//~ gmp_printf("r_gmp               : %Zd\n\n\n", E);
	//~ //goto fin;//*/

//~ chrono:
	
	/********************
	//Timings  !!!!!!!!!!!!!
	*********************/
	
	#define K_T 10
	
	
	printf("\t  /*********************/\n");
	printf("\t / Timings !!!!!!!!!!!!/\n");
	printf("\t/*********************/\n\n");

 
	unsigned long long int timer=0, timer4=0, timer5=0, timer6=0, timer7=0;

	
	// timer mult_mod_poly        :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini1 = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mult_mod_poly(pc, pa, pb);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			STAMP(START1)
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);//*/
			mult_mod_poly(pc, pa, pb);
			STAMP(STOP1)

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
		}
		
		timer += mini1;
	}

	printf("timer mult_mod_poly                : %llu\n",timer/(K_T*NSAMPLES));


	
	// timer gmp mpz                      :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer4 += mini;
	}
	
	printf("timer gmp mpz                      : %llu\n",timer4/(K_T*NSAMPLES));

	
	// timer gmp low level                :
		
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;

		for(int i=0;i<NTEST;i++)
		{
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer5 += mini;
	}
	
	printf("timer gmp low level                : %llu\n",timer5/(K_T*NSAMPLES));


	// timer gmp Montgomery               :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer6 += mini;
	}
	
	printf("timer gmp Montgomery               : %llu\n",timer6/(K_T*NSAMPLES));


	// timer gmp Montgomery bloc          :

	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer7 += mini;
	}
	printf("timer gmp Montgomery bloc          : %llu\n",timer7/(K_T*NSAMPLES));
	
/*********************************************************

	Instruction counting with rdpmc

*********************************************************/

	printf("\n\t  /**************************/\n");
	printf("\t / Instructions !!!!!!!!!!!!/\n");
	printf("\t/**************************/\n\n");
	
	
	timer=0, timer4=0, timer5=0, timer6=0, timer7=0;

	
	// #Instructions mult_mod_poly_n6_pmns        :
	
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini1 = (uint64_t)-1L;

		
		for(int i=0;i<NTEST;i++)
		{
			mult_mod_poly(pc, pa, pb);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START1 = rdpmc_instructions();
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);
			mult_mod_poly(pc, pa, pb);//*/
			mult_mod_poly(pc, pa, pb);
			STOP1 = rdpmc_instructions();

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
			
		}
		
		timer += mini1;
	}
		
	printf("#Instructions mult_mod_poly                : %llu\n",timer/(K_T*NSAMPLES));

	
	// #Instructions gmp mpz                      :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			mpz_mul (A, A, B);
			mpz_mod (A, A, modul_p);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer4 += mini;
	}
	
	printf("#Instructions gmp mpz                      : %llu\n",timer4/(K_T*NSAMPLES));

	
	// #Instructions gmp low level                :
		
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;

		for(int i=0;i<NTEST;i++)
		{
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			mpn_mod_mult(a_limbs, a_limbs, b_limbs, p_limbs, nb_limbs);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer5 += mini;
	}
	
	printf("#Instructions gmp low level                : %llu\n",timer5/(K_T*NSAMPLES));


	// #Instructions gmp Montgomery               :
	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			mpn_mont_mul_red_n(am2_limbs, am2_limbs, bm2_limbs, p_limbs, mip_limbs, nb_limbs);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer6 += mini;
	}
	
	printf("#Instructions gmp Montgomery               : %llu\n",timer6/(K_T*NSAMPLES));


	// #Instructions gmp Montgomery bloc          :

	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L;
		
		for(int i=0;i<NTEST;i++)
		{
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			START = rdpmc_instructions();
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			mpn_mont_mul_red_1(am1_limbs, am1_limbs, bm1_limbs, p_limbs, mip0, nb_limbs);
			STOP = rdpmc_instructions();
			
			if(mini>STOP-START) mini = STOP-START;

			
		}

		timer7 += mini;
	}
	
	printf("#Instructions gmp Montgomery bloc          : %llu\n",timer7/(K_T*NSAMPLES));


//~ fin:
	free(a_limbs);
	free(am1_limbs);
	free(am2_limbs);
	free(bm1_limbs);
	free(bm2_limbs);
	free(mip_limbs);
	free_data();

	printf("\n");
	printf("Au revoir et merci !\n\n");

}







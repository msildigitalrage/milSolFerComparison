/*
 ============================================================================
 Name        : 1.c
 Author      : 
 Version     :
 Copyright   : thomas.zafeiropoulos
 Description : cryptography Fermat, Sollovay-Strassen, Miller-Rabin C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <unistd.h>
gmp_randstate_t stat;
int main(void) {
	long sd = 0;
	int t = 20;
	int cs=0; int ds=0;//for 2o algorithmo
	int s,j_rab;//miller
	int result = 0; //miller
	mpz_t psi, d, n_minus_one;//miller
mpz_t n,three,two,a,r,b,k,one,ss,jac_ch;
mpz_t seed;
gmp_randinit(stat,GMP_RAND_ALG_LC,120);
mpz_init(n);//iniatialize
mpz_init(three);
mpz_init(a);
mpz_init(two);
mpz_init(one);
mpz_init(r);
mpz_init(k);
mpz_init(b);
mpz_init(ss);
mpz_init(jac_ch);
mpz_init(seed);
mpz_init(psi);//for miller-rabin
mpz_set_str(three, "3", 10);
mpz_set_str(two, "2", 10);
mpz_set_str(one, "1", 10);
srand((unsigned)getpid()); //initialize stat
sd = rand();
mpz_set_ui(seed,sd);
gmp_randseed(stat,seed);
int num=1000;//1000 numbers to create
int counter=0; int i = 0; int j;
int Fermat_primes=0; int SS_primes=0; int Miller_primes=0;
int Fermat_c=0; int SS_c=0; int Miller_c=0;
for(j=0; j<num; j++){

do{//RANDOM NUMBER
	mpz_urandomb(n,stat,50);//create a random number
}while((mpz_even_p(n))&& (n<=three));//checking n to be >=3 && n be odd
//Fermat Calculation starts
	for (i = 0; i < t; ++i) {
		do{
		mpz_urandomb(a,stat,50);//create random number a
		}while(!((a<=(n-two)) && (a>=two)));//checking a must be (2<=a<=n-2)
mpz_sub(k,n,one);//k=n-1
mpz_powm(r,a,k,n);//r=a^k mod n
if(mpz_cmp(r,one)!=0){ i=t-1; Fermat_c++;}//returns 0 if r=one | for r !=1 then mynumber is composite
}if(mpz_cmp(r,one)==0){Fermat_primes++;}    // for r = 1 mynumber is prime
//Sollovay-Strassen Compute
for (i = 0; i < t; ++i) {
		do{
		mpz_urandomb(a,stat,50);//create random number
		}while(!((a<=(n-two)) && (a>=two)));//checking a must be (2<=a<=n-2)
mpz_sub(b,n,one);//n-1=b
mpz_div(k,b,two); // k= n-1/2
mpz_powm(r,a,k,n);//r=a^k mod n
if((mpz_cmp(r,one)!=0) && (mpz_cmp(r,b)!=0)){ i=t-1; goto endme;
SS_c++;
} //if(r!=1 && !=n-1) then composite
else{
int jac= mpz_jacobi(a,n); //an de vrikes kati compute jacobi
if (jac==-1){ //make jac from -1 to 0 and add to r one
	mpz_add(r,r,one);
	jac=0;
}
	int bob = mpz_congruent_ui_p(r,jac,n);  // ine isotimo to r me to jacobi (jac)?
if(bob!=0){
	SS_c++;
	goto endme;
}
}
}
SS_primes++; goto endme;//an de vrethike composite tote ine prime

goto endme;
endme:
SS_primes=num-SS_c;
//Miller-Rabin Calculation
mpz_init(n_minus_one); //initialize
mpz_sub_ui(n_minus_one, n, 1);//n_minus_one = n-1
	s = 0;
	mpz_init_set(d, n_minus_one);
	while (mpz_even_p(d)) {// gia oso ine artios
		mpz_fdiv_q_2exp(d, d, 1); //shift right
		s++;//inc s
	}
	for (i = 0; i < t; ++i) {
		do{
		mpz_urandomb(a,stat,50);//create r%Zdandom number
		}while(!((a<=(n-two)) && (a>=two)));//checking a must be (2<=a<=n-2)
		 mpz_powm(psi,a,d,n);
		 if(mpz_cmp(psi,one)!=0 && mpz_cmp(psi,n_minus_one)){
			 j_rab=1;
			 while(j_rab<s &&  mpz_cmp(psi,n_minus_one)){
				 mpz_mul(psi,psi,psi); // y^2
				 mpz_mod(psi,psi,n); //psi= psi^2 mod n
				 	 if(mpz_cmp(psi,one)==0){//if y=1
				 		Miller_c++;  result = 1; goto exit_miller;
				 	 }
				 	 j_rab++;
			 }
		 	 if(mpz_cmp(psi,n_minus_one)!=0){//if y=n-1
		 		 Miller_c++; result = 1; goto exit_miller;
		 	 }
		 }//end external if
}//end for
	if(result!=1){ Miller_primes++;}
	exit_miller: result = 0;
counter++;
}
printf("\n%d Random numbers\n",counter);
printf("\n\tFermat: Primes = %d\tComposite = %d",Fermat_primes,Fermat_c);
printf("\n\tSollovay-Strassen: Primes = %d\tComposite = %d",SS_c,SS_primes);
printf("\n\tMiller-Rabin: Primes = %d\tComposite = %d",Miller_primes,Miller_c);
printf("\n");
mpz_clear(n);//clear
mpz_clear(three);
mpz_clear(a);
mpz_clear(two);
mpz_clear(r);
mpz_clear(seed);
mpz_clear(k);
mpz_clear(b);
mpz_clear(one);
mpz_clear(ss);
mpz_clear(jac_ch);
mpz_clear(d);
mpz_clear(n_minus_one);
mpz_clear(psi);
	return 0;
}

/* Calculate 2D, mq fast mas lineshapes, for any 1/2 integer spin and shear factor, Cq, iso cs distributions, multiple sites, delta m = 1,2,3 or 5

   OFFERED UNDER GNU GPL LICENSE. NO GUARANTEES OF FITNESS FOR USE, NO WARRANTIES IMPLIED OR GIVEN.

   likely compile command, x86 LINUX:

   gcc -o mqmas_opt mqmas_opt.c -lgslcblas -lm -lgsl -fopenmp

   WJB 12/07, 01/08, 06/08, 12/13 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#define M_PI 3.1415

//#include <gsl/gsl_qrng.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>




/*variable decs;
  loop counters */
int h, i, j, k, l;

/*  angle space */
int m, tt;
/* n samples cs/cq space */
int n;
/* dimensions f1/f2 */
int p, q;
/* Fib series*/

int M, N;

/*input parameters */
double *a_ptr, lm, dw, du, w, u, kk, sp, sh;
double mq, dt, sn, lg, s, r, rn, sn2;
int proc_n[3]={0,0,0};
double Q;

int step, Lstart, Lfinish, Sstart, Sfinish;
double f1A_der, f2A_der;
/*RNG */
double *rn_ptr, *gam_ptr, *rn_ptr2;

/* angles, functions etc     */

double a, b, *fab_ptr, *fib_ptr, *fabDer_ptr;
/* axes */
double *f1_ptr, *f2_ptr, *Gf2, *Lf2, *Gf22, *Gf21, *Lf22, *Lf21;
double *f1der_ptr, *f2der_ptr;

double f2m_der, f1m_der;

double cs_der, cs_stddev_der, cq_der, cq_stddev_der;

/* NMR parameters */
double *cq_ptr, *cs_ptr;

double clb[4], t, u;

double Omega;

/* output */
double *fn_ptr, *s_ptr;

/* temp values */
double cs, cq, f2m, f1m, nuq, nm, alp, ww;
double f2Gb, f1Gb, f2Lb, f1Lb, f2A, f1A, fA, NrM, LNrM;
double LNrM1, LNrM2, NrM1, NrM2;

double f1m_der, f2m_der; 

double *f1m_der_ptr, *f2m_der_ptr;

int inpFile=0;

double rho;

/*strings */
char str_2[]="gamma";
char str_1[]="Gaussian";
char str_3[]="log normal";

/*function decs, file i/o*/

int FileOfDoublesSize(FILE *fin);
void ErrorMessage(char *str);
double *BinRead(FILE *fin, int n);
void DispVector(double *ptr, int n);
void FileWrite(FILE *fout, double *ptr, int t);

/*math */
double *Fib(int r);
double Gam(int n);
int Mod(int g, int gg);

/*RNG's */
//double *GaussianRNG(double *ptr, double *ptr2, double *ptr3, int n, int m);
//double *GammaRNG(double *ptr, double *ptr2, double *ptr3, int n, int m);
//double *LogNormalRNG(double *ptr, double *ptr2, double *ptr3, int n, int m);


int main(int argc, char** argv)

{
	/* begin main */

	/*input parameters for operation modes */

	if (argc>1){

		for (i = 1; i < argc; i++) {

			/* These inputs may be used to specify running modes;

			   proc_n[0] 0==default, no output, 1==display inputs, 2==display site & samples with output, 3==deriv wrt lambda_2, 4==deriv wrt lambda_1,
			   5== deriv wrt lg, 6==deriv wrt cs mean, 7==deriv wrt cs std deviation, 8==deriv wrt cq mean, 9==deriv wrt cq st dev, 10==deriv rho (unimplemented
			   as 0f 0608), 11==deriv wrt eta, 12== deriv wrt Ampl

			   proc_n[1] is the processor number (default 0, single machine)
			   proc_n[2] is the total number of processors (ditto)
			   */


			proc_n[i-1]=atoi(argv[i]);

		}

		if ((proc_n[0]==1) && (proc_n[1]!=0))
		{
			printf("proc no./total procs= %i / %i\n",proc_n[1],proc_n[2]);
		}

	}

	/* end get input parameters */



	FILE *fptr1;  
	char filename1[]="input.bin";

	/*

	   FILE *fptr2, *fptr3;
	   char filename2[]="output.bin";
	   char filename3[]="samples.bin";

*/

	int count=0;
	int size;
	if ((fptr1 = fopen(filename1, "r+b")) == NULL){
		// ErrorMessage(filename1);
		float value=0.0;
		int foo=1;
		while (foo != -1){


			foo=scanf("%f ", &value);
			if (foo != -1){
				//printf("%f %i\n",value,foo);
				if (count==0) {

					size = (int) value;
					a_ptr = malloc(size * sizeof(double));
				} else{
					a_ptr[count-1]=(double) value;
				}
				count++;
			}
		}
		m=size;
	} else {

		inpFile=1;		 

		m=FileOfDoublesSize(fptr1);
		a_ptr = BinRead(fptr1,m);
	} 
	s=floor((m-17.0)/8.0+0.01);
	/* Initialize */
	sp=*(a_ptr); /*spin */
	mq=(int) *(a_ptr+1); /* transition; 1,2 (satellites) 3,5 (symmetric mq)*/
	sh=*(a_ptr+2); /*shear factor */
	lm=*(a_ptr+3); /*larmor */
	m=(int) *(a_ptr+4); /*powder steps, order Fib series*/
	n=(int) *(a_ptr+5); /*samples of cs, cq distributions */
	p=(int) *(a_ptr+6); /*# F2 frequency points */
	q=(int) *(a_ptr+7); /*# F1 frequency points */
	dw=*(a_ptr+8); /*step F2*/
	du=*(a_ptr+9); /*step F1*/
	w=*(a_ptr+10); /*start F2 */
	u=*(a_ptr+11); /*start F1 */
	dt=(int) *(a_ptr+12); /*distribution type; 1==gaussian, 2==gamma, 3==lognormal */
	sn=*(a_ptr+13); /*empty variable, unimplemented*/
	rn=(int) *(a_ptr+14); /*1==quasi random, else pseudo random */
	lg=0.5;   //*(a_ptr+15); Lorentz/Gauss broadening ratio 
	sn2=*(a_ptr+16); /* empty variable, unimplemented */



	/* remainder of input should read, for each consecutive chemical site

	   iso chemical shift mu (gaussian,log normal) OR k (gamma) *(a_ptr+17+j*8)  
	   iso chemical shift sigma (gaussian,log normal) OR theta (gamma) *(a_ptr+18+j*8)

	   quadrupole coupling 	... *(a_ptr+19+j*8)
	   "			... *(a_ptr+20+j*8)

	   asymmetry parameter	    *(a_ptr+21+j*8)

	   broadening F2			*(a_ptr+22+j*8)
	   broadening F1			*(a_ptr+23+j*8)

	   Amplitude			*(a_ptr+24+j*8)

*/



	/* Safety features */
	if(m > 30){
		m=30;
	}


	if(p > 1024){
		p=1024;
	}

	if(q > 1024){
		q=1024;
	}

	if(n > 1024){
		n=1024;
	}

	if(lg > 1.0){
		lg=1.0;
	}


	/*default loop limits*/

	//Lstart=0;
	//Lfinish=n;

	Sstart=0;
	Sfinish=s;

	/* set aside some memory */



	s_ptr=malloc(4*s*n*sizeof(double));

	f2_ptr=malloc(p*sizeof(double));
	f1_ptr=malloc(q*sizeof(double));

	f2der_ptr=malloc(p*q*sizeof(double));
	f1der_ptr=malloc(p*q*sizeof(double));

	Gf2=malloc(p*sizeof(double));
	Lf2=malloc(p*sizeof(double));


	Gf21=malloc(p*sizeof(double));
	Lf21=malloc(p*sizeof(double));


	Gf22=malloc(p*sizeof(double));
	Lf22=malloc(p*sizeof(double));


	*f2_ptr=w;
	*f1_ptr=u;

	for (i=1;i<p;i++){
		*(f2_ptr+i)=dw + *(f2_ptr+i-1);
	}



	for (i=1;i<q;i++){
		*(f1_ptr+i)=du + *(f1_ptr+i-1);
	}


	rho=0;

	/* display inputs */

	if (proc_n[0]==1){
		printf("No. of sites: %f\n",s);

		printf("points F2= %i\n",p);
		printf("points F1= %i\n",q);
		printf("dist. samples= %i\n",n);
		printf("larmor freq (Hz)= %e\n",lm);
		printf("step size F2 (Hz)= %f\n",dw);
		printf("step size F1 (Hz)= %f\n",du);
		printf("start freq F2 (Hz)= %f\n",w);
		printf("start freq F1 (Hz)= %f\n",u);
		printf("spin= %f\n",sp);
		printf("mq transition= %f\n",mq);
		printf("shear= %f\n",sh);

		if (dt==3) {printf("dist type= %s\n",str_3);}


		else if (dt==2) {printf("dist type= %s\n",str_2);}
		else {printf("dist type= %s\n",str_1);}

		printf("L/G ratio= %f\n",lg);
	}

	/*carve up sites iff proc[2]!=0 & n==1 & s>proc_n[2]*/

	if ((proc_n[2]!=0)&&(n==1)&&(s>proc_n[2])){


		step=floor(s/proc_n[2]);

		Sstart=proc_n[1]*step;
		Sfinish=(proc_n[1]+1)*step;

		if (proc_n[0]==1){
			printf("sites section= %i : %i\n",Sstart,Sfinish-1);
		}
	}






	/*set up sobol if need be*/




#if 0
	if (rn==1&n>1){
		gsl_qrng *qq = gsl_qrng_alloc(gsl_qrng_sobol, 2);
		rn_ptr=malloc(n*sizeof(double));
		rn_ptr2=malloc(n*sizeof(double));

		for (i = 0; i < n; i++)
		{ double vv[2];
			gsl_qrng_get(qq, vv);
			*(rn_ptr+i)=vv[0];
			*(rn_ptr2+i)=vv[1];

		}
		gsl_qrng_free(qq);


		if (proc_n[0]==1){
			printf("RNG type: quasi\n");
		}
		if (proc_n[2]!=0){
			/* carve up numbers for HTC/condor */

			step=floor(n/proc_n[2]);

			Lstart=proc_n[1]*step;
			Lfinish=(proc_n[1]+1)*step;

			if (proc_n[0]==1){
				printf("RNG sequence section= %i : %i\n",Lstart,Lfinish-1);
			}
		}


	}/* end QRN loop */

	/* or set up pseudo RN's */
	else if (rn!=1&n>1) {
		const  gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		rn_ptr=malloc(n*sizeof(double));
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		for (i = 0; i < n; i++) 
		{
			*(rn_ptr+i) = gsl_rng_uniform (r);

		}

		gsl_rng_free (r);

		if (proc_n[0]==1){
			printf("RNG type: pseudo\n");
		}



	}
#endif
	/* end PRN loop */
	/* DispVector(rn_ptr,n); */


	/* set up freq. axes */
	/* freq parameters */



	t=(double)(mq==5)*5/2+(double)(mq==3)*3/2+(double)(mq==2)*3/2+(double)(mq==1)*3/2;
	u=(double)(mq==5)*-5/2+(double)(mq==3)*-3/2+(double)(mq==2)*-1/2+(double)(mq==1)*1/2;

	clb[0]=-(sp*(sp+1.0)-3.0/4.0);
	clb[1]=-18.0*sp*(sp+1.0)+34.0/4.0+5.0;
	clb[2]=(t-u)*(sp*(sp+1.0)-3.0*(t*t+t*u+u*u));
	clb[3]=(t-u)*(18.0*sp*(sp+1.0)-34.0*(t*t+t*u+u*u)-5.0);

	//if (proc_n[0]==1){
	//printf("freq. consts= %f, %f, %f, %f\n",clb[0],clb[1],clb[2],clb[3]);

	//}


	fn_ptr=malloc(p*q*sizeof(double));





	/* DispVector(f2_ptr,p);
	   DispVector(f1_ptr,q);
	   */



	/* initialize powder steps */

	fib_ptr=Fib(20);



	M = *(fib_ptr+m-1);
	N = *(fib_ptr+m+1);  
	fab_ptr = malloc(N*sizeof(double));

	fabDer_ptr = malloc(N*sizeof(double));

	if (proc_n[0]==1){
		printf("powder steps: %i\n",N);
	}
	/* do for each site */  


	//fudge
	*(a_ptr+21+8*j) = *(a_ptr+15);

	/* simulate lineshape */
	if (proc_n[0]<3){

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */
#if 0
			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}
#endif
			//else {
			cs_ptr=malloc(sizeof(double));
			cq_ptr=malloc(sizeof(double));

			*cs_ptr=*(a_ptr+17+j*8); 

			*cq_ptr=*(a_ptr+19+j*8);
			/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			//} 

			//for (k=Lstart;k<(Lfinish-Lstart);k++){

			cs=*(cs_ptr+k);
			//cq=abs(*(cq_ptr+k));
			cq=*(a_ptr+14) * 1e6;
			//printf("%f %f\n",cq,a_ptr[21+8*j]);
			Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 

#if 0
			if ((proc_n[0]==2)&(n>1)){

				*(s_ptr+j*n+4*k)=cq;
				*(s_ptr+j*n+4*k+1)=cs;
				*(s_ptr+j*n+4*k+2)=(double) j;
				*(s_ptr+j*n+4*k+3)=*(a_ptr+21+8*j);


			} 
#endif

			fA = pow(Omega,2) / (1.0 * lm);

			f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



			f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



			NrM = lg * *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI * n * N );



			LNrM= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (n * N ); 

			f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
			f2Lb = pow(*(a_ptr+22+8*j),2); 

			f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
			f1Lb = pow(*(a_ptr+23+8*j),2);



			for (l=0; l<N; l++){


				f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
				f1m= f1A - fA* 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;




#pragma omp parallel for private(h,i)

				for (i=0;i<p;i++){



					*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
					*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




					for (h=0;h<q;h++){



						/* printf("%f , %f\n",f2m,f1m);*/


						/* Combination Lorentzian & Gaussian */
						*(fn_ptr+i*q+h) +=  NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i); 





					}/*end f1 loop */



				}/*end f2 loop */

			}/*end powder loop */



			//}

		}/* end sites loop */


	}/* end proc_n loop*/


	//normalize
	double sum=0.0;
	for (i=0; i<p*q; i++) sum += fn_ptr[i];

	for (i=0; i<p*q; i++) printf("%i:%f ",(i+1),fn_ptr[i]/sum);
	printf("\n");



#if 0

	/* derivative wrt lambda_2 */

	else if (proc_n[0]==3){

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 



				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM1 = -lg * *(a_ptr+24+8*j) / ( *(a_ptr+23+8*j)* pow(*(a_ptr+22+8*j),2) * 2.0 * M_PI * n *N);
				NrM2 = lg * * (a_ptr+24+8*j) / ( pow(*(a_ptr+22+8*j),4)* *(a_ptr+23+8*j) * 2.0 * M_PI * n *N);


				LNrM1= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+23+8*j) / (n* N); 
				LNrM2= -2.0 * (1.0-lg) * *(a_ptr+24+8*j) * pow(*(a_ptr+22+8*j),2)* *(a_ptr+23+8*j)/ (n *N); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf21+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);

						*(Gf22+i) = pow((*(f2_ptr+i)-f2m),2) *  *(Gf21 + i); 


						*(Lf21+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));

						*(Lf22+i) = pow(*(Lf21+i),2);




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) += NrM1 *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf21+i) + NrM2 *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf22+i) + LNrM1/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf21+i) + LNrM2/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf22+i); 





						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/

	/* derivative wrt lambda_1 */

	else if (proc_n[0]==4){


		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0));  


				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM1 = -lg * *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* pow(*(a_ptr+23+8*j),2) * 2.0 * M_PI * n *N);
				NrM2 = lg * * (a_ptr+24+8*j) / ( pow(*(a_ptr+23+8*j),4)* *(a_ptr+22+8*j) * 2.0 * M_PI * n * N);


				LNrM1= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j) / (n* N); 
				LNrM2= -2.0 * (1.0-lg) * *(a_ptr+24+8*j) * pow(*(a_ptr+23+8*j),2)* *(a_ptr+22+8*j)/ (n *N); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh *f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf21+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);

						*(Lf21+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) += NrM1 *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf21+i) + NrM2 *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf21+i) *  pow((*(f1_ptr+h)-f1m),2)+LNrM1/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf21+i) + LNrM2/ pow(f1Lb + pow((*(f1_ptr+h)-f1m),2),2) * *(Lf21+i); 





						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/


	/* derivative wrt lg */
	else if (proc_n[0]==5){
		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 



				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 20) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn *(t-u) * cs + sn2 * fA * ((clb[2] / 20) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM =  *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI * n *N);



				LNrM= -1.0 * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (n* N); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) += NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i); 





						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */






	}/* end proc_n loop*/  

	/* derivative wrt cs mean */

	else if (proc_n[0]==6){

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				cs_der = - rho * (cq - *(a_ptr+19+j*8)) / (*(a_ptr+20+j*8) * *(a_ptr+18+j*8) ) +  (cs - *(a_ptr+17+j*8)) / (*(a_ptr+18+j*8) * *(a_ptr+18+j*8) );




				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 



				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM = lg * *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI * n * N);



				LNrM= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (n* N); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) += cs_der* NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + cs_der * LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i); 




						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/

	/* derivative wrt cs std dev */

	else if (proc_n[0]==7){ 

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				cs_stddev_der = - rho * (cq - *(a_ptr+19+j*8))* (cs-*(a_ptr+17+j*8)) / (*(a_ptr+20+j*8) * pow(*(a_ptr+18+j*8),2) ) +  pow((cs - *(a_ptr+17+j*8)),2) / pow(*(a_ptr+18+j*8),3) - 1 / *(a_ptr+18+j*8);



				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 


				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM = lg * *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI * n * N);



				LNrM= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (n * N); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh* f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) += cs_stddev_der * NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + cs_stddev_der * LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i); 




						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/

	/* derivative wrt cq mean */                  


	else if (proc_n[0]==8){ 

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				cq_der = - rho * (cs - *(a_ptr+17+j*8)) / (*(a_ptr+20+j*8) * *(a_ptr+18+j*8) ) +  (cq - *(a_ptr+19+j*8)) / (*(a_ptr+20+j*8) * *(a_ptr+20+j*8) );

				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 




				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn* (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM = lg * *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI * N *n);



				LNrM= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j)/ (n* N); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) +=cq_der* NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + cq_der * LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i); 


						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/


	/* derivative wrt cq std dev */

	else if (proc_n[0]==9){ 

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				cq_stddev_der=- rho * (cq - *(a_ptr+19+j*8))* (cs-*(a_ptr+17+j*8)) / (*(a_ptr+18+j*8) * pow(*(a_ptr+20+j*8),2) ) +  pow((cq - *(a_ptr+19+j*8)),2) / pow(*(a_ptr+20+j*8),3) - 1 / *(a_ptr+20+j*8);




				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 


				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM = lg * *(a_ptr+24+8*j) / ( *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI * N *n);



				LNrM= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (N *n); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;



#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */
							*(fn_ptr+i*q+h) +=cq_stddev_der* NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + cq_stddev_der * LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i); 




						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/


	/* reserved for rho */

	/* derivative wrt eta */                
	else if (proc_n[0]==11){    

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 


				/* derivative anisotropy {d F(alpha,beta)/ d eta} */


				*(fabDer_ptr+l) = (1.0/11520.0)*(-6.0* *(a_ptr+21+8*j)+60.0* cos(2*a)-70.0* *(a_ptr+21+8*j)*cos(4.0*a)+(60.0* *(a_ptr+21+8*j)-480* cos(2*a)+140* *(a_ptr+21+8*j)*cos(4.0*a) )*pow(cos(b),2)+(-70.0* *(a_ptr+21+8*j) + 420* cos(2*a)-70.0* *(a_ptr+21+8*j)*cos(4*a))*pow(cos(b),4)); 

			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				Omega=  (cq )  / (2.0 * sp *(2.0 * sp -1 ));


				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));

				f2A_der = -fA * (clb[0]/5) * *(a_ptr+21+8*j); 




				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 


				f1A_der = sn2 * fA * (clb[2]/5) * *(a_ptr+21+8*j);


				NrM = lg * *(a_ptr+24+8*j) / (pow( *(a_ptr+22+8*j),3) * *(a_ptr+23+8*j) * 2.0 * M_PI * N * n);

				NrM2 = lg * *(a_ptr+24+8*j) / (pow( *(a_ptr+23+8*j),3) * *(a_ptr+22+8*j) * 2.0 * M_PI * N * n );




				LNrM= (1.0-lg) * *(a_ptr+24+8*j) * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (N *n); 






				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);

					f2m_der = f2A_der + fA* 2.0 * clb[1] * *(fabDer_ptr+l);

					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr + l) - sh * f2m;

					f1m_der = f1A_der - fA* 2.0 * clb[3] * *(fabDer_ptr + l) - sh* f2m_der;

#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / pow((f2Lb+ pow((*(f2_ptr+i)-f2m),2)),2);
						*(Lf22+i) = 1 /(f2Lb+ pow((*(f2_ptr+i)-f2m),2));


						for (h=0;h<q;h++){



							/* printf("%f , %f\n",f2m,f1m);*/

							/* Combination Lorentzian & Gaussian */
							*(f2der_ptr+i*q+h) = ( NrM * (*(f2_ptr+i)-f2m) * exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + 2.0 * LNrM * (*(f2_ptr+i)-f2m) / (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i) ); 

							*(f1der_ptr+i*q+h) = ( NrM2 * (*(f1_ptr+i)-f1m) * exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + 2.0 * LNrM * (*(f1_ptr+i)-f1m) / pow((f1Lb + pow((*(f1_ptr+h)-f1m),2)),2) * *(Lf22+i) ); 

							*(fn_ptr+i*q+h) +=  (*(f2der_ptr+i*q+h) * f2m_der +  *(f1der_ptr+i*q+h) * f1m_der);


						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/

	/* derivative wrt amp */

	else if (proc_n[0]==12){ 

		for (j=Sstart;j<(Sfinish-Sstart);j++){




			for (l=0;l<N;l++){

				a = 2*M_PI*( (double) (Mod(l*M,N)))/((double) N);


				b = acos(1-((2*((double) l)+1)/((double) N))); 



				*(fab_ptr+l) = (1.0/15120.0)*(-54.0-3.0*pow(*(a_ptr+21+8*j),2)+60.0*(*(a_ptr+21+8*j)*cos(2*a))-(35.0)*pow(*(a_ptr+21+8*j),2)*cos(4.0*a)+(540+30.0*pow(*(a_ptr+21+8*j),2)-480*(*(a_ptr+21+8*j))*cos(2*a)+70*pow(*(a_ptr+21+8*j),2)*cos(4.0*a))*pow(cos(b),2)+(-630-(35.0)*pow(*(a_ptr+21+8*j),2)+420*(*(a_ptr+21+8*j))*cos(2*a)-35.0*pow(*(a_ptr+21+8*j),2)*cos(4*a))*pow(cos(b),4)); 



			}/*end set-up powder loop */

			if (n>1&dt==1){

				cs_ptr=GaussianRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GaussianRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */


			}

			else if (n>1&dt==2){

				cs_ptr=GammaRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= GammaRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else if (n>1&dt==3){

				cs_ptr=LogNormalRNG((a_ptr+17+j*8),(a_ptr+18+j*8), rn_ptr, Lstart, Lfinish);
				cq_ptr= LogNormalRNG((a_ptr+19+j*8),(a_ptr+20+j*8), rn_ptr2, Lstart, Lfinish);
				/* DispVector(cq_ptr,Lfinish-Lstart); */
			}

			else {
				cs_ptr=malloc(sizeof(double));
				cq_ptr=malloc(sizeof(double));

				*cs_ptr=*(a_ptr+17+j*8); 

				*cq_ptr=*(a_ptr+19+j*8);
				/*printf("%f, %f\n",*cs_ptr,*cq_ptr); */
			} 

			for (k=Lstart;k<(Lfinish-Lstart);k++){

				cs=*(cs_ptr+k);
				cq=abs(*(cq_ptr+k));

				Omega=  (cq )  / (2.0 * sp * (2.0 * sp - 1.0)); 



				fA = pow(Omega,2) / (1.0 * lm);

				f2A= cs - fA * ((clb[0] / 10) * (3.0+pow(*(a_ptr+21+8*j),2)));



				f1A = sn * (t-u) * cs + sn2 * fA * ((clb[2] / 10) * (3.0+pow(*(a_ptr+21+8*j),2))); 



				NrM = lg  / (N * n * *(a_ptr+22+8*j)* *(a_ptr+23+8*j) * 2.0 * M_PI);



				LNrM= (1.0-lg) *  *(a_ptr+22+8*j)* *(a_ptr+23+8*j) / (N *n); 

				f2Gb = (2.0 * pow(*(a_ptr+22+8*j),2)); 
				f2Lb = pow(*(a_ptr+22+8*j),2); 

				f1Gb = (2.0 * pow(*(a_ptr+23+8*j),2)); 
				f1Lb = pow(*(a_ptr+23+8*j),2);



				for (l=0; l<N; l++){


					f2m = f2A + fA * 2.0 * clb[1] * *(fab_ptr+l);
					f1m = f1A - fA * 2.0 * clb[3] * *(fab_ptr+l) - sh * f2m;

#pragma omp parallel for private(h,i)

					for (i=0;i<p;i++){



						*(Gf2+i) = exp(-pow((*(f2_ptr+i)-f2m),2)/f2Gb);
						*(Lf2+i) = 1 / (f2Lb+ pow((*(f2_ptr+i)-f2m),2));




						for (h=0;h<q;h++){


							/* printf("%f , %f\n",f2m,f1m);*/


							/* Combination Lorentzian & Gaussian */

							*(fn_ptr+i*q+h) +=  (NrM *exp(-pow((*(f1_ptr+h)-f1m),2) / f1Gb) * *(Gf2+i) + LNrM/ (f1Lb + pow((*(f1_ptr+h)-f1m),2)) * *(Lf2+i) ); 




						}/*end f1 loop */



					}/*end f2 loop */

				}/*end powder loop */



			}/*end distributions loop */

		}/* end sites loop */


	}/* end proc_n loop*/


	if (proc_n[0]>1){
		DispVector(fn_ptr, p*q); 
		if (n>1){
			DispVector(s_ptr, 4*n*s);
		}
	}

#endif

	free(a_ptr);
	free(cs_ptr);
	free(cq_ptr);

	free(f1_ptr);
	free(f2_ptr);

	free(f1der_ptr);
	free(f2der_ptr);
	free(Gf2);
	free(Lf2);

	free(Gf21);
	free(Lf21);
	free(Gf22);
	free(Lf22);




	free(fn_ptr);
	free(s_ptr);

	free(fib_ptr);

	free(fab_ptr);
	free(fabDer_ptr);
	free(rn_ptr);
	free(rn_ptr2);




	if (inpFile==1) fclose(fptr1);

	return 0;
}/*end of main */


int FileOfDoublesSize(FILE *fin)
{
	int z=-1; 
	double x;

	while(!feof(fin)){
		fread(&x, sizeof(double), 1, fin); 
		z++;}

		fseek(fin, 0, SEEK_SET);
		return z;						    
}	      	 

double *BinRead(FILE *fin, int n)
{

	double *x_ptr;

	x_ptr = malloc(n * sizeof(double));


	fread(x_ptr, sizeof(double), n, fin);


	fseek(fin, 0, SEEK_SET);

	return x_ptr;
}



void ErrorMessage(char *str)
{
	printf("Can't open %s.\n",str);
	return;
}

void DispVector(double *ptr, int n)
{
	int i;
	for (i=0; i < n; i++){
		printf(" %e\n",*( ptr + i));
	}
	return;
}





void FileWrite(FILE *fout, double *ptr, int t)
{


	fwrite(ptr, sizeof(double), t, fout);

	return;
}

double *Fib(int r)
{


	double *f_ptr;

	f_ptr=malloc((r+2)*sizeof(double));
	*f_ptr=0;
	*(f_ptr+1)=1;

	for(i=0;i<r;i++){
		*(f_ptr+i+2) = *(f_ptr+i+1) + *(f_ptr+i);


	}

	return f_ptr;
}

double Gam(int n){

	if (n <= 1){
		return 1;
	}     
	else{
		return (double) (n-1)*Gam(n-1);

	}
}
int Mod(int g, int gg){

	int yy;

	yy= g - gg*floor(g/gg);

	return yy;

}

#if 0
double *GaussianRNG(double *ptr, double *ptr2, double *ptr3, int n, int m)
{
	double *ret_ptr=malloc((m-n)*sizeof(double));
	int i;

	for (i = n; i < m; i++)
	{
		*(ret_ptr+i-n)=gsl_cdf_gaussian_Qinv(*(ptr3+i),*ptr2)+*ptr; 
	}

	return ret_ptr;

}

double *GammaRNG(double *ptr, double *ptr2, double *ptr3, int n, int m)
{

	double *ret_ptr=malloc((m-n)*sizeof(double));
	int i;

	for (i = n; i < m; i++)
	{
		*(ret_ptr+i-n)=gsl_cdf_gamma_Qinv(*(ptr3+i),*ptr,*ptr2); 
	}

	return ret_ptr;

}

double *LogNormalRNG(double *ptr, double *ptr2, double *ptr3, int n, int m)
{
	double *ret_ptr=malloc((m-n)*sizeof(double));
	int i;

	for (i = n; i < m; i++)
	{
		*(ret_ptr+i-n)=gsl_cdf_lognormal_Qinv(*(ptr3+i),*ptr,*ptr2); 
	}

	return ret_ptr;

}
#endif

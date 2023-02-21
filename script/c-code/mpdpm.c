/*****************************************/

/* Mixture Prior for Dirichlet Process
   Mixture (MPDPM) of Gaussians */
/* mpdpm.c */

/* an R-callable function */
 
/* by Ken Ryan, Mark Culp */
/* implemented for prediction 
   in double-precision by Cory Lanker */

/*****************************************/
#include "standard.h"
#include "cdflib.h"
#include "randlib.h"
/*****************************************/

/*DYNAMIC MEMORY*/
int   *ivector(int length);
int   **imatrix(int n, int m);
double *fvector(int length);
double **fmatrix(int n, int m);
/*****************************************/
/*RANDOM NUMBERS*/
extern float genbet(float a, float b);
extern float gennor(float av, float sd);
extern float genunf(float low, float high);
extern long  ignuin(long low, long high);

void mvn(int k, double *mu, double *Schol, double *z, double *x);
double iwish(int k, int kk, int rho, double *Sinvchol, double *z, double *x, double *sigmaINV, double *sigma);
/*****************************************/

/*TOLERANCE FOR CHOLESKY*/
#define eta 1.0E-18  // from 1.0E-9, to accommodate double

/*TOLERANCE FOR INVERSE*/
#define tol 6.0E-309 // below this and inverse is inf

/*INDEX FXNS: LOWER TRIANGULAR MATRIX, ROW-WISE*/
#define ISWAP(x,y) {int tmp=(x); (x)=(y); (y)=(tmp);}
int ind(int i, int j){if(i<j) ISWAP(i,j); return i*(i+1)/2+j;}

/*MATRIX OPERATIONS*/
double absf(double x){if(x>0.0) return x; else return -1.0*x;}
double chol(int N, int NDIM, double *A, double *U);
double inv(int n, int nn, double *wn, double *A, double *U);
/*****************************************/
void delete_return(char name[50]){int i; i=strlen(name)-1; if(name[i]=='\n') name[i]='\0';}
void computep(int n, double *phi, double *p);
/*****************************************/
void mpdpm(int *dimint, double *dimdbl, double *datain, double *priorin, int *indout, double *phiout, 
	   double *muout, double *sigmaout, /*int *Mout,*/ char **seedstring)
{
  int i, j, l, k, kk, n, N, B, Burn, Thin, loop, draw=0, *index, /* **M,*/ row, col, *sized, rho; 
  double **x, *phi, **etaMU, **etaSIGMA, **etaINV, *etaRDET, *p, *mu, *xbar;
  double **TEMP, **TEMP2, *size, sizes, **total, *SS, *V, *Vinv, *Vchol, QF;
  double u, a, b, *z, *G, *Gchol, *Ginv, Grdet, *wk, tmpf;
  long seed1, seed2;
  //char name[50], fname[60];
  //FILE *fp, *fp2;

  /*MARCH 2013*/
  int L;
  double  *wl, *cumwl, **S, **Sinvchol, *SRDET, **SSS, *SSSRDET, *TEMPKK;

  /*MISSING DATA*/
  int Ncomplete, Nincomplete, *complete, *rcomplete, *rincomplete, **missing, *nmissing, l1, l2;
  double *mu1, *mu2, *sigma11, **sigma21, *sigma22, *x1;
  //FILE *fpinc;

  /*OCTOBER 2013: conversion of some 'float' variables to type 'double'*/
  //int loopprint=0;
  int iterprint=dimint[8], printflag=dimint[9];
  double wdiv, *TEMPL, *wN, *wn, wvara, wvarb;

  /*PROMPT USER FOR N, k*/
  N = dimint[0];
  if(N<4) exit(0); 
  k = dimint[1];
  if(k<1) exit(0); 
  x = (double **)calloc(N, sizeof(double *));
  for(j = 0; j < N; j++)
    x[j] = datain + k * j;
  index=ivector(N); /*M=imatrix(N,N);*/ mu=fvector(k); xbar=fvector(k); kk=k*(k+1)/2; 
  G=fvector(kk); Gchol=fvector(kk); Ginv=fvector(kk); wk=fvector(k); z=fvector(k);
  SS=fvector(kk); TEMP=fmatrix(k,k); TEMP2=fmatrix(k,k);
  V=fvector(kk); Vchol=fvector(kk); Vinv=fvector(kk);

  /*MISSING DATA*/
  Ncomplete=Nincomplete=0; complete=ivector(N);

  /*PROMPT USER FOR DATA FILE*/
  for(i=0; i<N; ++i)
    {
      complete[i]=1;
      for(j=0; j<k; ++j)
	{
	  if(x[i][j]<-9998.9) complete[i]=0;
	}
      if(complete[i]==0) ++Nincomplete;
      else
	{
	  ++Ncomplete;
	  for(j=0; j<k; ++j) xbar[j]+=x[i][j];
	}
    }

  /*+1 IN FOLLOWING ALLOCATIONS TO CIRCUMVENT ISSUES WHEN NO INCOMPLETE RECORDS*/
  rcomplete   = ivector(Ncomplete+1);     /*ROW INDICES WITH NO   MISSING VALUE(S)*/
  rincomplete = ivector(Nincomplete+1);   /*ROW INDICES WITH SOME MISSING VALUE(S)*/
  missing     = imatrix(Nincomplete+1,k); /*INDICATORS OF MISSING COMPONENT(S)*/
  nmissing    = ivector(Nincomplete+1);   /*ROW SUMS OF MISSING*/

  l1=l2=0; /*l1 INCOMPLETE RECORD, l2 COMPLETE RECORD ENCOUNTERED*/
  for(i=0; i<N; ++i)
    {
      if(complete[i]==1)
	{
	  rcomplete[l2]=i;
	  ++l2;
	}
      else
	{
	  rincomplete[l1]=i;
	  for(j=0; j<k; ++j)
	    {
	      if(x[i][j]<-9998.9){missing[l1][j]=1; ++nmissing[l1];}
	      else missing[l1][j]=0;
	    }
	  ++l1;
	}
    }
  mu1=fvector(k); mu2=fvector(k); x1=fvector(k); /*ALLOCATE FOR MISSING DATA UPDATES*/
  sigma11=fvector(kk); sigma21=fmatrix(k,k); sigma22=fvector(kk);
  if(printflag > 0){
    Rprintf("\nK-Mix program input: N=%d (Ncomplete=%d, Nincomplete=%d), d=%d-variate observations,", 
	    N, Ncomplete, Nincomplete, k);
  }

  /*COMPUTE Gchol, Ginv*/
  for(j=0; j<k; ++j) xbar[j]/=(1.0*Ncomplete);
  for(i=0; i<Ncomplete; ++i)
    {
      l=rcomplete[i];
      for(row=0; row<k; ++row)
	for(col=0; col<=row; ++col)
	  G[ind(row,col)]+=(x[l][row]-xbar[row])*(x[l][col]-xbar[col]);
    }
  u=pow(1.0*Ncomplete,-1.0*(1.0*k+6.0)/(1.0*k+4.0));
  for(i=0; i<kk; ++i) G[i]*=u;
  Grdet=chol(k, kk, G, Gchol);
  Grdet=inv(k, kk, wk, G, Ginv);

  /*PROMPT USER FOR STICK BREAKING PRIOR*/
  n=dimint[2];
  if(n <= 0) exit(0); 
  if(printflag > 0) Rprintf(" n=%d.\n", n);
  phi=fvector(n-1); p=fvector(n); etaMU=fmatrix(n,k); etaSIGMA=fmatrix(n,kk); etaINV=fmatrix(n,kk);
  etaRDET=fvector(n); size=fvector(n); sized=ivector(n); total=fmatrix(n,k);
  a=dimdbl[0];
  if(a<0.0){exit(0);}
  b=dimdbl[1];
  if(b<0.0){exit(0);}
  if(printflag > 1) Rprintf("    Parameters:  a=%f, b=%f, ", a, b);

  /*PROMPT USER FOR PRIOR PARAMETERS*/
  rho=dimint[3];
  if(rho<k) exit(0);
  if(printflag > 1) Rprintf("rho=%d.  ", rho);

  /*MARCH 2013*/
  L=dimint[4];
  S=fmatrix(L,kk); Sinvchol=fmatrix(L,kk); SRDET=fvector(L); wl=fvector(L); cumwl=fvector(L); 
  SSS=fmatrix(L,kk); SSSRDET=fvector(L); TEMPKK=fvector(kk); 

  /*OCTOBER 2013*/
  TEMPL=fvector(L); wN=fvector(N); wn=fvector(n); 
  double *qfs;
  qfs=fvector(n); 

  //tmpf=1.0*(k+rho+1); // THIS IS MODE OF INV-WISH DSN
  // tmpf=1.0*(rho-k-1);  // changed NOV 2014 
  // perhaps (rho-k-1) is too extreme given that the default rho=k+2. Try just rho instead
  //  (rho is average of that which yields the mode and the mean...)
  tmpf=1.0*rho;  // Final mpdpm program value.
  for(l=0; l<L; ++l)
    {
      for(j=0; j<kk; ++j)
	{
	  S[l][j]=priorin[l*(kk+1)+j]*tmpf;
	}
      wl[l]=priorin[l*(kk+1)+j];
      if(wl[l]<=0)
	{ 
	  REprintf("ERROR 1: Variance kernel file has nonpositive weight.\n\n"); 
	  exit(0);
	}
      cumwl[l]=wl[l]; if(l>0) cumwl[l]+=cumwl[l-1];
      SRDET[l]=chol(k,kk,S[l],Sinvchol[l]); /*CHOL WILL TERMINATE PROGRAM IF MATRIX NOT POSITIVE DEFINITE*/
      inv(k,kk,wk,S[l],Sinvchol[l]);
      chol(k,kk,Sinvchol[l],Sinvchol[l]);
    }
  /*normalize cumwl to [0,1], 10/13*/
  for(l=0; l<L; ++l)
    cumwl[l] = (1.0*cumwl[l])/cumwl[L-1];
  if(printflag > 1) Rprintf("Read (SIGMA_l,w_l) for l=1..%d.  ", L);

  /*PROMPT USER FOR SEED AND LENGTH OF MCMC*/
  phrtsd(*seedstring, &seed1, &seed2); setall(seed1, seed2);
  if(printflag > 1) Rprintf("Seed=%s\n", *seedstring);
  Burn = dimint[5];
  if(Burn<0) exit(0); 
  Thin = dimint[6];
  if(Thin<=0) exit(0); 
  B = dimint[7];
  if(B<=0) exit(0); 
  if(printflag > 1) Rprintf("    With Burn=%d, Thin=%d, and Samples=%d, there will %d draws.\n", Burn, Thin, B, Burn+Thin*B);

  /*DEFINITIONS COMPLETE. 
    STARTING METHOD.*/

  /*INITIALIZE SAMPLER USING PRIOR*/
  for(j=0; j<n-1; ++j) {phi[j]=(double)genbet((float)a,(float)b);}
  computep(n, phi, p);

  /* NOT SURE WHY I ADDED THIS SECTION IN OCT 2013... REMOVING. DEC 2014 */
  /* /\*TO CORRECT THE p, 10/13*\/  */
  /* for (j=0; j<n; j++) */
  /*   { */
  /*     p[j]=(double)genunf(0.0,1.0); */
  /*     if (j>0) p[j]+=p[j-1]; */
  /*   } */
  /* wdiv = 1.0/p[n-1]; */
  /* for (j=0; j<n; j++) p[j] *= wdiv; */

  /*INITIALIZE MEANS etaMU and COVARIANCES etaSIGMA*/
  for(j=0; j<n; ++j)
    {
      mvn(k, x[rcomplete[ignuin(0,Ncomplete-1)]], Gchol, z, etaMU[j]); /*MEAN*/
      u=cumwl[L-1]+1;
      while (u >= cumwl[L-1]){
	u=(double)genunf(0.0,(float)cumwl[L-1]); /*cumwl[L-1]=1 now*/
      }
      l=0; while(cumwl[l]<u) ++l; 
      etaRDET[j]=iwish(k, kk, rho, Sinvchol[l], z, wk, etaINV[j], etaSIGMA[j]); /*etaSIGMAj ITS RDET and INVERSE*/
    }

  /*STARTING GIBBS*/
  for(loop=0; loop<Burn+Thin*B; ++loop) /*BEGIN GIBBS*/
    {
      for(j=0; j<n; ++j) /*INITIALIZE SAMPLE SIZES AND TOTALS TO ZERO*/
	{
	  sized[j]=0;  /*SIZES STORED AS INTEGERS*/
	  size[j]=0.0; /*SIZES STORED AS DOUBLES*/
	  for(i=0; i<k; ++i)
	    total[j][i]=0.0;
	}

      /*UPDATE MISSING VALUES: SKIPPING THIS (DO NOT USE MISSING VALUES FOR KBGM)*/
      /*
      for(i=0; i<Nincomplete; ++i) 
	{
	  l1=l2=0; //l1 INDEX OF MISSING VECTOR, l2 INDEX OF COMPLETE VECTOR
	  for(row=0; row<k; ++row) //FILL mu1, mu2, sigma11, sigma21, sigma22
	    {
	      if(missing[i][row]==1)
		{
		  mu1[l1]=etaMU[index[rincomplete[i]]][row];
		  l=0;
		  for(col=0; col<=row; ++col)
		    {
		      if(missing[i][col]==1)
			{
			  sigma11[ind(l1,l)]=etaSIGMA[index[rincomplete[i]]][ind(row,col)];
			  ++l;
			}
		      else 
			sigma21[col-l][l1]=etaSIGMA[index[rincomplete[i]]][ind(row,col)];
		    }
		  ++l1;
		}
	      else 
		{
		  mu2[l2]=x[rincomplete[i]][row]-etaMU[index[rincomplete[i]]][row];
		  l=0;
		  for(col=0; col<=row; ++col)
		    if(missing[i][col]==1)
		      {
			sigma21[l2][l]=etaSIGMA[index[rincomplete[i]]][ind(row,col)];
			++l;
		      }
		    else
		      sigma22[ind(l2,col-l)]=etaSIGMA[index[rincomplete[i]]][ind(row,col)];
		  ++l2;
		}
	    }

	  inv(k-nmissing[i], (k-nmissing[i])*(k-nmissing[i]+1)/2, wk, sigma22, sigma22); //OVERWRITE SIGMA22 WITH SIGMA22INV
	  for(row=0; row<nmissing[i]; ++row) //FILL TEMP WITH t(SIGMA21)*SIGMA22INV
	    for(col=0; col<k-nmissing[i]; ++col)
	      {
		TEMP[row][col]=0.0;
		for(j=0; j<k-nmissing[i]; ++j)
		  TEMP[row][col]+=sigma21[j][row]*sigma22[ind(j,col)];
		mu1[row]-=TEMP[row][col]*mu2[col]; //CONDITIONAL MEAN
	      }
	  for(row=0; row<nmissing[i]; ++row) //FILL SIGMA11 WITH CONDITIONAL VARIANCE
	    for(col=0; col<=row; ++col)
	      {
		l=ind(row,col);
		for(j=0; j<k-nmissing[i]; ++j)
		  sigma11[l]-=TEMP[row][j]*sigma21[j][col];
	      }

	  chol(nmissing[i], nmissing[i]*(nmissing[i]+1)/2, sigma11, sigma11); //GET CHOL(sigma11)
	  mvn(nmissing[i], mu1, sigma11, z, x1); //UPDATE MISSING VALUES
	  l=0;
	  for(j=0; j<k; ++j) //IMPUT
	    if(missing[i][j]==1)
	      {
		x[rincomplete[i]][j]=x1[l];
		++l;
	      }
	      }*/
      /*END SKIPPED MISSING DATA CODE*/

      /* OBSERVATION-BY-OBSERVATION  */
      for(i=0; i<N; ++i) 
	{

	  /*UPDATE INDEX FOR OBS i TO SOME j*/
	  for(j=0; j<n; ++j)
	    {
	      qfs[j]=0.0;
	      l=0; /*index of sigmaINV*/
	      for(row=0; row<k; ++row)
		{
		  tmpf=0.0;
		  for(col=0; col<row; ++col) /*TERMS FROM BELOW MAIN DIAGONAL*/
		    {
		      tmpf+=(x[i][col]-etaMU[j][col])*etaINV[j][l];
		      ++l;
		    }
		  tmpf*=2.0; /*SYMMETRY WITH TERMS FROM ABOVE MAIN DIAGONAL*/

		  u=x[i][row]-etaMU[j][row]; /*MAIN DIAGONAL TERMS*/
		  tmpf+=u*etaINV[j][l];
		  ++l;

		  qfs[j]+=tmpf*u;
		}
	    }
	  tmpf=1.0;
	  wn[n-1]=0.0;
	  while(wn[n-1] < tol)
	    {
	      for(j=0; j<n; ++j)
		{
		  wn[j]=exp(-0.5*tmpf*qfs[j])*(p[j]/etaRDET[j]); /*OCTOBER 2013*/ /*RDET is |Sigma|^0.5 */
		  /* confirmed October 2013 change is needed, before about 20% were 0, now rarely are any 0. */
		  if(j!=0) wn[j]+=wn[j-1];
		}
	      tmpf*=0.5;  /* a correction to find *any* reasonable mixture for this observation
			       rather than random selection with equal probability */
	    }

	  /* /\*normalize to prevent zeroing of small values in conversion to float in genunf*\/ */
	  /* if(wn[n-1] < tol){  /\*OCTOBER 2013*\/ */
	  /*   /\* if ((loopprint % 2) == 0) *\/ */
	  /*   /\*   loopprint = loopprint + 1; *\/ */
	  /*   for(j=0; j<n; ++j) */
	  /*     wn[j] = (j+1.0)/n; */
	  /* } else {   */
	    /* reason for this: wn[n-1] could be less than lowest allowed float value, 
	       causing problems with genunf for u. To prevent, rescale to [0,1] */
	  wdiv=1.0/wn[n-1]; /*normalize wn to [0,1], using wdiv to save flops*/
	  for(j=0; j<n; ++j)
	    wn[j] = wn[j]*wdiv; 
	  //} /*end OCTOBER 2013*/

	  u=wn[n-1]+1.0;  // force while loop to run once
	  while (u >= wn[n-1]){
	    u=(double)genunf(0.0,(float)wn[n-1]); /* wn[n-1]=1 now */
	  }
	  index[i]=0; while(wn[index[i]]<u) ++index[i]; /*UPDATE index[i]*/

	  /*UPDATE SAMPLE SIZE AND TOTAL BASED ON index[i]*/
	  size[index[i]]+=1.0;
	  ++sized[index[i]];
	  for(j=0; j<k; ++j)
	    total[index[i]][j]+=x[i][j];
	}

      /*UPDATE p, phi*/
      sizes=1.0*N; 
      for(j=0; j<n-1; ++j) {sizes-=size[j]; phi[j]=(double)genbet((float)(a+size[j]),(float)(b+sizes));}
      computep(n, phi, p); 

      /*UPDATE ETAj=(muETAj,sigmaETAj) for i=0,1, j=0..n*/
      for(j=0; j<n; ++j)
	{
	  if(sized[j]==0)
	    { /*EMPTY CLUSTER*/

	      mvn(k, x[rcomplete[ignuin(0,Ncomplete-1)]], Gchol, z, etaMU[j]); /*UPDATE etaMUj*/
	      u=cumwl[L-1]+1;
	      while (u >= cumwl[L-1]){
		u=(double)genunf(0.0,(float)cumwl[L-1]); /*cumwl[L-1]=1 now*/
	      }
	      l=0; while(cumwl[l]<u) ++l;
	      etaRDET[j]=iwish(k, kk, rho, Sinvchol[l], z, wk, etaINV[j], etaSIGMA[j]); /*UPDATE etaSIGMAj ITS RDET and INVERSE*/
	    }
	  else 
	    { /*NON-EMPTY CLUSTER*/

	      /*UPDATE VARIANCE MATRIX V FOR MVN UPDATE OF etaMUj*/
	      for(i=0; i<kk; ++i) Vinv[i]=Ginv[i]+size[j]*etaINV[j][i];
	      inv(k, kk, wk, Vinv, V);
	      chol(k, kk, V, Vchol); 
	      for(row=0; row<k; ++row) /*WRITE Ginv*V INTO TEMP*/ 
		for(col=0; col<k; ++col)
		  {
		    TEMP[row][col]=0.0;
		    for(i=0; i<k; ++i)
		      TEMP[row][col]+=Ginv[ind(row,i)]*V[ind(i,col)];
		  }
	      for(row=0; row<k; ++row) /*WRITE TEMP*etaINV INTO TEMP2*/
		for(col=0; col<k; ++col)
		  {
		    TEMP2[row][col]=0.0;
		    for(i=0; i<k; ++i)
		      TEMP2[row][col]+=TEMP[row][i]*etaINV[j][ind(i,col)];
		  }
	      for(i=0; i<k; ++i) xbar[i]=total[j][i]/size[j];
	      for(i=0; i<N; ++i)
		{
		  QF=0.0;
		  for(col=0; col<k; ++col) /*COL OF TEMP2*/
		    {
		      tmpf=0.0;
		      for(row=0; row<k; ++row) /*ROW OF TEMP2*/
			tmpf+=(x[i][row]-xbar[row])*TEMP2[row][col];
		      QF+=tmpf*(x[i][col]-xbar[col]);
		    }
		  wN[i]=exp(-0.5*size[j]*QF);
		  if(i!=0) wN[i]+=wN[i-1];
		}
	      if(wN[N-1] < tol){  /*OCTOBER 2013*/
		/* if (((int)floor(loopprint/2.0) % 2) == 0) */
		/*   loopprint = loopprint + 2; */
		for(i=0; i<N; ++i)
		  wN[i] = (i+1.0)/N; /* here I don't expect this to be a problem. The component is not empty
					so I assume one point will be close enough to xbar 
					(some tests confirm this) */
	      } else {
		wdiv=1/wN[N-1]; /*normalize wN to [0,1], using wdiv to save flops*/
		for(i=0; i<N; ++i)
		  wN[i] = wN[i]*wdiv; 
	      } /*end OCTOBER 2013*/
	      u = wN[N-1]+1.0;
	      while(u >= wN[N-1]){
		u=(double)genunf(0.0,(float)wN[N-1]); 
	      }
	      i=0; while(wN[i]<u) ++i; /*UPDATE INDEX i OF KERNEL MEAN x_i*/

	      l=0;
	      for(row=0; row<k; ++row)
		{
		  wk[row]=0.0;
		  for(col=0; col<row; ++col) /*TERMS FROM BELOW MAIN DIAG AND SYMMETRY WITH ABOVE IT*/
		    {
		      wk[row]+=Ginv[l]*x[i][col]+etaINV[j][l]*total[j][col]; /*BELOW*/
		      wk[col]+=Ginv[l]*x[i][row]+etaINV[j][l]*total[j][row]; /*ABOVE*/
		      ++l;
		    }
		  wk[row]+=Ginv[l]*x[i][row]+etaINV[j][l]*total[j][row]; /*MAIN DIAG row=col*/
		  ++l;
		}
	      l=0;
	      for(row=0; row<k; ++row)
		{
		  mu[row]=0.0;
		  for(col=0; col<row; ++col)
		    {
		      mu[row]+=V[l]*wk[col];
		      mu[col]+=V[l]*wk[row];
		      ++l;
		    }
		  mu[col]+=V[l]*wk[col];
		  ++l;
		}

	      mvn(k, mu, Vchol, z, etaMU[j]); /*UPDATE MEAN VECTOR etaMUj*/

	      /*UPDATE VARIANCE*/
	      for(i=0; i<kk; ++i) SS[i]=0.0;
	      for(i=0; i<N; ++i) 
		if(index[i]==j)
		  {
		    l=0;
		    for(row=0; row<k; ++row)
		      for(col=0; col<=row; ++col)
			{
			  SS[l]+=(x[i][row]-etaMU[j][row])*(x[i][col]-etaMU[j][col]);
			  ++l;
			}
		  }
	      for(l=0; l<L; ++l)
		for(i=0; i<kk; ++i) 
		  SSS[l][i]=S[l][i]+SS[i];

	      for(l=0; l<L; ++l) SSSRDET[l]=chol(k,kk,SSS[l],TEMPKK);

	      /*OCTOBER 2013*/
	      row=0;
	      TEMPL[L-1]=0.0;  // forces loop to run once
	      while(TEMPL[L-1] < tol && row++ < 10)
		{
		  for(l=0; l<L; ++l)
		    {
		      wvara = pow(SRDET[l]/SSSRDET[l],1.0*rho);
		      wvarb = pow(SSSRDET[l],1.0*size[j]);
		      TEMPL[l] = (wvara/wvarb)*wl[l];
		      if(l>0) TEMPL[l]+=TEMPL[l-1];
		    }
		  if (TEMPL[L-1] < tol)
		    {
		      /* this part is only a correction if the S_l matrices yield sum(TEMPL) < tol */
		      for(l=0; l<L; ++l)
			SSSRDET[l]*=0.5;	      
		    }
		}
	      if(TEMPL[L-1] < tol){  /*OCTOBER 2013*/
		/* This is in case TEMPL[L-1] is still less than tol after 10 correction runs above...*/
		/* if (((int)floor(loopprint/4.0) % 2) == 0) */
		/*   loopprint = loopprint + 4; */
		for(l=0; l<L; ++l)
		  TEMPL[l] = (l+1.0)/L; 
	      } else {
		wdiv=1/TEMPL[L-1]; /*normalize TEMPL to [0,1], using wdiv to save flops*/
		for(l=0; l<L; ++l)
		  TEMPL[l] = TEMPL[l]*wdiv;
	      } /*end OCTOBER 2013*/
	      u=TEMPL[L-1]+1.0;
	      while(u >= TEMPL[L-1]){
		u=(double)genunf(0.0,(float)TEMPL[L-1]);
	      }
	      l=0; while(TEMPL[l]<u) ++l;

	      /*INVERT Sl+SSm FOR iwish COVARIANCE MATRIX C CODE*/
	      inv(k, kk, wk, SSS[l], SSS[l]);
	      chol(k, kk, SSS[l], SSS[l]);
	      /*UPDATE etaSIGMAj ITS RDET and INVERSE*/
	      etaRDET[j]=iwish(k, kk, rho+sized[j], SSS[l], z, wk, etaINV[j], etaSIGMA[j]);
	    }
	}

      /*OUTPUT DRAW*/
      if(loop>=Burn && (loop-Burn)%Thin==0)
	{
	  for(i=0; i<N; ++i) indout[draw*N+i]=index[i];
	  for(j=0; j<n-1; ++j) phiout[draw*(n-1)+j]=phi[j];
	  for(j=0; j<n; ++j) 
	    {
	      for(i=0; i<k; ++i)  muout[draw*n*k+j*k+i]=etaMU[j][i]; 
	      for(i=0; i<kk; ++i) sigmaout[draw*n*kk+j*kk+i]=etaSIGMA[j][i];
	    }

	  /* IGNORE NEXT 4 LINES FOR M */
	  /* for(i=0; i<N-1; ++i) /\*UPDATE INDICENCE MATRIX TOTAL*\/ */
	  /*   for(j=i+1; j<N; ++j) */
	  /*     if(index[i]==index[j])  */
	  /* 	++Mout[i*N+j]; */

	  draw++; 
	  /* 	/\*MISSING DATA*\/ */
	  /* 	for(i=0; i<N; ++i) */
	  /* 	  for(j=0; j<k; ++j)  */
	  /* 	    fprintf(fpinc, " %f", x[i][j]); */
	  /* 	fprintf(fpinc, "\n");	 */
	}
      /* if((loop % 500) == 0) */
      /*   if((printflag > 1) && (loopprint > 0)){ printf("e%d:%d ", loopprint, loop); loopprint = 0; } // OCTOBER 2013 */
      if(((loop+1) % iterprint) == 0 & printflag > 1) Rprintf(" -- finished iteration # %d --\n",loop+1); // OCTOBER 2013
    }

  /* IGNORE NEXT 4 LINES FOR M */
  /* for(i=0; i<N-1; ++i) /\*OUTPUT TOTAL INCIDENCE MATRIX*\/ */
  /*   for(j=i+1; j<N; ++j)  */
  /*     Mout[j*N+i]=Mout[i*N+j]; /\*SYMMETRY*\/ */
  /* for(i=0; i<N; ++i) Mout[i*N+i]=B; /\*MAIN DIAGONAL*\/ */

  if(printflag > 0) Rprintf("C program finished normally.\n");
  //exit(0);
}
/*****************************************/
int *ivector(int length)
{
  int i, *x;

  if((x = (int*)malloc((unsigned)(length * sizeof(int))))==NULL){
    REprintf("\nError 2A: allocation failure in ivector().\n\n"); exit(-1);}

  for(i=0; i<length; ++i) x[i]=0;
 
  return x;
}
/*****************************************************/
int **imatrix(int n, int m)
{
  int i, j, **x, *xnow;

  if((x = (int**)malloc((unsigned)(n * sizeof(int*))))==NULL){
    REprintf("\nError 2B: allocation failure 1 in matrix().\n\n"); exit(-1);}

  if((xnow = (int*)malloc((unsigned)(n * m * sizeof(int))))==NULL){
    REprintf("\nError 2B: allocation failure 2 in matrix().\n\n"); exit(-1);}

  for(i=0;i<n;i++,xnow += m)x[i] = xnow;

  for(i=0; i<n; ++i)
    for(j=0; j<m; ++j)
      x[i][j]=0;

  return x;
}
/*****************************************************/
double *fvector(int length)  /*OCTOBER 2013*/
{
  int i;
  double *x;

  if((x = (double*)malloc((unsigned)(length * sizeof(double))))==NULL){
    REprintf("\nError 2C: allocation failure in fvector().\n\n"); exit(-1);}

  for(i=0; i<length; ++i) x[i]=0.0;
 
  return x;
}
/*****************************************************/
double **fmatrix(int n, int m)  /*OCTOBER 2013*/
{
  int i, j;
  double **x, *xnow;

  if((x = (double**)malloc((unsigned)(n * sizeof(double*))))==NULL){
    REprintf("\nError 2D: allocation failure 1 in fmatrix().\n\n"); exit(-1);}

  if((xnow = (double*)malloc((unsigned)(n * m * sizeof(double))))==NULL){
    REprintf("\nError 2D: allocation failure 2 in fmatrix().\n\n"); exit(-1);}

  for(i=0;i<n;i++,xnow += m)x[i] = xnow;

  for(i=0; i<n; ++i)
    for(j=0; j<m; ++j)
      x[i][j]=0.0;

  return x;
}
/*****************************************************/
void computep(int n, double *phi, double *p)
{
  int j;

  p[0]=1.0;
  for(j=0; j<n-1; ++j)
    {
      p[j+1]=p[j]*(1.0-phi[j]);
      p[j]*=phi[j];
    }
}
/*****************************************************/
double chol(int N, int NDIM, double *A, double *U)
{
  int irow, icol, i, j=0, k=-1, l, m;
  double w, rdet;

  /*MODIFIED AS6 FROM STATLIB FOR A POSITIVE DEFINITE INPUT MATRIX A*/

  /*THIS FXN TERMINATES THE PROGRAM IF THE SYMMETRIC MATRIX A (GIVEN 
    IN LOWER TRIANGULAR FORM BY ROWS) IS NOT POSITIVE SEMI DEFINITE 
    WITHIN TOLERANCE ETA>0 OR IF THE DIMENSION OF ITS NULL SPACE IS
    NOT 0.*/

  /*OTHERWISE, LOWER TRIANGULAR MATRIX U (GIVEN IN LOWER TRIANGULAR 
    FORM BY ROWS) IS THE CHOLESKY DECOMPOSITION OF A, I.E., A=UU'.
    A AND U MAY COINCIDE IF YOU WANT TO OVERWRITE A WITH U.
    THE FXN RETURNS sqrt(det(A))=det(U).*/

  rdet=1.0;
  for(icol=0; icol<N; ++icol)
    {
      l=-1;
      for(irow=0; irow<=icol; ++irow)
	{
	  ++k;
	  w=A[k];
	  m=j;
	  for(i=0; i<irow; ++i)
	    {
	      ++l;
	      w-=U[l]*U[m]; 
	      ++m;
	    }
	  ++l;
	  if(irow!=icol)
	    {if(absf(U[l])>eta) U[k]=w/U[l]; else U[k]=0.0;}
	}

      if(w >= absf(eta*A[k])) U[k]=sqrt(w);  /*Matrix MAY be NND and invertible*/ /*This line is different from other files.*/
      else if(absf(w)<absf(eta*A[k]))
	{
	  REprintf("\nERROR 3A: Call to fxn chol with a matrix that is not invertible.\n");
	  exit(-1);
	}
      else 
	{
	  REprintf("\nERROR 3B: Call to fxn chol with a matrix that is not NND.\n");
	  exit(-1);
	}
      j+=icol+1;
      rdet*=U[k];
    }

  return rdet;
}
/*****************************************************/
double inv(int n, int nn, double *wn, double *A, double *U)
{
  int irow, icol, jcol, i, j, k, ndiag, mdiag, go;
  double x, rdet;

  /*MODIFIED AS7 FROM STATLIB FOR A POSITIVE DEFINITE INPUT MATRIX A*/

  /*n MATRIX DIMENSION, nn=n+(n choose 2), wn A WORK VECTOR OF
    LENGTH AT LEAST n*/

  /*THIS FXN TERMINATES THE PROGRAM IF THE SYMMETRIC MATRIX A (GIVEN 
    IN LOWER TRIANGULAR FORM BY ROWS) IS NOT POSITIVE SEMI DEFINITE 
    WITHIN TOLERANCE ETA>0 OR IF THE DIMENSION OF ITS NULL SPACE IS
    NOT 0 BY WAY OF CALLING MODIFIED AS 6.*/

  /*OTHERWISE, THE SYMMETRIC MATRIX U (GIVEN IN LOWER TRIANGULAR FORM 
    BY ROWS) IS THE INVERSE OF A.  A AND U MAY COINCIDE IF YOU WANT TO 
    OVERWRITE A WITH U. THE FXN RETURNS sqrt(det(A))=1/sqrt(det(U)).*/

  rdet=chol(n, nn, A, U);
  --n;
  --nn;
  irow=n;
  ndiag=nn;
  while(irow!=-1)
    {
      j=ndiag;
      if(absf(U[ndiag])>eta)
	{
	  for(i=irow; i<=n; ++i){wn[i]=U[j]; j+=i+1;}
	  icol=n;
	  jcol=mdiag=nn;
	  go=1;
	  while(go)
	    {
	      j=jcol;
	      x=0.0;
	      if(icol==irow) x=1.0/wn[irow];
	      k=n;
	      while(k!=irow)
		{
		  x-=wn[k]*U[j];
		  --k;
		  --j;
		  if(j>mdiag) j-=k;
		}
	      U[j]=x/wn[irow];
	      if(icol==irow) go=0;
	      else {mdiag-=icol+1; --icol; --jcol;}
	    }
	}
      else for(i=irow; i<=n; ++i){U[j]=0.0; j+=i+1;}
      ndiag-=irow+1;
      --irow;
    }

return rdet;
}
/*****************************************************/
void mvn(int k, double *mu, double *Schol, double *z, double *x)
{
  int i, j, l;

  l=0;
  for(j=0; j<k; ++j)
    {
      z[j]=(double)gennor(0.0,1.0);
      x[j]=mu[j]; 
      for(i=0; i<=j; ++i) 
	{
	  x[j]+=Schol[l]*z[i];
	  ++l;
	}
    }
}
/*****************************************************/
double iwish(int k, int kk, int rho, double *Sinvchol, double *z, double *x, double *sigmaINV, double *sigma)
{
  int i,j,l,loop;

  for(i=0; i<kk; ++i) sigmaINV[i]=0.0;
  for(loop=0; loop<rho; ++loop)
    {
      l=0;
      for(j=0; j<k; ++j) /*PUT MVNk(0,S) DRAW INTO x*/
	{
	  z[j]=(double)gennor(0.0,1.0);
	  x[j]=0.0; 
	  for(i=0; i<=j; ++i)
	    {
	      x[j]+=Sinvchol[l]*z[i];
	      ++l;
	    }
	}

      l=0;
      for(i=0; i<k; ++i) /*ADD x*x' TO sigmaINV*/
	for(j=0; j<=i; ++j)
	  {
	    sigmaINV[l]+=x[i]*x[j];
	    ++l;
	  }
    }

  return 1.0/inv(k, kk, x, sigmaINV, sigma); /*CONVERT ROOT DET OF INVERSE TO ROOT DET*/
}
/*****************************************************/

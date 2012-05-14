#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

double d1_getLambda(double t,double lambda_0,double alpha,double beta,SEXP history)
{
	double res = lambda_0;
	SEXP tail = history;
	for (int i=0;i<length(history);i++)
	{
		if (REAL(CAR(tail))[0] <= t)
		{
			res = res + alpha * exp(-beta * (t - REAL(CAR(tail))[0]));
			tail = CDR(tail);
		}
	}
	return res;
}

SEXP d1_Simulate(SEXP Lambda_0,SEXP Alpha,SEXP Beta,SEXP TT)
{
	double lambda_0 = REAL(coerceVector(Lambda_0,REALSXP))[0];
	double alpha = REAL(coerceVector(Alpha,REALSXP))[0];
	double beta = REAL(coerceVector(Beta,REALSXP))[0];
	double T = REAL(coerceVector(TT,REALSXP))[0];

	GetRNGstate(); 
	SEXP block;
	SEXP rootHistory = PROTECT(list1(allocVector(REALSXP, 1)));
	SEXP tailHistory = rootHistory;
	double lambda_star = lambda_0;
	//first event
	double U = unif_rand();
	double s = -(1.0 / lambda_star) * log(U);
	if (s <= T)
	{
		REAL(CAR(tailHistory))[0]=s;
		tailHistory = SETCDR(tailHistory, list1(allocVector(REALSXP, 1)));

	}
	else
	{
		PutRNGstate();
		int n = length(rootHistory);
		SEXP res = PROTECT(allocVector(REALSXP,n-1));
		double *xres=REAL(res);
		for(int i=0;i<n-1;i++)
		{
			xres[i] = REAL(CAR(rootHistory))[0];
			rootHistory = CDR(rootHistory);

		}
		UNPROTECT(2);
		return res;
	}
	//general routine
	while (1)
	{
		lambda_star = d1_getLambda(s,lambda_0,alpha,beta,rootHistory);
		U = unif_rand();
		s = s - (1.0 / lambda_star) * log(U);
		if (s > T)
		{
			PutRNGstate();
			int n = length(rootHistory);
			SEXP res = PROTECT(allocVector(REALSXP,n-1));
			double *xres=REAL(res);
			for(int i=0;i<n-1;i++)
			{
				xres[i] = REAL(CAR(rootHistory))[0];
				rootHistory = CDR(rootHistory);

			}
			UNPROTECT(2);
			return res;
		}
		double D = unif_rand();
		if (D <= d1_getLambda(s,lambda_0,alpha,beta,rootHistory) / lambda_star)
		{
			REAL(CAR(tailHistory))[0]=s;
			tailHistory = SETCDR(tailHistory, list1(allocVector(REALSXP, 1)));

		}
	}

}

struct H{int n;double *history;};
double d1_likelihoodCalc(int N,double* parameters,void *ex)// minus likelihood
{
	double lambda_0 = parameters[0];
	double alpha = parameters[1];
	double beta = parameters[2];

	struct H * h = (struct H *)ex;	
	int nh = h->n;
	double * history = h->history;
	double T = history[nh-1];
	//////////////////////////////////////
	/////////////////////////////////////

	double sum = 0.0;
	for (int i = 0; i < nh; i++)
	{
		sum = sum+   (1 - exp(-beta * (T - history[i]))) ;
	}
	sum = (alpha / beta) * sum;
	double integratedDensity = lambda_0 * T + sum;
	double res = T - integratedDensity;

	for (int i = 0; i < nh; i++)
	{
		sum = lambda_0;
		for (int k = 0; k < i; k++)
		{
			sum = sum + alpha*exp(-beta*(history[i]-history[k]));
		}
		res = res + log(sum);
	}
	res = -1.0*res;
	return (res);
}        
SEXP d1_likelihood(SEXP parameters,SEXP History)// minus likelihood
{
	SEXP params = PROTECT(allocVector(REALSXP,3));
	double *xparams = REAL(params);

	xparams[0] = REAL(coerceVector(parameters,REALSXP))[0];
	xparams[1] = REAL(coerceVector(parameters,REALSXP))[1];
	xparams[2] = REAL(coerceVector(parameters,REALSXP))[2];

	PROTECT(History = coerceVector(History,REALSXP));	
	int n = length(History);
	struct H h;
	h.n = n;
	h.history = REAL(History);			
	SEXP res = PROTECT(allocVector(REALSXP,1));
	double *xres = REAL(res);
	xres[0] =  d1_likelihoodCalc(3,REAL(params),&h);
	UNPROTECT(3);
	return (res);
}
SEXP d1_MLE_calibrate(SEXP History,SEXP initialParams)
{
	SEXP params = PROTECT(allocVector(REALSXP,3));
	double *xparams = REAL(params);

	xparams[0] = REAL(coerceVector(initialParams,REALSXP))[0];
	xparams[1] = REAL(coerceVector(initialParams,REALSXP))[1];
	xparams[2] = REAL(coerceVector(initialParams,REALSXP))[2];

	PROTECT(History = coerceVector(History,REALSXP));	
	int n = length(History);
	struct H h;
	h.n = n;
	h.history = REAL(History);
	SEXP res = PROTECT(allocVector(REALSXP,6));
	double *xres = REAL(res);
	double Fmin=0.0;int fail=0;double abstol = -1e16;double intol = 1e-6;
	int fncount=0;int maxit = 1000;
	double alpha=1.0; double beta=0.5; double gamma=2.0; int trace=0;

	nmmin(3, xparams, xres, &Fmin, d1_likelihoodCalc,
		&fail, abstol,  intol, (void *)&h,
		alpha, beta,  gamma,  trace,
		&fncount, maxit);

	/*
	dyn.load('Hawkes.so')
	a =  .Call('d1_Simulate',as.double(1.0),as.double(0.5),as.double(0.8),as.double(100))
	b= .Call('d1_likelihood',c(1,0.5,0.8),a)
	d=.Call('d1_hawkes1d',a,c(1.1,0.6,0.9))
	*/
	xres[3] = Fmin;
	xres[4] = fncount;
	xres[5] = fail;
	UNPROTECT(3);
	return res;
}


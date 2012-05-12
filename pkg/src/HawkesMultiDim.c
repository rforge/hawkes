#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

struct HHH{double T;int dim;int *lengths;double **history;};
double dn_getLambda(double t, int index,int dimension,double *lambda_0,double **alpha,double **beta,SEXP history1)
{
	double res = lambda_0[index];
	double **history = (double **)R_alloc(dimension,dimension*sizeof(double *));
	int *lengths =(int *)R_alloc(dimension,dimension*sizeof(int));
	//history1 = coerceVector(VECSXP
	for(int i=0;i<dimension;i++)
	{
		int n = length(VECTOR_ELT(history1,i));
		lengths[i] = n;
		history[i] = (double*)R_alloc(n,n*sizeof(double));
		SEXP block = VECTOR_ELT(history1,i);
		for (int j=0;j<n;j++)
		{
        	history[i][j] = REAL(CAR(block))[0];
			block = CDR(block);		
		}
	}

	for (int m = 0; m < dimension; m++)
	{
		//SEXP tail = history[m];
		for (int i=0;i<lengths[m];i++)
		{
			if (history[m][i] <= t)
			{
				res = res + alpha[index][m] * exp(-beta[index][m] * (t - history[m][i]));
			}
		}
	}
	return res;
}

int dn_attribute(double alea, double t, double I_star,int dimension,double *lambda_0,double **alpha,double **beta,SEXP history)
{
	int index = 0;
	double cumul = dn_getLambda(t, 0,dimension,lambda_0,alpha,beta,history);
	while (alea > (cumul / I_star))
	{
		index = index + 1;
		cumul = cumul + dn_getLambda(t, index,dimension,lambda_0,alpha,beta,history);
	}
	return (index);
}

double dn_partialLikelihood(int m,double T,int dimension,double *lambda_0,double **alpha,double **beta,int *lengths,double**history)
{
	double sum = 0.0;
	for (int n = 0; n < dimension; n++)
	{
		for (int k = 0; k < lengths[n]; k++)
		{

			sum = sum + (alpha[m][n] / beta[m][n]) *
				(1-exp(-beta[m][n] * (T - history[n][k])));

		}
	}

	double integratedDensity = lambda_0[m] * T + sum;
	double res = T - integratedDensity;



	for (int i = 0; i < lengths[m]; i++)
	{
		sum = lambda_0[m];
		for (int n = 0; n < dimension; n++)
		{
			for (int k = 0; k < lengths[n]; k++)
			{
				if (history[n][k] < history[m][i])
				{
					sum = sum + alpha[m][n] * exp(-beta[m][n] * (history[m][i] - history[n][k]));
				}
			}
		}
		res = res + log(sum);
	}

	return res;
}


double dn_likelihood(int N,double* parameters,void *ex)// minus likelihood
{
	struct HHH* h = (struct HHH *)ex;
	int dimension = h->dim;
	double * lambda_0 = (double *) R_alloc(dimension, sizeof(double));
	double **alpha = (double **) R_alloc(dimension, sizeof(double*));
	double **beta = (double **) R_alloc(dimension, sizeof(double*));
	for (int i = 0; i < dimension; i++)
	{
			lambda_0[i] = parameters[i];
			alpha[i] = (double *) R_alloc(dimension, sizeof(double));
			beta[i] = (double *) R_alloc(dimension, sizeof(double));
			
			for (int n = 0; n < dimension; n++)
			{
				alpha[i][n] = parameters[dimension + i * dimension + n];
				beta[i][n] = parameters[dimension + dimension * dimension + i * dimension + n];
			}
	}
	double res = 0.0;
	for (int i = 0; i < dimension; i++)
	{
		res = res + dn_partialLikelihood(i,h->T,dimension,lambda_0,alpha,beta,h->lengths,h->history);
	}
	return (-res);
}

SEXP  dn_Simulate(SEXP Lambda_0,SEXP Alpha,SEXP Beta,SEXP TT)
{
	GetRNGstate();
	double *lambda_0 = REAL(coerceVector(Lambda_0,REALSXP));
	double *lalpha = REAL(coerceVector(Alpha,REALSXP));
	double *lbeta = REAL(coerceVector(Beta,REALSXP));
	double T = REAL(coerceVector(TT,REALSXP))[0];
	int dimension = length(Lambda_0);
	double **alpha = (double **) R_alloc(dimension, sizeof(double*));
	double **beta = (double **) R_alloc(dimension, sizeof(double*));
	for(int i=0;i<dimension;i++)
	{
		alpha[i] = (double *) R_alloc(dimension, sizeof(double));
		beta[i] = (double *) R_alloc(dimension, sizeof(double));
		for(int j=0;j<dimension;j++)
		{
			alpha[i][j] = lalpha[i*dimension+j];
			beta[i][j] = lbeta[i*dimension+j];
		}
	}


	GetRNGstate();
	SEXP block;
	SEXP rootHistory;
	SEXP tailHistory;
	PROTECT(rootHistory=allocVector(VECSXP,dimension));
	PROTECT(tailHistory=allocVector(VECSXP,dimension));
	SEXP xrootHistory;
	for(int i=0;i<dimension;i++)
	{
		SET_VECTOR_ELT(rootHistory, i, block = list1(allocVector(REALSXP, 1)));
		SET_VECTOR_ELT(tailHistory, i, block);
	}

	double lambda_star = 0;
	for(int i=0;i<dimension;i++)
	{
		lambda_star += lambda_0[i];
	}
	//first event
	double U = unif_rand();
	double s = -(1.0 / lambda_star) * log(U);
	if (s <= T)
	{
		double D = unif_rand();
		int n0 = dn_attribute(D, 0, lambda_star,dimension,lambda_0,alpha,beta,rootHistory);
		REAL(CAR(VECTOR_ELT(tailHistory,n0)))[0]=s;
		SET_VECTOR_ELT(tailHistory,n0,SETCDR(VECTOR_ELT(tailHistory,n0), list1(allocVector(REALSXP, 1))));
	}
	else
	{
		PutRNGstate();
		SEXP res = PROTECT(allocVector(VECSXP,dimension));
		for(int i=0;i<dimension;i++)
		{
			int n = length(VECTOR_ELT(rootHistory,i));
			SET_VECTOR_ELT(res,i,allocVector(REALSXP,n-1));
			double *xres=REAL(VECTOR_ELT(res,i));
			for(int i=0;i<n-1;i++)
			{
				xres[i] = REAL(CAR(VECTOR_ELT(rootHistory,i)))[0];
				SET_VECTOR_ELT(rootHistory,i,CDR(VECTOR_ELT(rootHistory,i)));

			}
		}
		UNPROTECT(2);
		return res;
	}
	//general routine
	while (1)
	{
		//	REAL(CAR(tailHistory))[0]=s;
		//	tailHistory = SETCDR(tailHistory, list1(allocVector(REALSXP, 1)));

		lambda_star = 0;
		for (int i = 0; i < dimension; i++)
		{//dn_getLambda(double t, int index,int dimension,double *lambda_0,double **alpha,double **beta,SEXP history)
			lambda_star = lambda_star + dn_getLambda(s, i,dimension,lambda_0,alpha,beta,rootHistory);
		}
		U = unif_rand();
		s = s - (1.0 / lambda_star) * log(U);
		if (s <= T)
		{
			double D = unif_rand();
			double I_M = 0.0;
			for (int i = 0; i < dimension; i++)
			{
				I_M = I_M + dn_getLambda(s, i,dimension,lambda_0,alpha,beta,rootHistory);
			}
			if (D <= I_M / lambda_star)
			{
				int n0 = dn_attribute(D, s, lambda_star,dimension,lambda_0,alpha,beta,rootHistory);
				REAL(CAR(VECTOR_ELT(tailHistory,n0)))[0]=s;
				SET_VECTOR_ELT(tailHistory,n0,SETCDR(VECTOR_ELT(tailHistory,n0), list1(allocVector(REALSXP, 1))));
			}
		}
		else
		{
			PutRNGstate();
		SEXP res = PROTECT(allocVector(VECSXP,dimension));
		for(int i=0;i<dimension;i++)
		{
			int n = length(VECTOR_ELT(rootHistory,i));
			SET_VECTOR_ELT(res,i,allocVector(REALSXP,n-1));
			double *xres=REAL(VECTOR_ELT(res,i));
			for(int i=0;i<n-1;i++)
			{
				xres[i] = REAL(CAR(VECTOR_ELT(rootHistory,i)))[0];
				SET_VECTOR_ELT(rootHistory,i,CDR(VECTOR_ELT(rootHistory,i)));

			}
		}
		UNPROTECT(2);
		return res;
		}
	}
}
SEXP dn_hawkesnd(SEXP dim,SEXP TT,SEXP History,SEXP lengths,SEXP initialParams)
{
	int dimension = INTEGER(coerceVector(dim,INTSXP))[0];
	double T = REAL(coerceVector(TT,REALSXP))[0];
	SEXP params = PROTECT(allocVector(REALSXP,dimension+2*dimension*dimension));
	double *xparams = REAL(params);
 	for (int i = 0; i < dimension+2*dimension*dimension; i++)
	{
    	xparams[i] = REAL(coerceVector(initialParams,REALSXP))[i];
    }
	PROTECT(lengths = coerceVector(lengths,INTSXP));
	PROTECT(History = coerceVector(History,REALSXP));	
	int n = length(History);
	struct HHH h;
	h.T = T;
	h.dim = dimension;
	h.lengths=INTEGER(lengths);
	double **history = (double **)R_alloc(dimension,dimension*sizeof(double *));
	int doneSoFar=0;
	for (int i = 0; i < dimension+2*dimension*dimension; i++)
	{
		int n = REAL(coerceVector(lengths,REALSXP))[i];
		history[i] = (double*)R_alloc(n,n*sizeof(double));
		for (int j=0;j<n;j++)
		{
        	history[i][j] = REAL(coerceVector(History,REALSXP))[doneSoFar+i];
			doneSoFar++;
			
		}
	}
	h.history = history;
	SEXP res = PROTECT(allocVector(REALSXP,dimension+2*dimension*dimension+3));
	double *xres = REAL(res);
	double Fmin=0.0;int fail=0;double abstol = -1e16;double intol = 1e-6;
	int fncount=0;int maxit = 1000;
	double alpha=1.0; double beta=0.5; double gamma=2.0; int trace=0;

	nmmin(dimension+2*dimension*dimension, xparams, xres, &Fmin, dn_likelihood,
		&fail, abstol,  intol, (void *)&h,
		alpha, beta,  gamma,  trace,
		&fncount, maxit);


	xres[dimension+2*dimension*dimension] = Fmin;
	xres[dimension+2*dimension*dimension+1] = fncount;
	xres[dimension+2*dimension*dimension+2] = fail;
	UNPROTECT(4);
	return res;
}



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;




int Attribute(double alea, double t, double I_star,const arma::vec& m_lambda)
{
  int index = 0;
	double cumul = m_lambda[0];
	while (alea > (cumul / I_star))
	{
		index = index + 1;
		cumul = cumul + m_lambda[index];
	}
	return (index);
}

// [[Rcpp::export]]
std::vector<std::vector<double> > SimulateHawkes(SEXP lambda0,SEXP alpha,SEXP beta,SEXP horizon)
{
  Rcpp::NumericVector lambda0_internal(lambda0);  
  int dimension = lambda0_internal.size();
  double m_horizon = as<double>(horizon);
  
  std::vector<std::vector<double> > history;
  for (int i=0;i<dimension;i++)
  {
    std::vector<double> a;
    history.push_back(a);
  }
  if (dimension == 1)
  {
    double m_lambda0 = as<double>(lambda0);
    double m_alpha = as<double>(alpha);
    double m_beta = as<double>(beta);
    
    double lambda_star = m_lambda0;
	  double dlambda = 0.0,t=0;
	  //first event
	  double U = arma::randu(1)[0];
	  double s = -(1.0 / lambda_star) * log(U);
	  if (s <= m_horizon)
	  {
		  history[0].push_back(s);
		  dlambda = m_alpha;
		  t = s;
	  }
  	else
  	{
  		return (history);
  	}
  	//general routine
  	while (true)
  	{
  		lambda_star = m_lambda0+dlambda*exp(-m_beta*(s-t));
  		U = arma::randu(1)[0];
  		s = s - (1.0 / lambda_star) * log(U);
  		if (s > m_horizon)
  		{
  			return (history);
  		}
  		double D = arma::randu(1)[0];
  		if (D <= (m_lambda0+dlambda*exp(-m_beta*(s-t))) / lambda_star)
  		{
  			history[0].push_back(s);
  			dlambda = dlambda*exp(-m_beta*(s-t)) + m_alpha;
  			t=s;
  		}
  	}
  }
  else{
    arma::mat dlambda(dimension,dimension);
    
    Rcpp::NumericMatrix alpha_internal(alpha);
    Rcpp::NumericVector beta_internal(beta);
    
    arma::vec m_lambda0(lambda0_internal.begin(),dimension,false);
    arma::mat m_alpha(alpha_internal.begin(),dimension,dimension,false);
    arma::vec m_beta(beta_internal.begin(),dimension,false);
    arma::vec m_lambda(dimension);
  	double lambda_star = 0.0;
  	
  	double t=0;
  	for (int i = 0; i < dimension; i++)
  	{
  		lambda_star += m_lambda0[i];
  		m_lambda[i] = m_lambda0[i];
  	}
    
  	//first event
  	double U = arma::randu(1)[0];
  	double s = -(1.0 / lambda_star) * log(U);
    
    
  	if (s <= m_horizon)
  	{
  		double D = arma::randu(1)[0];
  		int n0 = Attribute(D, 0, lambda_star,m_lambda);
  		history[n0].push_back(s);
  
  		for (int i=0;i<dimension;i++)
  		{
  			dlambda(i,n0) = m_alpha(i,n0);
  			m_lambda[i] = m_lambda0[i]+m_alpha(i,n0);
  		}
  	}
  	else
  	{
  		return (history);
  	}
  	t=s;
    //general routine
  	lambda_star = 0;
  	for (int i = 0; i < dimension; i++)
  	{
  		lambda_star = lambda_star + m_lambda[i];
  	}
  	while (TRUE)
  	{
  		U = arma::randu(1)[0];
  		s = s - (1.0 / lambda_star) * log(U);
  		if (s <= m_horizon)
  		{
  			double D = arma::randu(1)[0];
  			double I_M = 0.0;
  			for (int i = 0; i < dimension; i++)
  			{
  				double dl = 0.0;
  				for (int j = 0; j < dimension; j++)
  				{
  					dl += dlambda(i,j)*exp(-m_beta(i)*(s-t));
  				}
  				m_lambda[i] = m_lambda0[i]+dl;
  				I_M = I_M + m_lambda[i];
  			}
  			if (D <= (I_M / lambda_star))
  			{
  				int n0 = Attribute(D, s, lambda_star,m_lambda);
  				history[n0].push_back(s);
  				lambda_star=0.0;
  				for (int i=0;i<dimension;i++)
  				{
  					double dl=0.0;
  					for (int j = 0; j < dimension; j++)
  					{
  						dlambda(i,j) = dlambda(i,j)*exp(-m_beta(i)*(s-t));
  						if (n0==j)
  						{
  							dlambda(i,n0) += m_alpha(i,n0);
  						}
  						dl +=dlambda(i,j);
  					}
  					lambda_star+= m_lambda0[i]+dl;
  				}
  				t=s;
  			}
  			else
  			{
  				lambda_star = I_M;
  			}
  
  		}
  		else
  		{
  			return (history);
  		}
  	}
  }
}



/*
library(Rcpp)
compileAttributes('C:\\Users\\riadh\\Dropbox\\dev_R\\hawkes\\pkg\\')
library(devtools)
install(pkg='C:\\Users\\riadh\\Dropbox\\dev_R\\hawkes\\pkg')
library(hawkes)
*/
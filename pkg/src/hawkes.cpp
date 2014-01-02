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
std::vector<std::vector<double> > SimulateHawkes(arma::vec& a_lambda0,arma::mat& a_alpha,arma::vec& a_beta,double horizon)
{
  
  int m_dimension = a_alpha.n_elem;
  std::vector<std::vector<double> > m_history;
  for (int i=0;i<m_dimension;i++)
  {
    std::vector<double> a;
    m_history.push_back(a);
  }
  
  arma::mat dlambda(m_dimension,m_dimension);
  arma::vec m_lambda0(a_lambda0);
  arma::mat m_alpha(a_alpha);
  arma::vec m_beta(a_beta);
  arma::vec m_lambda(m_dimension);
	double lambda_star = 0.0;
	
	double t=0;
	for (int i = 0; i < m_dimension; i++)
	{
		lambda_star += m_lambda0[i];
		m_lambda[i] = m_lambda0[i];
	}
  
	//first event
	double U = arma::randu(1)[0];
	double s = -(1.0 / lambda_star) * log(U);
  
  
	if (s <= horizon)
	{
		double D = arma::randu(1)[0];
		int n0 = Attribute(D, 0, lambda_star,m_lambda);
		m_history[n0].push_back(s);

		for (int i=0;i<m_dimension;i++)
		{
			dlambda(i,n0) = m_alpha(i,n0);
			m_lambda[i] = m_lambda0[i]+m_alpha(i,n0);
		}
	}
	else
	{
		return (m_history);
	}
	t=s;
  //general routine
	lambda_star = 0;
	for (int i = 0; i < m_dimension; i++)
	{
		lambda_star = lambda_star + m_lambda[i];
	}
	while (TRUE)
	{
		U = arma::randu(1)[0];
		s = s - (1.0 / lambda_star) * log(U);
		if (s <= horizon)
		{
			double D = arma::randu(1)[0];
			double I_M = 0.0;
			for (int i = 0; i < m_dimension; i++)
			{
				double dl = 0.0;
				for (int j = 0; j < m_dimension; j++)
				{
					dl += dlambda(i,j)*exp(-m_beta(i)*(s-t));
				}
				m_lambda[i] = m_lambda0[i]+dl;
				I_M = I_M + m_lambda[i];
			}
			if (D <= (I_M / lambda_star))
			{
				int n0 = Attribute(D, s, lambda_star,m_lambda);
				m_history[n0].push_back(s);
				lambda_star=0.0;
				for (int i=0;i<m_dimension;i++)
				{
					double dl=0.0;
					for (int j = 0; j < m_dimension; j++)
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
			return (m_history);
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
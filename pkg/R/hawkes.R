Simulate <- function(dim,lambda0,alpha,beta,T)
{
	if(!is.integer(dim) || (dim < 1))
	{
	  stop("'dim' must be a strictly positive integer")
	}
	if(dim==1)
	{
		t = .Call("d1_Simulate",lambda0,alpha,beta,T)
	}
	if(dim > 1)
	{
		t = .Call("dn_Simulate",lambda0,alpha,beta,T)
	}
	t
}

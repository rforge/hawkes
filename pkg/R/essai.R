
.First.lib <- function(libname, pkgname)
{
	library.dynam("Hawkes")
   
    cat("Hawkes process simulation and calibration toolkit.\n")
   


}


Simulate <- function(dim,lambda0,alpha,beta,T)
{
	if (dim < 1)
	{
		cat("Make sure the dimension is a positive integer.\n")
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


.First.lib <- function(libname, pkgname)
{
	library.dynam("Hawkes")
   
    cat("Hawkes process simulation and calibration toolkit.\n")
   


}


f <- function(x)
{
	.C("Simulate")
}

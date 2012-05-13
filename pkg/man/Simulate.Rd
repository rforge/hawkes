\name{Simulate}
\alias{Simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The function simulates a Hawkes process in a one or multi dimensional setting.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Simulate(dim, lambda0, alpha, beta, T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{lambda0}{
%%     ~~Describe \code{lambda0} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dim, lambda0, alpha, beta, T) 
{
    if (!is.integer(dim) || (dim < 1)) {
        stop("'dim' must be a strictly positive integer")
    }
    if (dim == 1) {
        t = .Call("d1_Simulate", lambda0, alpha, beta, T)
    }
    if (dim > 1) {
        t = .Call("dn_Simulate", lambda0, alpha, beta, T)
    }
    t
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

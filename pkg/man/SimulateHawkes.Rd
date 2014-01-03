\name{SimulateHawkes}
\alias{SimulateHawkes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hawkes process simulation Function
}
\description{
The function simulates a Hawkes process for the given parameter, and until a time horizon.
}
\usage{
SimulateHawkes(lambda0, alpha, beta, horizon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda0}{Vector of initial intensity.}
  \item{alpha}{Matrix of excitation.}
  \item{beta}{Vector of Betas. This is a vector even in the multidimensional case, as we consider a unique beta by hawkes process component.}
  \item{horizon}{Time horizon until which the simulation is to be conducted.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector of jump times for every component.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Y. Ogata. (1981)
On Lewis simulation method for point processes.
\emph{IEEE Transactions on Information Theory}, \bold{31}
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
lambda0<-c(0.02,0.02)
alpha<-matrix(c(0.05,0,0,0.05),byrow=TRUE,nrow=2)
beta<-c(0.06,0.06)
horizon<-3600#one hour
h<-SimulateHawkes(lambda0,alpha,beta,horizon)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ hawkes }
%\keyword{ simulation }

\name{A_eff}
\alias{A_eff}
\title{A-efficiency of a binary incomplete block design}
\description{This function computes lower bound to A-efficiency of a binary incomplete block design. Treatment by block incidence matrix of the design is to be supplied as input to the function.}
\usage{A_eff(N)}
\arguments{
	\item{N}{Treatment by block incidence matrix}
	}
\value{
	\item{Aeff}{A-efficiency}
       }
\examples{
N=matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,
1,0,1,1,0,0),nrow=7,byrow=TRUE)
A_eff(N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{A-efficiency}

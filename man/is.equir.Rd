\name{is.equir}
\alias{is.equir}
\title{Equi-replicateness a binary incomplete block design}
\description{This function checks whether an incomplete block design is equi-replicated or not. 
Incidence matrix of the design is to be supplied as input to the function. If the design is 
equir-eplicated, it returns a value of 1 else it returns 0.}
\usage{is.equir(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{equir}{equi-replicated}
       }
\examples{N=matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,
0,0,1,0,0,1,0,1,1,0,0),nrow=7,byrow=TRUE)
is.equir(N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{equi-replicated}

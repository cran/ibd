\name{LIP}
\alias{LIP}
\title{Linear integer programming formulation to update incidence matrix}
\description{This function update an incidence matrix by adding a row through linear integer 
programming formulation for given number of treatments, number of blocks, block sizes and 
concurrence matrix and the previous rows of incidence matrix}
\usage{LIP(v,b,kvec,NNPo,N1,T,rownum,relaxed,binary=TRUE)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{kvec}{block sizes}
	\item{NNPo}{specified concurrence matrix}
	\item{N1}{incidence matrix upto previous rows}
	\item{T}{tabu list of deleted rows}
	\item{relaxed}{whether concurrence constraints are to be relaxed? value 1 for yes, 0 for no.
 Default is 0.}
	\item{binary}{All decision variables are binary to get a binary block. Default is TRUE.}
	}
\value{return the updated incidence matrix}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{internal}
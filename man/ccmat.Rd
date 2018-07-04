\name{ccmat}
\alias{ccmat}
\title{Concurrence matrix for given v, b and k}
\description{This function generates a nearly balanced concurrence matrix from given number
 of treatments (v), number of blocks (b) and block size (k) by trial and error.}
\usage{ccmat(v,b,k)}
\arguments{
	\item{v}{Number of treatments}
	\item{b}{Number of blocks}
	\item{k}{Block size}
	}
\value{
	A v by v matrix is returned if a desired concurrence matrix is found, else it returns a v by v 
matrix of zeros. If the parameters are infeasible for a nearly balanced concurrence matrix, the 
function returns the parameter values.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{internal}

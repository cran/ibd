\name{ccmat}
\alias{ccmat}
\title{Concurrence Matrix for Given v, b and k}
\description{Generates a nearly balanced concurrence matrix from given number
 of treatments (v), number of blocks (b) and block size (k) by trial and error.}
\usage{ccmat(v,b,k)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	}
\value{
	A v by v matrix is returned if a desired concurrence matrix is found, else it returns a v by v 
matrix of zeros. If the parameters are infeasible for a nearly balanced concurrence matrix, the 
function returns the parameter values.}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}

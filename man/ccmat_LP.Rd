\name{ccmat_LP}
\alias{ccmat_LP}
\title{Concurrence matrix for given v, b and k}
\description{This function generates the nearly balanced concurrence matrix from given number 
of treatments (v), number of blocks (b) and block size (k)}
\usage{ccmat_LP(v,b,k)}
\arguments{
	\item{v}{Number of treatments}
	\item{b}{Number of blocks}
	\item{k}{Block size}
	}
\value{
	A v by v matrix is returned if a desired concurrence matrix is found, else it returns 0.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{internal}

\name{ccmat_LP}
\alias{ccmat_LP}
\title{Concurrence Matrix for Given v, b and k}
\description{Generates the nearly balanced concurrence matrix from given number 
of treatments (v), number of blocks (b) and block size (k)}
\usage{ccmat_LP(v,b,k)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	}
\value{
	A v by v matrix is returned if a desired concurrence matrix is found, else it returns 0.}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}

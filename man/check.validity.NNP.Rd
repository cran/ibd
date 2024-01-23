\name{check.validity.NNP}
\alias{check.validity.NNP}
\title{Whether a Given Matrix is Concurrence Matrix or Not}
\description{Checks whether a v by v matrix is concurrence matrix or not. Mainly it
 checks whether the sum of off-diagonal elements of the matrix is (k-1) times the diagonal element for each row of the given matrix. Applicable for proper binary incomplete block design only. If the condition is satisfied, it returns a value of 1 else it returns 0.}
\usage{check.validity.NNP(NNP,k)}
\arguments{
	\item{NNP}{a v by v matrix}
	\item{k}{block size}
	}
\value{
	A value of 1 for valid concurrence matrix or 0 for non valid matrix. 
       }
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}
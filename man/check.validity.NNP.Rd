\name{check.validity.NNP}
\alias{check.validity.NNP}
\title{Whether a given matrix is concurrence matrix or not}
\description{This function checks whether a v by v matrix is concurrence matrix or not. Mainly it
 checks whether the sum of off-diagonal elements of the matrix is (k-1) times the diagonal element for each row of the given matrix. Applicable for proper binary incomplete block design only. If the condition is satisfied, it returns a value of 1 else it returns 0.}
\usage{check.validity.NNP(NNP,k)}
\arguments{
	\item{NNP}{a v by v matrix}
	\item{k}{block size}
	}
\value{
	A value of 1 for valid or 0 for not valid. 
       }
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{internal}
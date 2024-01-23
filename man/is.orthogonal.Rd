\name{is.orthogonal}
\alias{is.orthogonal}
\title{Orthogonality a Block Design}
\description{Checks whether an incomplete block design is orthogonal or not. 
Treatment by block incidence matrix of the design is to be supplied as input to the function. If the design is 
orthogonal, it returns a value of 1 else it returns 0.}
\usage{is.orthogonal(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{orthogonal}{orthogonal}
       }
\examples{N=matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,
0,0,1,0,0,1,0,1,1,0,0),nrow=7,byrow=TRUE)
is.orthogonal(N)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{block design}
\keyword{orthogonal}

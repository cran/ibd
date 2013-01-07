\name{Cmat}
\alias{Cmat}
\title{information matrix from given incidence matrix of a block design and given constant block size}
\description{This function gives the information matrix from a given incidence matrix of a block design and block size}
\usage{Cmat(k,N)}
\arguments{
	\item{k}{block size}
	\item{N}{incidence matrix}
	}
\value{
	\item{C}{v by v information matrix }
       }
\examples{
k=3
N=matrix(c(1	,	0	,	0	,	0	,	1	,	0	,	1	,
0	,	0	,	0	,	1	,	0	,	1	,	1	,
0	,	0	,	1	,	0	,	1	,	1	,	0	,
1	,	0	,	1	,	1	,	0	,	0	,	0	,
0	,	1	,	1	,	0	,	0	,	0	,	1	,
1	,	1	,	0	,	0	,	0	,	1	,	0	,
0	,	1	,	0	,	1	,	1	,	0	,	0	
),nrow=7,byrow=TRUE)
Cmat(k,N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
\keyword{information matrix}

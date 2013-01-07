\name{NNPmat}
\alias{NNPmat}
\title{concurrence matrix from given incidence matrix of a block design}
\description{This function gives the concurrence matrix from a given incidence matrix of a block design}
\usage{NNPmat(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{NNP}{A v by v matrix with diagonal elements as replications of treatments and off-diagonal elements as concurrences between two treatments. Here v is number of treatments.}
       }
\examples{
N=matrix(c(1	,	0	,	0	,	0	,	1	,	0	,	1	,
0	,	0	,	0	,	1	,	0	,	1	,	1	,
0	,	0	,	1	,	0	,	1	,	1	,	0	,
1	,	0	,	1	,	1	,	0	,	0	,	0	,
0	,	1	,	1	,	0	,	0	,	0	,	1	,
1	,	1	,	0	,	0	,	0	,	1	,	0	,
0	,	1	,	0	,	1	,	1	,	0	,	0	
),nrow=7,byrow=TRUE)
NNPmat(N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
\keyword{concurrence matrix}

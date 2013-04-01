\name{C}
\alias{C}
\title{Information matrix from given incidence matrix of a block design}
\description{This function gives the information matrix from a given incidence matrix of a block design}
\usage{C(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{C}{v by v information matrix where v is number of treatments}
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
C(N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
\keyword{information matrix}

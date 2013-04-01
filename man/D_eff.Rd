\name{D_eff}
\alias{D_eff}
\title{D-efficiency of a binary incomplete block design}
\description{This function computes D-efficiency of a binary incomplete block design. Incidence matrix of the design is to be supplied as input to the function.}
\usage{D_eff(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{Deff}{D-efficiency}
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
D_eff(N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{D-efficiency}

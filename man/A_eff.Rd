\name{A_eff}
\alias{A_eff}
\title{A-efficiency of a block design given the parameters v,b,k and incidence matrix}
\description{This function computes A-efficiency of an incomplete block design with given number of treaments(v), number of blocks(b), block size(k) and incidence matrix(N).}
\usage{A_eff(v,b,k,N)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{N}{incidence matrix}
	}
\value{
	\item{Aeff}{A-efficiency}
       }
\examples{
v=7
b=7
k=3
N=matrix(c(1	,	0	,	0	,	0	,	1	,	0	,	1	,
0	,	0	,	0	,	1	,	0	,	1	,	1	,
0	,	0	,	1	,	0	,	1	,	1	,	0	,
1	,	0	,	1	,	1	,	0	,	0	,	0	,
0	,	1	,	1	,	0	,	0	,	0	,	1	,
1	,	1	,	0	,	0	,	0	,	1	,	0	,
0	,	1	,	0	,	1	,	1	,	0	,	0	
),nrow=7,byrow=TRUE)
A_eff(v,b,k,N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{A-efficiency}

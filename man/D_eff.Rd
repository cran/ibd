\name{D_eff}
\alias{D_eff}
\title{D-efficiency of a block design given the parameters v,b,k and incidence matrix}
\description{This function computes D-efficiency of an incomplete block design with given number of treaments(v), number of blocks(b), block size(k) and incidence matrix(N).}
\usage{D_eff(v,b,k,N)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{N}{incidence matrix}
	}
\value{
	\item{Deff}{D-efficiency}
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
D_eff(v,b,k,N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{D-efficiency}

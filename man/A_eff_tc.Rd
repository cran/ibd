\name{A_eff_tc}
\alias{A_eff_tc}
\title{A-efficiency of incomplete block design for test vs control(s) comparisons}
\description{This function computes A-efficiency of incomplete block design for test vs control(s) comparisons}
\usage{A_eff_tc(N,v1,v2,b,k)}
\arguments{	
	\item{N}{incidence matrix}
	\item{v1}{number of test treatments}
	\item{v2}{number of control treatments}
	\item{b}{number of blocks}
	\item{k}{block size}	
	}
\value{
	\item{Aeff}{A-efficiency of the design}
         }
\examples{
N=matrix(c(1	,	1	,	0	,	0	,	0	,
0	,	0	,	0	,	1	,	1	,
0	,	1	,	1	,	0	,	1	,
1	,	0	,	1	,	1	,	0	,
1	,	1	,	1	,	1	,	1	
),nrow=5,byrow=TRUE)
A_eff_tc(N,4,1,5,3)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{test vs control(s) comparison}
\keyword{A-efficiency}
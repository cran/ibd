\name{design_to_N}
\alias{design_to_N}
\title{block design to incidence matrix}
\description{This function generates incidence matrix from a given block design}
\usage{design_to_N(design)}
\arguments{
	\item{design}{design}
	}
\value{
	\item{N}{A treatment by block incidence matrix of order v by b with elements as 0 and 1 where v is number of treatments and b is number of blocks}
       }
\examples{
d=matrix(c(1	,	4	,	6	,
5	,	6	,	7	,
3	,	4	,	5	,
2	,	4	,	7	,
1	,	3	,	7	,
2	,	3	,	6	,
1	,	2	,	5	
),nrow=7,byrow=TRUE)
design_to_N(d)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
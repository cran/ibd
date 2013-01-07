\name{N_to_design}
\alias{N_to_design}
\title{incidence matrix from given block design}
\description{This function generates the block contents from a given incidence matrix}
\usage{N_to_design(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{design}{A matrix with number of rows equal to number of blocks and number of columns equal to block size. Constant block size is assumed. Treatments are numbered as 1, 2, ..., v}
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
N_to_design(N)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}

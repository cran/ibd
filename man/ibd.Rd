\name{ibd}
\alias{ibd}
\title{Binary incomplete block design for given v, b and k and optionally, with a specified concurrence 
matrix}
\description{This function generates an A- and D- efficient binary incomplete block design with given number of 
treaments(v), number of blocks(b) and block size(k) and optionally with a specified concurrence 
matrix(NNP).}
\usage{ibd(v,b,k,NNPo,ntrial,pbar=FALSE)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{NNPo}{optionally, desired concurrence matrix. If not specified, a nearly balanced 
concurrence matrix is obtained automatically.}
	\item{ntrial}{number of trials}
	\item{pbar}{progress bar}
	}
\value{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{NNP}{specified concurrence matrix}
	\item{N}{incidence matrix of the generated design}
	\item{design}{block contents in a b by k matrix}
	\item{conc.mat}{concurrence matrix of the generated design}
	\item{A.efficiency}{lower bound to A-efficiency of the generated design}
	\item{D.efficiency}{lower bound to D-efficiency of the generated design}
	\item{time.taken}{time taken to generate the design}
  }
\examples{
v=9
b=12
k=3
ibd(v,b,k,pbar=FALSE)
}
\references{Mandal, B. N., Gupta, V. K. and Parsad, R. (2013). Application of optimization techniques 
for construction of incomplete block designs. Project report, IASRI, New Delhi.

Mandal, B. N., Gupta, V. K., & Parsad, R. (2014). Efficient Incomplete Block Designs Through 
Linear Integer Programming. American Journal of Mathematical and Management Sciences, 33(2),
 110-124.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{incidence matrix}
\keyword{concurrence matrix}
\keyword{A-efficiency}
\keyword{D-efficiency}
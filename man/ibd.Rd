\name{ibd}
\alias{ibd}
\title{Binary Incomplete Block Design for Given v, b and k and Optionally, with a Specified Concurrence 
Matrix}
\description{Generates an A- and D- efficient binary incomplete block design with given number of 
treaments(v), number of blocks(b) and block size(k) and optionally with a specified concurrence 
matrix(NNP).}
\usage{ibd(v,b,k,NNPo,ntrial=5,pbar=FALSE)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{NNPo}{optionally, desired concurrence matrix. If not specified, a nearly balanced 
concurrence matrix is obtained automatically.}
	\item{ntrial}{number of trials. Default is 5.}
	\item{pbar}{progress bar. Default is FALSE.}
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
ibd(v = 7,b = 7,k = 4, pbar=FALSE)
}
\references{Mandal, B. N., Gupta, V. K. and Parsad, R. (2013). Application of optimization techniques 
for construction of incomplete block designs. Project report, IASRI, New Delhi.

Mandal, B. N., Gupta, V. K., & Parsad, R. (2014). Efficient Incomplete Block Designs Through 
Linear Integer Programming. American Journal of Mathematical and Management Sciences, 33(2),
 110-124.}
 
\note{This function works best for values of number of treatments (v) up to 30 and block size (k) up to 10. However, for block size (k) up to 3, much larger values of number of treatments (v) may be used.}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{incidence matrix}
\keyword{concurrence matrix}
\keyword{A-efficiency}
\keyword{D-efficiency}
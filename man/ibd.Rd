\name{ibd}
\alias{ibd}
\title{Incomplete block design for given v,b and k and optionally, with a specified concurrence matrix}
\description{This function generates an efficient incomplete block design with given number of treaments(v), number of blocks(b) and block size(k) and optionally with a specified concurrence matrix(NNP).}
\usage{ibd(v,b,k,ntrial,NNPo,pbar=FALSE)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{ntrial}{number of trials}
	\item{NNPo}{desired concurrence matrix}
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
	\item{A.efficiency}{A-efficiency of the generated design}
	\item{D.efficiency}{D-efficiency of the generated design}
	\item{time.taken}{time taken to generate the design}
       }
\examples{
v=9
b=12
k=3
ibd(v,b,k,pbar=TRUE)
}
\references{Mandal, BN, Gupta, VK and Parsad, R. (2013). Application of optimization techniques for construction of incomplete block designs. Project report, IASRI, New Delhi.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{incidence matrix}
\keyword{concurrence matrix}
\keyword{A-efficiency}
\keyword{D-efficiency}
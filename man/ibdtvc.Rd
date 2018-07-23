\name{ibdtvc}
\alias{ibdtvc}
\title{incomplete block design for test vs control(s) comparions}
\description{This function generates an incomplete block design for test vs control(s) comparisons 
with specified parameters and concurrence matrix.}
\usage{ibdtvc(v1,v2,b,k,NNPo,ntrial=5,pbar=FALSE)}
\arguments{	
	\item{v1}{number of test treatments}
	\item{v2}{number of control treatments}
	\item{b}{number of blocks}
	\item{k}{block size}		
	\item{NNPo}{desired concurrence matrix}
	\item{ntrial}{number of trials, default is 5}
	\item{pbar}{Logical value indicating whether progress bar will be displayed. Default
 is FALSE.}
	}
\value{
	v1=v1,v2=v2,b=b,k=k,design=design,N=N, NNP=NNP,Aeff=Aeff)
	\item{v1}{number of test treatments}
	\item{v2}{number of control treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{design}{generated block design}
	\item{N}{treatment by block incidence matrix of the generated block design}
	\item{NNP}{concurrence matrix of the generated design}
	}
\examples{NNPo=matrix(c(7,3,3,3,3,3,3,3,3,7,3,3,3,3,3,3,3,3,7,3,3,3,3,3,3,3,3,7,3,3,3,3,3,3,3,3,7,
3,3,3,3,3,3,3,3,7,3,3,3,3,3,3,3,3,9,9,3,3,3,3,3,3,9,9),nrow=8,byrow=TRUE)
ibdtvc(6,2,15,4,NNPo)
}
\references{Mandal, B. N., Gupta, V.K. and Parsad, R. (2013). Binary Incomplete Block Designs with a 
Specified Concurrence Matrix through Integer Programming, to be submitted for publication}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{test vs control(s) comparison}
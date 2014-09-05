\name{bibd}
\alias{bibd}
\title{Balanced incomplete block design for given parameters}
\description{This function generates a balanced incomplete block design with given number of
 treaments(v), number of blocks(b), number of replications (r), block size(k) and number of 
concurrences (lambda).}
\usage{bibd(v,b,r,k,lambda,ntrial,pbar=FALSE)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{r}{number of replications}
	\item{k}{block size}
	\item{lambda}{number of concurrences}
	\item{ntrial}{number of trials}
	\item{pbar}{logical value indicating whether progress bar will be displayed or not. Default is 
FALSE}
	}
\value{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{r}{number of replications}
	\item{k}{block size}
	\item{lambda}{number of replications}
	\item{design}{block contents in a b by k matrix}
	\item{N}{incidence matrix of the generated design}
	\item{NNP}{concurrence matrix of the generated design}
	\item{Aeff}{Lower bound to the A-efficiency of the generated design}
	\item{Deff}{Lower bound to the D-efficiency of the generated design}
       }
\examples{bibd(7,7,3,3,1,pbar=FALSE)}
\references{
Mandal, BN, Gupta, VK and Parsad, R. (2013). Application of optimization techniques for 
construction of incomplete block designs. Project report, IASRI, New Delhi.
Mandal, B. N., Gupta, V. K., & Parsad, R. (2014). Efficient Incomplete Block Designs Through
 Linear Integer Programming. American Journal of Mathematical and Management Sciences, 33(2), 
110-124.
Mandal, B. N. (2014). Linear integer programming approach to construction of balanced incomplete
 block designs. Communications in Statistics-Simulation and Computation, (in print).
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{balanced incomplete block design}
\keyword{incidence matrix}
\keyword{concurrence matrix}
\keyword{A-efficiency}
\keyword{D-efficiency}
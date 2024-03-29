\name{bibd}
\alias{bibd}
\title{Balanced Incomplete Block Design for Given Parameters}
\description{Generates a balanced incomplete block design with given number of
 treaments (v), number of blocks (b), number of replications (r), block size (k) and number of 
concurrences (lambda).}
\usage{bibd(v,b,r,k,lambda,ntrial=5,pbar=FALSE)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{r}{number of replications}
	\item{k}{block size}
	\item{lambda}{number of concurrences}
	\item{ntrial}{number of trials. Default value is 5.}
	\item{pbar}{logical value indicating whether progress bar will be displayed or not. Default is 
FALSE}
	}
\value{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{r}{number of replications}
	\item{k}{block size}
	\item{lambda}{number of concurrences}
	\item{design}{block contents in a b by k matrix}
	\item{N}{treatments by blocks incidence matrix of the generated design}
	\item{NNP}{concurrence matrix of the generated design}
	\item{Aeff}{Lower bound to the A-efficiency of the generated design}
	\item{Deff}{Lower bound to the D-efficiency of the generated design}
       }
\examples{
bibd(7,7,3,3,1)
bibd(9,12,4,3,1)}

\references{
Mandal, B. N., Gupta, V. K. and Parsad, R. (2013). Application of optimization techniques for 
construction of incomplete block designs. Project report, IASRI, New Delhi.

Mandal, B. N., Gupta, V. K., & Parsad, R. (2014). Efficient Incomplete Block Designs Through
 Linear Integer Programming. American Journal of Mathematical and Management Sciences, 33(2), 
110-124.

Mandal, B. N. (2015). Linear integer programming approach to construction of balanced incomplete
 block designs. Communications in Statistics-Simulation and Computation, 44:6, 1405-1411.
}

\note{The function works best for values of number of treatments (v) up to 30 and block size (k) up to 10. However, for block size (k) up to 3, much larger values of number of treatments (v) may be used.}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{balanced incomplete block design}
\keyword{incidence matrix}
\keyword{concurrence matrix}
\keyword{A-efficiency}
\keyword{D-efficiency}
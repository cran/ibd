\name{btib}
\alias{btib}
\title{balanced treatment incomplete block designs}
\description{This function generates a balanced treatment incomplete block design for specified
 parameters.}
\usage{btib(v,b,r,r0,k,lambda,lambda0,ntrial=5,pbar=FALSE)}
\arguments{	
	\item{v}{number of test treatments}	
	\item{b}{number of blocks}
	\item{r}{number of replications of test treatments}
	\item{r0}{number of replications of the control treatment}
	\item{k}{block size}	
	\item{lambda}{number of concurrences among test treatments}
	\item{lambda0}{number of concurrences between test treatments and the control treatment}
	\item{ntrial}{number of trials. Default is 5.}
	\item{pbar}{Logical value indicating whether progress bar will be displayed or not. Default is
 FALSE.}
	}
\value{
	\item{v}{number of test treatments}
	\item{b}{number of blocks}
	\item{r}{number of replications of test treatments}
	\item{r0}{number of replications of the control treatment}
	\item{k}{block size}
	\item{lambda}{number of concurrences among test treatments}
	\item{lambda0}{number of concurrences between test treatments and the control treatment}
	\item{design}{generated block design}
	\item{N}{treatment by block incidence matrix of the generated block design}
	\item{NNP}{concurrence matrix of the generated design}
	\item{Aeff}{A-efficiency of the generated design}
  }
\references{Mandal, B. N., Gupta, V. K. and Parsad, R. (2013). Application of optimization techniques 
for construction of incomplete block designs. Project report, IASRI, New Delhi.

Mandal, B. N., Gupta, V. K., & Parsad, R. (2014). Balanced treatment
incomplete block designs through integer programming, Communications in Statistics - Theory 
and Methods, 46:8, 3728-3737.
}
\note{The function works best for values of number of treatments (v) up to 30 and block size (k) up to 10. However, for block size (k) up to 3, much larger values of number of treatments (v) may be used.}
\examples{btib(4,6,3,6,3,1,3,10)}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{test vs control(s) comparison}
\keyword{A-efficiency}
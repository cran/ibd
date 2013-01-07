\name{btib}
\alias{btib}
\title{balanced treatment incomplete block design}
\description{This function generates balanced treatment incomplete block design for specified parameters.}
\usage{btib(v,b,r,r0,k,lambda,lambda0,ntrial)}
\arguments{	
	\item{v}{number of test treatments}	
	\item{b}{number of blocks}
	\item{r}{number of replications of test treatments}
	\item{r0}{number of replications of control treatment}
	\item{k}{block size}	
	\item{lambda}{number of concurrences among test treatments}
	\item{lambda0}{number of concurrences between test treatments and control treatment}
	\item{ntrial}{number of trials}
	}
\value{
	\item{v}{number of test treatments}
	\item{b}{number of blocks}
	\item{r}{number of replications of test treatments}
	\item{r0}{number of replications of control treatment}
	\item{k}{block size}
	\item{lambda}{number of concurrences among test treatments}
	\item{lambda0}{number of concurrences between test treatments and control treatment}
	\item{design}{generated block design}
	\item{N}{incidence matrix of the generated block design}
	\item{NNP}{concurrence matrix of the generated design}
	\item{Aeff}{A-efficiency of the generated design}
         }
\examples{btib(4,6,3,6,3,1,3,10)}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{test vs control(s) comparison}
\keyword{A-efficiency}
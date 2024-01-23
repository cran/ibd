\name{btib1}
\alias{btib1}
\title{Balanced Treatment Incomplete Block Designs}
\description{Generates a balanced treatment incomplete block design for specified
 parameters by searching all possible combinations.}
\usage{btib1(v,b,r,r0,k,lambda,lambda0)}
\arguments{	
	\item{v}{number of test treatments}	
	\item{b}{number of blocks}
	\item{r}{number of replications of test treatments}
	\item{r0}{number of replications of the control treatment}
	\item{k}{block size}	
	\item{lambda}{number of concurrences among test treatments}
	\item{lambda0}{number of concurrences between test treatments and control treatment}
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
	\item{N}{treatment by block incidence matrix of the generated block design}
	\item{NNP}{concurrence matrix of the generated design}
	\item{Aeff}{A-efficiency of the generated design}
         }
\references{Mandal, B. N., Gupta, V. K. and Parsad, R. (2013). Application of optimization techniques 
for construction of incomplete block designs. Project report, IASRI, New Delhi.

MANDAL, B. N., GUPTA, V. K. and PARSAD, R. (2012). Generation of Binary Incomplete Block Design with
 a Specified Concurrence Matrix. Journal of Statistics & Applications, 7.
}
\note{The function works best for values of number of treatments (v) up to 30 and block size (k) up to 10. However, for block size (k) up to 3, much larger values of number of treatments (v) may be used.}
\examples{btib(4,6,3,6,3,1,3)}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{test vs control(s) comparison}
\keyword{A-efficiency}
\name{design_to_N}
\alias{design_to_N}
\alias{N}
\title{Block Design to Treatment by Block Incidence Matrix}
\description{Generates treatment by block incidence matrix from a given block design}
\usage{design_to_N(design)
N(design)}
\arguments{
	\item{design}{design}
	}
\value{
	\item{N}{A treatment by block incidence matrix of order v by b with elements as 0 and 1 where
 v is the number of treatments and b is the number of blocks}
       }
\examples{design = matrix(c(1,4,6,5,6,7,3,4,5,2,4,7,1,3,7,2,3,6,1,2,5),nrow=7,byrow=TRUE)
design_to_N(design)
# or alternatively
N(design)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
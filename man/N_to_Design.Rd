\name{N_to_design}
\alias{N_to_design}
\alias{design}
\title{Block Design from Given Treatment by Block Incidence Matrix }
\description{Generates the block contents from a given treatment by block incidence matrix}
\usage{N_to_design(N)
design(N)}
\arguments{
	\item{N}{treatment by block incidence matrix}
	}
\value{
	\item{design}{a matrix with number of rows equal to number of blocks and number of columns 
equal to block size. Constant block size is assumed. Treatments are labelled as 1, 2, ..., v.}
       }
\examples{N=matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,
1,1,0,0,0,1,0,0,1,0,1,1,0,0),nrow=7,byrow=TRUE)
N_to_design(N)
design(N)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}

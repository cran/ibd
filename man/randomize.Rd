\name{randomize}
\alias{randomize}
\title{Randomize a block design}
\description{Randomize a given block design}
\usage{randomize(design)}
\arguments{
	\item{design}{design}
	}
\value{
	\item{design}{Block design with a constant block size}
       }
\examples{design = matrix(c(1,4,6,5,6,7,3,4,5,2,4,7,1,3,7,2,3,6,1,2,5),nrow=7,byrow=TRUE)
randomize(design)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
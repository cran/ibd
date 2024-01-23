\name{is.connected}
\alias{is.connected}
\title{Connctedness of a Binary Incomplete Block Design}
\description{Checks whether an incomplete block design is connected or not. Treatment by block incidence 
matrix of the design is to be supplied as input to the function. If the design is connected, it 
returns a value of 1 else it returns 0.}
\usage{is.connected(N)}
\arguments{
	\item{N}{incidence matrix}
	}
\value{
	\item{connected}{connctedness}
       }
\examples{N=matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,
0,0,1,0,0,1,0,1,1,0,0),nrow=7,byrow=TRUE)
is.connected(N)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{Connectedness}

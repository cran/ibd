\name{Cmatrix}
\alias{Cmatrix}
\alias{Information.Matrix}
\title{Information Matrix of a Block Design}
\description{Gives the information matrix from a given treatment by block incidence matrix of a block
 design}
\usage{Cmatrix(N)}
\arguments{
	\item{N}{treatment by block incidence matrix}
	}
\value{
	\item{Cmatrix}{v by v information matrix where v is number of treatments}
       }
\examples{
N = matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0
,1,0,1,1,0,0),nrow=7,byrow=TRUE)
Cmatrix(N)
Information.Matrix(N)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incidence matrix}
\keyword{incomplete block design}
\keyword{information matrix}

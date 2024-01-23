\name{D_eff}
\alias{D_eff}
\title{D-efficiency of a Binary Incomplete Block Design}
\description{Computes lower bound to D-efficiency of a binary incomplete block design}
\usage{D_eff(N)}
\arguments{
	\item{N}{treatment by block incidence matrix}
	}
\value{
	\item{Deff}{lower bound to D-efficiency}
       }
\examples{
N=matrix(c(1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0
,1,0,1,1,0,0),nrow=7,byrow=TRUE)
D_eff(N)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{D-efficiency}

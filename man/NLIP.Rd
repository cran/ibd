\name{NLIP}
\alias{NLIP}
\title{An alternative approach to obtain a block design with specified concurrence matrix}
\description{This function generates a block design with specified parameters and concurrence matrix through an alternative linear integer programming approach}
\usage{NLIP(v,b,k,NNPo)}
\arguments{
	\item{v}{number of test treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{NNPo}{specified concurrence matrix}
	}
\value{Returns a design or a text saying that no design found}
\references{MANDAL, B., GUPTA, V., & PARSAD, R. (2012). Generation of Binary Incomplete Block Design with a Specified Concurrence Matrix. Journal of Statistics & Applications, 7.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{internal}
\name{btibgen}
\alias{btibgen}
\title{Generate a treatment by block incidence matrix of a BTIB design for given v,b,k, and concurrence matrix}
\description{This function generate an incidence matrix of a BTIB design for given v,b,k and concurrence matrix}
\usage{ibdgen(v,b,k,NNPo,ntrial,pbar)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{NNPo}{specified concurrence matrix}
	\item{ntrial}{number of trials, default is 5}
	\item{pbar}{logical, progress bar to be displayed. Default is FALSE}
	}
\value{returns a v by b incidence matrix or a text that design was not found}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{internal}
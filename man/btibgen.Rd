\name{btibgen}
\alias{btibgen}
\title{Treatment by Block Incidence Matrix of a BTIB Design}
\description{This function generate an incidence matrix of a BTIB design for given v,b,k and concurrence matrix}
\usage{btibgen(v,b,k,NNPo,ntrial,pbar)}
\arguments{
	\item{v}{number of treatments}
	\item{b}{number of blocks}
	\item{k}{block size}
	\item{NNPo}{specified concurrence matrix}
	\item{ntrial}{number of trials, default is 5}
	\item{pbar}{logical, progress bar to be displayed. Default is FALSE}
	}
\value{returns a v by b incidence matrix or a text that design was not found}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}
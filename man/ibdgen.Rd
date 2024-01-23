\name{ibdgen}
\alias{ibdgen}
\title{Generate a Treatment by Block Incidence Matrix for Given v,b,k, and Concurrence Matrix}
\description{Generates an incidence matrix for given v,b,k and concurrence matrix}
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
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}
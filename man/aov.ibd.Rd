\name{aov.ibd}
\alias{aov.ibd}
\title{Analysis of variance of data from a block design}
\description{This function performs analysis of variance of data from experiments using a block design. Data should be specified with 3 columns, first column should have treatments, second column should have blocks and third column should have response variable. Treatments and blocks should be numbered as 1, 2,... There is also option for estimating and testing a treatment contrast.}
\usage{aov.ibd(data,contrast)}
\arguments{
	\item{data}{data}
	\item{contrast}{treatment contrast}
	}
\value{
	\item{aov.table}{ANOVA Table}
	\item{contrast.table}{Estimate and test of significance of the contrast}
       }
\references{Dey,A. (1986). Theory of block designs.Wiley Eastern Limited, New Delhi.}
\examples{
trt=c(1,2,3,7,9,8,6,4,5,1,7,4,8,5,2,3,6,9,1,8,6,4,2,9,3,5,7,1,9,5,2,6,7,8,3,4)
blk=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12)
y=c(77,65,75,47,61,60,54,60,65,70,62,62,72,55,55,50,40,60,63,67,54,62,53,57,68,67,66,69,62,52,61,63,79,65,65,38)
data=data.frame(trt,blk,y)
contrast=c(1,-1,0,0,0,0,0,0,0)
aov.ibd(data,contrast)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{block design}
\keyword{analysis of variance}
\keyword{treatment contrast}

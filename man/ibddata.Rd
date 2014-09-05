\name{ibddata}
\alias{ibddata}
\docType{data}
\title{
Data from an experiment using incomplete block design
}
\description{
Data from an experiment using incomplete block design
}
\usage{data("ibddata")}
\format{
  A data frame with 36 observations on the following 3 variables.
  \describe{
    \item{\code{trt}}{
Treatments}
    \item{\code{blk}}{Blocks}
    \item{\code{y}}{The response variable}
  }
}
\details{
The experiment used a balanced incomplete block design.
}
\references{Dey,A. (1986). Theory of block designs.Wiley Eastern Limited, New Delhi.}
\examples{
data(ibddata)
}
\keyword{datasets}

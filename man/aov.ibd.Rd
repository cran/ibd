\name{aov.ibd}
\alias{aov.ibd}
\title{Analysis of variance, least square means and contrast analysis of data from a block design}
\description{This function performs intrablock analysis of variance of data from experiments using
 a block design. It also computes least square means of the factor variables (e.g. treatments) 
and optionally estimates and tests the contrasts of factor variables (e.g treatments). }
\usage{aov.ibd(formula,specs,data,contrast,joint=FALSE,details=FALSE,sort=TRUE,by=NULL,
alpha=0.05,Letters = "ABCDEFGHIJ",...)}
\arguments{
	\item{formula}{A formula specifying the model of the form response~treatment+block or 
	response~block+treatment. Make sure the treatment and blocks are factor variables.}
	\item{specs}{A character vector specifying the names of the factors over which Least 
square means are desired}
	\item{data}{A data frame in which the variables specified in the formula will be found.
 If missing, the variables are searched for in the standard way.}
	\item{contrast}{A matrix whose rows are contrasts of factors (e.g. treatments)}
	\item{joint}{If contrast argument has more than one row, then whether a joint test of 
the contrasts will be performed. Default is FALSE. If joint=TRUE, a check is performed
 whether the contrasts are pairwise orthogonal or not and then if orthoghonal, joint 
test is performed.}
	\item{details}{Logical, if details=TRUE then all objects including lm object from lm(),
 lsm object from lsmeans() are returned. Default is FALSE.}
	\item{sort}{Logical value determining whether the least square means are sorted before the
 comparisons are produced. Default is TRUE.}
	\item{by}{Character value giving the name or names of variables by which separate families of 
comparisons are tested. If NULL, all means are compared.}
	\item{alpha}{Numeric value giving the significance level for the comparisons}
	\item{Letters}{Characters to be used for compact letter display of groups of factor variables
 over which least square means are computed. Default is english alphabet capital letters "ABCDEFGHIJ"}
	\item{...}{Not used}
} 
\value{Returns a list with following components
	\item{lm.obj}{An object of class lm if details=TRUE}
	\item{ANOVA.table}{ANOVA table from the fitted lm object}
	\item{LSMEANS}{Least square means with compact letter display}
	\item{contrast.analysis}{Contrast analysis result if contrast matrix was supplied}
}
\details{The function makes use of lm() function in R and Anova() function in car package with 
specification of Type III sum of squares and lsmeans(), cld(), contrast() functions in lsmeans()
 package and combines the results in a single place.}
\examples{
data(ibddata)
aov.ibd(y~factor(trt)+factor(blk),data=ibddata)
contrast=matrix(c(1,-1,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0),nrow=2,byrow=TRUE)
aov.ibd(y~factor(trt)+factor(blk),specs="trt",data=ibddata,contrast=contrast)
}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{analysis of variance}
\keyword{contrast analysis}
\keyword{balanced incomplete block design}
\keyword{Least square means}
\keyword{LSmeans}
\keyword{Type III sum of squares}
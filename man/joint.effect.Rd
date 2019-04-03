\name{joint.effect}
\alias{joint.effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of the joint effect of a list of mediators
}
\description{
Compute the estimations, standard deviations and confidence intervals of the joint effect of the list of variables identified by vari. 
}
\usage{
joint.effect(object,vari,alpha=0.05) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
a mma.big object created initially call to mma.big.
}
 \item{vari}{
a list of mediators whose joint effect is to be estimated. vari can alsobe the name of a multicategorical mediator.
}
  \item{alpha}{
the significance level for the confidence interval.
}
}
\details{
joint.effect gives a list of the estimations and summary statistics for the joint mediation effect of a list of mediators based on the bootstrap results.
}
\value{
The function joint.effect return a list. The first item list all the included variables. The second item, effect, is the list for the joint mediation effects. The third item, relative.effect, is the list of relative effects. Under them, est is the estimation of the corresponding (relative) joint mediation effect based on the whole data; mean is the average estimated (relative) effects from the bootstrap samples; and sd is the standard deviation of the estimate from the bootstrap sample. upbd and lwbd are the upper and lower bounds of the confidence interval of the estimation using parametric method from the bootstrap sample, upbd_q and lwbd_q are the corresponding quantiles of the estimation from the bootstrap sample.
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23. 

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11. 

Li, B. and Yu, Q. (2018). "Mediation Analysis with Large Data Sets". 
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}
}

\seealso{
\code{"\link[=mma.big]{mma.big}"} 
}
\examples{
#See help(mma.big)
}

%\keyword{ Summary }

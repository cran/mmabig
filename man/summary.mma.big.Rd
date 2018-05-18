\name{summary.mma.big}
\alias{summary.mma.big}
\alias{print.summary.mma.big}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of an mma.big project
}
\description{
Compute the estimations, standard deviations and confidence intervals of the mediation effects.
}
\usage{
\method{summary}{mma.big}(object,...,alpha=0.05,plot=TRUE,RE=FALSE,quant=T,ball.use=T) 
\method{print}{summary.mma.big}(x,...,digit=3) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
a mma.big object created initially call to mma.big.
}
 \item{x}{
a summary.mma.big object created initially call to summary.mma.big.
}
  \item{\dots}{
other arguments passed to the print function.
}
  \item{alpha}{
the alpha level for confidence interval.}
  \item{plot}{
default is TRUE, if ture, draw a barplot of the mediation effects with confidence intervals.}
  \item{RE}{
default is FALSE, if ture, show the inferences on relative mediation effects.}
  \item{quant}{
default is TRUE, if ture and ball.use is F, draw the confidence intervals of relative effects using quantile.}
  \item{ball.use}{
default is TRUE, if ture, draw the confidence intervals of relative effects using the confidence ball. If both quant and ball.use are false, draw the confidence intervals based on the standard deviaitons from bootstrap estimates.}
  \item{digit}{
the number of decimal digits to keep.
}
}
\details{
summary.mma.big gives a list of the estimations and summary statistics based on the bootstrap results. If plot=T, draw a barplot of the estimated direct effect of the predictor and indirect effects of the mediators. Relative effects is defined as the (in)direct effect divided by the total effect. The plot is arranged in order from the largest to the smallest mediation effects.
}
\value{
The function summary.mma.big return a list. The first item, result, is the list for mediation effects, and the second item, re, is the list for relative effects. Under them, est is the estimation of the corresponding (relative) mediation effects based on the whole data; mean is the average estimated (relative) effects from the bootstrap samples; and sd is the standard deviation of the estimates from the bootstrap sample. upbd and lwbd are the upper and lower bound of the confidence interval of the estimation using parametric method from the bootstrap sample, upbd_q and lwbd_q are the corresponding quantiles of the estimation from the bootstrap sample; upbd_b and lwbd_b are the corresponding ball confidence interval of the estimation from the bootstrap sample.
If plot=T, draw a barplot of the mediation effects with confidence intervals.
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

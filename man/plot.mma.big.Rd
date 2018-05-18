\name{plot.mma.big}
\alias{plot.mma.big}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relative effects plot of the fitted mma.big object
}
\description{
Plot the marginal effect of the selected variable(mediator) on the outcome, and the marginal effect of the predictor on the selected variable.
}
\usage{
\method{plot}{mma.big}(x,vari,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a mma.big object created initially call to mma.big.
}
  \item{vari}{
  an indices or the name of the mediator to plot.
  }
  \item{\dots}{
other arguments passed to the plot function.
}
}
\details{
plot.mma.big plots on the upper panel, the boxplot of fitted coefficients of the selected mediator on each outcome from bootstrap samples. The lower panel is for each predictor, the marginal change in the mediator when the predictor changes. If the predictors are binary, show the boxplot of the mean changes in mediator when the corresponding predictor changes from 0 to 1. If the predictors are continuous, show the smoothed changes of the mean mediator with the predictor, with 95\% confidence bounds.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{"\link[=mma.big]{mma.big}"}
}
\examples{
#See help(mma.big)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Plot }

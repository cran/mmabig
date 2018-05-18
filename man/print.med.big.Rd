\name{print.med.big}
\alias{print.med.big}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print an med.big object
}
\description{
Print the estimation of mediation effects from an med.big object: from the function med.big.
}
\usage{
\method{print}{med.big}(x,...,pred.new=NULL,w.new=NULL,print=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a med object created initially call to med.
}
  \item{\dots}{
other arguments passed to the print function.
}
  \item{pred.new}{
the new set of predictors.
}
  \item{w.new}{
the new weight with the new set of predictors.
}
  \item{print}{
print the results if TRUE.
}
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
\code{"\link[=med.big]{med.big}"}
}
\examples{
#See help(mma.big)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ Print }

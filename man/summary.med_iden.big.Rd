\name{summary.med_iden.big}
\alias{summary.med_iden.big}
\alias{print.summary.med_iden.big}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for class "med_iden.big".
}
\description{
Compute the estimations, standard deviations and confidence intervals of the mediation effects.
}
\usage{
\method{summary}{med_iden.big}(object,...,only=F) 
\method{print}{summary.med_iden.big}(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{
a med_iden object created initially call to data.org.big.
}
 \item{x}{
a summary.med_iden.big object created initially call to summary.med_iden.big
}
 \item{\dots}{
other arguments passed to the print function.
}
 \item{only}{
if only=T, show test results for selected covariates and mediators only.
}
}
\details{
summary.med_iden gives a list of identified mediators, covariates and the test results. 
}
\value{
The function summary.med_iden return a list of covariates and mediators as identified by a series tests.
  \item{mediator }{variable names of the identified mediators, either siginicant in both full model and in relate to the predictor, or being a member of the pre-identified joint mediators.}
  \item{covariate }{variable names of covariates: being significant in the full model but not significantly relate with the predictor.}
  \item{tests }{statistical test results.}
  \item{results }{the original object.}
}
\references{
Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}
}

\seealso{
\code{"\link[=mma.big]{mma.big}"} 
}
\examples{
#see help(mma.big)
}

%\keyword{ Summary }

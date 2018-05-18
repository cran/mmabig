\name{data.org.big}
\alias{data.org.big}
\title{
Data Organization and Identify Potential Mediators
}
\description{
Do a preliminary data analysis to identify potential mediators and covariates. Each variable listed in jointm is forced in the final estimation model as a mediator. Also organize the data into a format that can be directly used for the mediation analysis functions.
}
\usage{
data.org.big(x,y,pred,mediator=NULL,contmed=NULL,binmed=NULL,binref=NULL,catmed=NULL,
                       catref=NULL,jointm=NULL, 
                       family1=as.list(rep(NA,ncol(data.frame(y)))),
                       predref=NULL,alpha=1,alpha1=0.01,alpha2=0.01,testtype=2, w=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates
}
  \item{y}{
the vector of outcome variable.The outcome can be binary, continuous, multicategorical or of "Surv" class (see survival package for help). For a factor, the last level in alphabetical order is the target class.
}
  \item{pred}{
the column or matrix of predictor(s): the predictor is the exposure variable, it can be a binary or multi-categorical factor or one/a few contiuous variable(s).
}
  \item{mediator}{
the list of mediators (column numbers in x or by variable names). The mediators to be checked can be identified by "contmed", "binmed" and "catmed", or by this argument, "mediator", where binary and categorical mediators in x are identified as factors or characters, the reference group is the first level of the factor or factorized character. if a mediator has only two unique values, the mediator is identified as binary. If the reference groups need to be changed, the binary or categorical mediators can be listed in binmed or catmed, and the corresponding reference group in binref or catref.
}
  \item{contmed}{
a vector of variable names or column numbers that locate the potential continuous mediators in x.
}
  \item{binmed}{
a vector of column numbers that locate the potential binary mediators in x.
}
  \item{binref}{
the defined reference groups of the binary potential mediators in binmed. The first levels of the mediators if is null.
}
  \item{catmed}{
a vector of variable names or column numbers that locate the potential categorical mediators in x. The first levels of the mediators if is null.
}
  \item{catref}{
the defined reference groups of the categorical potential mediators in catmed.
}
  \item{jointm}{
a vector of variables that need to be forced in as mediators.
}
%  \item{refy}{
%if y is binary, the reference group of y. The default is the last level of as.factor(y). 
%}
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component in generalized linear model.  The default value of family1 is binomial(link = "logit") for binary y, and gaussian(link="identity") for continuous y.
}
  \item{predref}{
if the predictor is binary, identify the reference group of the binary predictor. The default is the first level of the predictor. The value of the predictor is 0 for the reference grop.
}
  \item{alpha}{
the elasticnet mixing parameter such that 0<=alpha<=1, with alpha=1 is the lasso penalty, and alpha=0 be the ridge penalty. By default, alpha=1.
}
  \item{alpha1}{
the significance level at which to test if the potential mediators (identified by contmed, binmed, and catmed) can be used as a covariate or mediator in estimating y when all variables in x are included in the model.  The default value is alpha1=0.01.
}
  \item{alpha2}{
the significant level at which to test if a potential mediator is related with the predictor. The default value is alpha2=0.01.
}  
  \item{testtype}{if the testtype is 1, covariates/mediators are identified using full model; if the testtype is 2 (by default), covariates/mediators are tested one by one in models with the predictor only.
}
  \item{w}{the weight for data analysis, by default is rep(1,length(y)).
}
}

\value{
data.org.big returns a "med_iden" class (see mma package), a list with the organized data and identifiers of the potential mediators in the organized data set.
  \item{x }{the organized data frame that include all potential mediators and covariates that should be used to estimate the outcome.}
  \item{dirx }{the vector/matrix of predictor(s)/exposure variable(s).}
  \item{contm }{the column numbers of x that locate the potential continuous mediators.}
  \item{binm }{when the predictor is continuous, binm gives the column numbers of x that locate the potential binary mediators.}
  \item{catm }{when the predictor is binary, catm gives the column numbers of x that locate the potential binary or categorical mediators; when the predictor is continuous, catm gives a list where the first item is the number of potential categorical mediators, and the following items give the column numbers of each binarized categorical mediator in x.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators, and each of the following items identifies the column numbers of the mediators in the newly organized x for each group of joint mediators.}
  \item{y }{the vector/matrix of outcome(s).}
  \item{y_type }{the variable type of outcome(s): 1 is continuous, 2 is binary, 3 is reserved for multi-categorical (no 3 would show in y_type, since all categorical responses are binarized), and 4 is survival.}
  \item{fullmodel}{a list with each item the full linear model fitted with all potential mediators and covariates for each response.}
  \item{rela}{p-values of tests on the realtionship between the predictor(s) and each potiential mediator.}
  \item{P1}{If testtype=2, P1 gives the p-value of the corresponding variables in predicting the outcome(s) when only the variable and predictor are covariates in the model.}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi:10.1037/0022-3514.51.6.1173>. The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations. J. Pers Soc Psychol, 51(6), 1173-1182.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23. 

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11. 

Li, B. and Yu, Q. (2018). "Mediation Analysis with Large Data Sets". 
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}
}
\note{
All other variables in x but not identified by mediator, contmed, binmed, or catmed are forced in the final model as covariates.  Compared with data.org, joint mediators are considered in this function. Every variable in the jointm should be listed in contmed, binmed, or catmed, and these variables are forced to be included as mediators for further mediation analysis. A variables can be included in more than one groups of joint mediators in jointm.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#See help(mma.big)
}
\keyword{ Mediator Tests }

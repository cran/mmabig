\name{med.big}
\alias{med.big}
\title{
Mediation Analysis with Binary or Continuous Predictor
}
\description{
To estimate the mediation effects when the predictor is binary or continuous.
}
\usage{
med.big(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm,  
                  contm = data$contm, catm = data$catm, jointm = data$jointm, 
                  allm = c(contm, binm, unique(unlist(catm)[-1])),margin=1,
                  df=1,family1=data$family1,refy=rep(NA,ncol(y)), 
                  binpred=data$binpred,type=NULL,w=NULL,alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
the list of result from data.org that organize the covariates, mediators, predictor and outcome. If data is FALSE, then need to set x1, y1, dirx, contm, catm, and jointm.
}
   \item{x}{
a data frame contains all mediators and covariates.  Need to set up only when data is FALSE.
}
  \item{y}{
the vector of outcome variable.  Need to set up only when data is FALSE.
}
  \item{dirx }{the vector or matrix of predictor(s). The reference group is set to be 0. Need to set up only when data is FALSE.}
  \item{binm}{
the variable names or the column number of x that locates the binary mediators. Need to set up only when data is FALSE.
}
  \item{contm }{the variable names or the column numbers of x that locate the potential continuous mediators. Need to set up only when data is FALSE.}
  \item{catm }{categorical mediators should be binarized and be presented as a list, where the first item is the number of categorical variables and the following items are the names or the column numbers of each binarized categorical variable in x. data.org organizes the categorical mediators in this format after they pass the mediator tests. Need to set up only when data is FALSE.}
  \item{jointm }{a list where the first item is the number of groups of joint mediators to be considered, and each of the following items identifies the names or the column numbers of the mediators in x for each group of joint mediators. Need to set up only when data is FALSE.}
  \item{allm}{
the column numbers of all mediators.  Need to set up only when data is FALSE. The default value of allm is c(contm,catm).
}
  \item{margin}{
the change in predictor when calculating the mediation effects, see Yu et al. (2014).
}
  \item{df}{
Natural cubic spline will be used to fit the relationship between the predictor and each mediator when the predictors are continuous.  The df is the degree of freedom in the ns() function, the default is 1 for linear relationship.
}
  \item{family1}{a list with the ith item 
define the conditional distribution of y[,i] given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is gaussian(link="identity") for contiuous y[,i], and binomial(link = "logit") for binary y[,i].
}
  \item{refy}{
if y is binary, the reference group of y.
}
  \item{binpred}{
if TRUE, the predict variable is binary.  
}
  \item{type}{
the type of prediction when y is class Surv. By default, type is "risk".
}
  \item{w}{
the weight for each case in x.
}
  \item{alpha}{
the elasticnet mixing parameter such that 0<=alpha<=1, with alpha=1 is the lasso penalty, and alpha=0 be the ridge penalty. By default, alpha=1.
}
}
\details{
The mediators are not tested in this function.  data.org.big should be used first for the tests and data organizing, and then the resulting list from data.org.big can be used directly to define the arguments in this function. med considers all variables in x as mediators or covariates in the final model and all variables identified by contm, binm, catm, or jointm as mediators.
}
\value{The result is an med.big object with:
  \item{denm }{a matrix gives the estimated direct effects. Each row is for a different response variable, and each column is for difference predictor.}
  \item{model }{a list, where full.model item is the dth item is the fitted final full model where y[,d] is the outcome and all predictor, covariates, and transformed mediators are the explanatory variables. Other items including family, surv, and biny.}
  \item{coef.m }{a matrix of fitted coefficients for transformed mediators where each row is for a different response variable.} 
  \item{coef.model }{a matrix of fitted coefficients from the full model where each row is for a different response variable.}
  \item{deltaM }{The transformation vector for mediators.}
  \item{dm }{If the predictors are binary, dm is a matrix showing the changes in the mean of mediators when the predictor changes from the reference group to the corresponding group. Each row is for a differenct predictor. If the predictors are continuous, dm is a list. Each item is a matrix for a different predictor. The matrix shows the changes in the mean of mediators when the corresponding predictor changes by one margin. }
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
\code{"\link[=mma.big]{mma.big}"} to make inferences on the estimated mediation effects using bootstrap method.
}
\examples{
#See help(mma.big)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Mediation Analysis }

\name{mma.big}
\alias{mma.big}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple Mediation Analysis for mediaiton analysis with big data sets
}
\description{
Test for mediators and do statistical inferences on the identified mediation effects.
}
\usage{
mma.big(data=NULL,x=data$x, y=data$y,pred=data$dirx, mediator=NULL, binm=data$binm,
        contm=data$contm,catm=data$catm,jointm=data$jointm,margin=1,df=1,binref=NULL,
        catref=NULL,predref=NULL,alpha=1,alpha1=0.01,alpha2=0.01,
        family1=data$family1,n2=50,w=rep(1,nrow(x)),
        refy=NULL,pred.new=NULL,binpred=data$binpred,type=NULL,w.new=NULL,
        lambda=exp(seq(log(0.001), log(5), length.out=15)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{the organized dataset from data.org.big. The default is NULL}
  \item{x}{
a data frame contains the predictor, all potential mediators and covariates.
}
  \item{y}{
the vector of outcome variable.
}
  \item{pred}{
the vector/matrix of the predictor(s).
}
  \item{mediator}{
the list of mediators (column numbers in x or by variable names). The mediators to be checked can be identified by "contmed", "binmed" and "catmed", or by this argument, "mediator", where binary and categorical mediators in x are identified by factors, the reference group is the first level of the factor.  
}
  \item{binm}{
a vector of column numbers that locate the potential binary mediators in x.
}
  \item{contm}{
a vector of column numbers that locate the potential continuous mediators in x.
}
  \item{catm}{
a vector of column numbers that locate the potential categorical mediators in x.
}
  \item{jointm}{
a list that identifies the mediators that need to be forced as mediators.
}
  \item{margin}{
if binpred is FALSE, define the change in predictor when calculating the mediation effects, see Yu et al. (2014).
}
  \item{df}{
Natural cubic spline will be used to fit the relationship between the countinuous predictor and each mediator.  The df is the degree of freedom in the ns() function, the default is 1.
}
  \item{binref}{
the defined reference groups of the binary potential mediators in binm.
}
  \item{catref}{
the defined reference groups of the categorical potential mediators in catm.
}
  \item{predref}{
if predictor is binary, identify the reference group of the binary predictor. By default, the reference group will be the first level of the predictor.
}
  \item{alpha}{
the elasticnet mixing parameter such that 0<=alpha<=1, with alpha=1 is the lasso penalty, and alpha=0 be the ridge penalty. By default, alpha=1.
}
  \item{alpha1}{
the significance level at which to test if the potential mediators (identified by contmed, binmed, and catmed) can be used as a covariate or mediator in estimating y when all variables in x are included in the model.  The default value is alpha=0.1
}
  \item{alpha2}{
the significant level at which to test if a potential mediator is related with the predictor. The default value is alpha2=0.1.
}  
  \item{family1}{
define the conditional distribution of y given x, and the linkage function that links the mean of y with the system component if generalized linear model is used as the final full model.  The default value of family1 is binomial(link = "logit") for binary y, gaussian(link="identity") for continuous y.
}
  \item{n2}{
the number of times of bootstrap resampling.  The default value is n2=50.
}
  \item{w}{
the weight for each observation.
}
  \item{refy}{
if y is binary, the reference group of y.By default, the reference group will be the first level of as.factor(y).
}
  \item{pred.new}{
A new set of predictor(s), of the same format as x (after data.org), on which to calculate the mediation effects. For continuous predictor only.
}
  \item{binpred}{
indicator of whether the predictor is binary or not.
}
  \item{type}{
the type of prediction when y is class Surv. Is "risk" if not specified.
}
  \item{w.new}{
the weights for new.x.
}
  \item{lambda}{
the lambda sequence when fit cv.glmnet. By default, it is exp(seq(log(0.001), log(5), length.out=15)).
}
}
\details{
mma first tests if the potential mediators defined by binm, contm, and catm should be treated as mediators or covariates (if none, the variable will be deleted from further analysis). All variables identified by jointm are treated as mediators. All other variables in x that are not tested are treated as covariates. Then mma does the mediation effects estimation and inference on the selected variables.  
}
\value{
Returns an mma.big object.
  \item{dm }{The matrix of changes in mediator when the predictor change by a margin. Each item of the list is for a different predictor. }  
  \item{coef }{The coefficients of mediators fitted by elastic net. Each item of the list is for a different response variable.} 
  \item{deltaM }{a matrix of deltaM: each row is for a bootstrap sample, and each column is for a different mediator identified by the column name.}
  \item{bootresults }{estimated mediation effects from each bootstrap sample. Each row is for a different mediation effects, identified by the rownames, the columnes is by predictors then by the bootstrap sample order. Each item of the list is for a different response variable.}
  \item{results}{the med.big results using all data.}
  \item{data }{a list that contains all the used data: x=x, y=y, dirx=dirx, binm=binm, contm=contm, catm=catm, jointm=jointm, binpred=F.}
}
\references{
Baron, R.M., and Kenny, D.A. (1986) <doi:10.1037/0022-3514.51.6.1173>. "The moderator-mediator variable distinction in social psychological research: conceptual, strategic, and statistical considerations," J. Pers Soc Psychol, 51(6), 1173-1182.

Yu, Q., Fan, Y., and Wu, X. (2014) <doi:10.4172/2155-6180.1000189>. "General Multiple Mediation Analysis With an Application to Explore Racial Disparity in Breast Cancer Survival," Journal of Biometrics & Biostatistics,5(2): 189. 

Yu, Q., Scribner, R.A., Leonardi, C., Zhang, L., Park, C., Chen, L., and Simonsen, N.R. (2017) <doi:10.1016/j.sste.2017.02.001>. "Exploring racial disparity in obesity: a mediation analysis considering geo-coded environmental factors," Spatial and Spatio-temporal Epidemiology, 21, 13-23. 

Yu, Q., and Li, B. (2017) <doi:10.5334/hors.160>. "mma: An r package for multiple mediation analysis," Journal of Open Research Software, 5(1), 11. 

Li, B. and Yu, Q. (2018). "Mediation Analysis with Large Data Sets". 
}
\author{
Qingzhao Yu  \email{qyu@lsuhsc.edu} and Bin Li \email{bli@lsu.edu}
}
\seealso{
\code{"\link[=data.org.big]{data.org.big}"} is for mediator tests, and \code{"\link[=med.big]{med.big}"} for mediation analysis and inferences.
}
\examples{
#A. binary predictor
set.seed(1)
n=100
pred=rbinom(n,1,0.5)
m1=matrix(rnorm(n*10),n,10)
m2<-matrix(rnorm(n*10,mean=pred,sd=1),n,10)
m3.1=m2[,6:10]
m3=m3.1
m2=m2[,1:5]
m3[m3.1<=0.1]=0
m3[0.1<m3.1 & m3.1<=1]=1
m3[m3.1>1]=2
m3<-apply(m3,2,as.factor)
m<-data.frame(m1,m2,m3)
colnames(m)<-c(paste("m0",1:9,sep=""),paste("m",10:20,sep=""))

lu<--0.5363+0.701*pred+0.801*m[,1]+0.518*m[,2]+1.402*m[,11]+0.773*m[,12]+
    ifelse(m[,16]=="2",2.15,0)+ifelse(m[,16]=="1",0.201,0)

#a. continuous y
y<-rnorm(n,lu,1)

data.e1<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),
                      pred=data.frame(pred),testtype=1)
summary(data.e1,only=TRUE)
data.e1.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred))
summary(data.e1.2,only=TRUE)
med.e1<-med.big(data.e1)
med.e1
med.e1.2<-med.big(data.e1.2)
med.e1.2
mma.e1<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                alpha=1,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
mma.e1.2<-mma.big(data=data.e1.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.e1)
summary(mma.e1.2,RE=TRUE)
plot(mma.e1.2,vari="m16")
plot(mma.e1.2,vari="m11")
joint.effect(mma.e1.2,vari=c("m11","m16")) #find the joint effect of m11 and m16
joint.effect(mma.e1.2,vari="m16") #find the joint effect of m16

\donttest{
#b. binary y
y<-rbinom(n,1,exp(lu)/(1+exp(lu)))
data.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),
                      pred=data.frame(pred),testtype=1)
summary(data.e2,only=TRUE)
data.e2.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                        alpha1=0.05,alpha2=0.05)
summary(data.e2.2,only=TRUE)
med.e2<-med.big(data.e2)
med.e2
med.e2.2<-med.big(data.e2.2)
med.e2.2
mma.e2<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                alpha=1,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
mma.e2.2<-mma.big(data=data.e2.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.e2)
summary(mma.e2.2,RE=TRUE)
plot(mma.e2.2,vari="m16")
plot(mma.e2.2,vari="m11")

#c.time-to-event y
lambda=1/500
survt=-log(runif(n))/lambda/exp(lu)
st=round(runif(n,1,500),0)
time=ifelse(st+survt>600,600,st+survt)-st
cen=ifelse(st+survt>600,0,1)
y=Surv(time,cen)
data.e3<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.e3,only=TRUE)
data.e3.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                        alpha1=0.05,alpha2=0.05)
summary(data.e3.2,only=TRUE)
med.e3<-med.big(data.e3)
med.e3
med.e3.2<-med.big(data.e3.2)
med.e3.2
mma.e3<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha=1,alpha1=0.05,
                alpha2=0.05,n2=3)  #use only the test results.
mma.e3.2<-mma.big(data=data.e3.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.e3)
summary(mma.e3.2,RE=TRUE)
plot(mma.e3.2,vari="m16")
plot(mma.e3.2,vari="m11")

#B. continuous predictor
set.seed(1)
n=100
pred=runif(n,-1,1)
m1=matrix(rnorm(n*10),n,10)
m2<-matrix(rnorm(n*5,mean=0.701*pred,sd=1),n,5)
m3.1=matrix(rnorm(n*5,mean=0.518*pred,sd=1),n,5)
m3=m3.1
m3[m3.1<=0]=0
m3[0<m3.1 & m3.1<=1.28]=1
m3[m3.1>1.28]=2
m3<-apply(m3,2,as.factor)
m<-data.frame(m1,m2,m3)
colnames(m)<-c(paste("m0",1:9,sep=""),paste("m",10:20,sep=""))

lu<--0.6852+0.701*pred+0.801*m[,1]+0.518*m[,2]+1.402*m[,11]+0.773*m[,12]+
    ifelse(m[,16]=="2",2.15,0)+ifelse(m[,16]=="1",0.201,0)

#a. continuous y
y<-rnorm(n,lu,1)

data.c.e1<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                        testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.c.e1,only=TRUE)
data.c.e1.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          alpha1=0.05,alpha2=0.05)
summary(data.c.e1.2,only=TRUE)
med.c.e1<-med.big(data.c.e1)
med.c.e1
med.c.e1.2<-med.big(data.c.e1.2)
med.c.e1.2
mma.c.e1<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha1=0.05,
                  alpha2=0.05,n2=3)
mma.c.e1.2<-mma.big(data=data.c.e1.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.c.e1)
summary(mma.c.e1.2,RE=TRUE)
plot(mma.c.e1,vari="m11")
plot(mma.c.e1,vari="m16")

#b. binary y
y<-rbinom(n,1,exp(lu)/(1+exp(lu)))
data.c.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                        testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.c.e2,only=TRUE)
data.c.e2.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          alpha1=0.05,alpha2=0.05)
summary(data.c.e2.2,only=TRUE)
med.c.e2<-med.big(data.c.e2)
med.c.e2
med.c.e2.2<-med.big(data.c.e2.2)
med.c.e2.2
mma.c.e2<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha1=0.05,
                  alpha2=0.05,n2=3)
mma.c.e2.2<-mma.big(data=data.c.e2.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.c.e2)
summary(mma.c.e2.2,RE=TRUE)
plot(mma.c.e2,vari="m11")
plot(mma.c.e2,vari="m16")

#c.time-to-event y
lambda=1/500
survt=-log(runif(n))/lambda/exp(lu)
st=round(runif(n,1,500),0)
time=ifelse(st+survt>600,600,st+survt)-st
cen=ifelse(st+survt>600,0,1)
y=Surv(time,cen)
data.c.e3<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.c.e3,only=TRUE)
data.c.e3.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          alpha1=0.05,alpha2=0.05)
summary(data.c.e3.2,only=TRUE)
med.c.e3<-med.big(data.c.e3)
med.c.e3
med.c.e3.2<-med.big(data.c.e3.2)
med.c.e3.2
mma.c.e3<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha1=0.05,
                  alpha2=0.05,n2=3)
mma.c.e3.2<-mma.big(data=data.c.e3.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.c.e3)
summary(mma.c.e3.2,RE=TRUE,quant=FALSE)
plot(mma.c.e3,vari="m11")
plot(mma.c.e3,vari="m16")


#C. multiple predictors
#A) categorical predictor
set.seed(1)
n=100
pred=rmultinom(100,1,c(0.5, 0.3, 0.2))
pred=pred[1,]*0+pred[2,]*1+pred[3,]*2
m1=matrix(rnorm(n*10),n,10)
m2<-matrix(rnorm(n*10,mean=pred,sd=1),n,10)
m3.1=m2[,6:10]
m2=m2[,1:5]
m3=m3.1
m3[m3.1<=0.1]=0
m3[0.1<m3.1 & m3.1<=1]=1
m3[m3.1>1]=2
m3<-apply(m3,2,as.factor)
m<-data.frame(m1,m2,m3)
colnames(m)<-c(paste("m0",1:9,sep=""),paste("m",10:20,sep=""))
pred<-as.factor(pred)

lu<--0.5363+ifelse(pred=="1",0.3,0)+ifelse(pred=="2",0.7,0)+0.801*m[,1]+0.518*m[,2]+
    1.402*m[,11]+0.773*m[,12]+ifelse(m[,16]=="2",2.15,0)+ifelse(m[,16]=="1",0.201,0)

#a. continuous y
y<-rnorm(n,lu,1)

data.m.e1<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.m.e1,only=TRUE)
data.m.e1.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred))
summary(data.m.e1.2,only=TRUE)
med.m.e1<-med.big(data.m.e1)
med.m.e1
med.m.e1.2<-med.big(data.m.e1.2)
med.m.e1.2
mma.m.e1<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha=1,
                  alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
mma.m.e1.2<-mma.big(data=data.m.e1.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.m.e1)
summary(mma.m.e1.2,RE=TRUE)
summary(mma.m.e1.2,use.ball=TRUE)
plot(mma.m.e1,vari="m11")
plot(mma.m.e1,vari="m16")

#b. binary y
y<-rbinom(n,1,exp(lu)/(1+exp(lu)))
data.m.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.m.e2,only=TRUE)
data.m.e2.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          alpha1=0.05,alpha2=0.05)
summary(data.m.e2.2,only=TRUE)
med.m.e2<-med.big(data.m.e2)
med.m.e2
#med.m.e2.2<-med.big(data.m.e2.2)
#med.m.e2.2
#mma.m.e2<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha=1,
#                  alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
#mma.m.e2.2<-mma.big(data=data.m.e2.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
#summary(mma.m.e2)
#summary(mma.m.e2.2,RE=TRUE)
#plot(mma.m.e2,vari="m11")
#plot(mma.m.e2,vari="m16")

#c.time-to-event y
lambda=1/500
survt=-log(runif(n))/lambda/exp(lu)
st=round(runif(n,1,500),0)
time=ifelse(st+survt>600,600,st+survt)-st
cen=ifelse(st+survt>600,0,1)
y=Surv(time,cen)
data.m.e3<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.m.e3,only=TRUE)
data.m.e3.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          alpha1=0.05,alpha2=0.05)
summary(data.m.e3.2,only=TRUE)
med.m.e3<-med.big(data.m.e3)
med.m.e3
med.m.e3.2<-med.big(data.m.e3.2)
med.m.e3.2
mma.m.e3<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha=1,
                  alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
mma.m.e3.2<-mma.big(data=data.m.e3.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.m.e3)
summary(mma.m.e3.2,RE=TRUE)
plot(mma.m.e3,vari="m11")
plot(mma.m.e3,vari="m16")

#B) continuous predictor
set.seed(1)
n=100
pred=cbind(runif(n,-1,1),rnorm(n))
m1=matrix(rnorm(n*10),n,10)
m2<-matrix(rnorm(n*5,mean=0.3*pred[,1]+0.4*pred[,2],sd=1),n,5)
m3.1=matrix(rnorm(n*5,mean=0.7*pred[,1]+0.8*pred[,2],sd=1),n,5)
m3=m3.1
m3[m3.1<=0]=0
m3[0<m3.1 & m3.1<=1.28]=1
m3[m3.1>1.28]=2
m3<-apply(m3,2,as.factor)
m<-data.frame(m1,m2,m3)
colnames(m)<-c(paste("m0",1:9,sep=""),paste("m",10:20,sep=""))
colnames(pred)=c("x1","x2")

lu<--0.6852+0.3*pred[,1]+0.7*pred[,2]+0.801*m[,1]+0.518*m[,2]+1.402*m[,11]+0.773*m[,12]+
     ifelse(m[,16]=="2",2.15,0)+ifelse(m[,16]=="1",0.201,0)

#a. continuous y
y<-rnorm(n,lu,1)

data.m.c.e1<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.m.c.e1,only=TRUE)
data.m.c.e1.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                            alpha1=0.05,alpha2=0.05)
summary(data.m.c.e1.2,only=TRUE)
med.m.c.e1<-med.big(data.m.c.e1)
med.m.c.e1
med.m.c.e1.2<-med.big(data.m.c.e1.2)
med.m.c.e1.2
mma.m.c.e1<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                    alpha1=0.05,alpha2=0.05,n2=3)
mma.m.c.e1.2<-mma.big(data=data.c.e1.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.m.c.e1)
summary(mma.m.c.e1.2,RE=TRUE)
plot(mma.m.c.e1,vari="m11")
plot(mma.m.c.e1,vari="m16")

#b. binary y
y<-rbinom(n,1,exp(lu)/(1+exp(lu)))
data.m.c.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.m.c.e2,only=TRUE)
data.m.c.e2.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                            alpha1=0.05,alpha2=0.05)
summary(data.m.c.e2.2,only=TRUE)
med.m.c.e2<-med.big(data.m.c.e2)
med.m.c.e2
med.m.c.e2.2<-med.big(data.m.c.e2.2)
med.m.c.e2.2
pred.new=cbind(runif(n,-1,1),rnorm(n))
colnames(pred.new)=c("x1","x2")
print(med.m.c.e2.2,pred.new=pred.new)
mma.m.c.e2<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha1=0.05,
                    alpha2=0.05,n2=3)
mma.m.c.e2.2<-mma.big(data=data.m.c.e2.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.m.c.e2)
summary(mma.m.c.e2.2,RE=TRUE)
plot(mma.m.c.e2,vari="m11")
plot(mma.m.c.e2,vari="m16")


#c.time-to-event y
lambda=1/500
survt=-log(runif(n))/lambda/exp(lu)
st=round(runif(n,1,500),0)
time=ifelse(st+survt>600,600,st+survt)-st
cen=ifelse(st+survt>600,0,1)
y=Surv(time,cen)
data.m.c.e3<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.m.c.e3,only=TRUE)
data.m.c.e3.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                            alpha1=0.05,alpha2=0.05)
summary(data.m.c.e3.2,only=TRUE)
med.m.c.e3<-med.big(data.m.c.e3)
med.m.c.e3
med.m.c.e3.2<-med.big(data.m.c.e3.2)
med.m.c.e3.2
mma.m.c.e3<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha1=0.05,
                    alpha2=0.05,n2=3)
mma.m.c.e3.2<-mma.big(data=data.m.c.e3.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.m.c.e3)
summary(mma.m.c.e3.2,RE=TRUE)
plot(mma.m.c.e3,vari="m11")
plot(mma.m.c.e3,vari="m16")


#d. multivariate y
y<-cbind(rnorm(n,lu,1),rbinom(n,1,exp(lu)/(1+exp(lu))))
colnames(y)=c("y1","y2")
data.m.m.c.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                            testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.m.m.c.e2,only=TRUE)
data.m.m.c.e2.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                              alpha1=0.05,alpha2=0.05)
summary(data.m.m.c.e2.2,only=TRUE)
med.m.m.c.e2<-med.big(data.m.m.c.e2)
med.m.m.c.e2
med.m.m.c.e2.2<-med.big(data.m.m.c.e2.2)
med.m.m.c.e2.2
pred.new=cbind(runif(n,-1,1),rnorm(n))
colnames(pred.new)=c("x1","x2")
print(med.m.m.c.e2.2,pred.new=pred.new)
mma.m.m.c.e2<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha1=0.05,
                      alpha2=0.05,n2=3)
mma.m.m.c.e2.2<-mma.big(data=data.m.m.c.e2.2,alpha1=0.05,alpha2=0.05,n2=3)  
summary(mma.m.m.c.e2)
summary(mma.m.m.c.e2.2,RE=TRUE)
plot(mma.m.m.c.e2,vari="m11")
plot(mma.m.m.c.e2,vari="m16")
mma.m.m.c.e2
mma.m.m.c.e2.2}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Mediation Analysis }
\keyword{ Mediator Tests }% __ONLY ONE__ keyword per line

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error=TRUE,
  warning=FALSE
)

## ---- include=F---------------------------------------------------------------
library(mmabig)

## -----------------------------------------------------------------------------
# a binary predictor
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

# a continuous y
y<-rnorm(n,lu,1)

## -----------------------------------------------------------------------------
data.e1<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),
                      pred=data.frame(pred),testtype=1)
summary(data.e1,only=TRUE)

## -----------------------------------------------------------------------------
data.e1.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred))
summary(data.e1.2,only=TRUE)

## -----------------------------------------------------------------------------
lambda=1/500
survt=-log(runif(n))/lambda/exp(lu)
st=round(runif(n,1,500),0)
time=ifelse(st+survt>600,600,st+survt)-st
cen=ifelse(st+survt>600,0,1)
y=Surv(time,cen)

data.e3<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.e3,only=TRUE)

## ---- include=F---------------------------------------------------------------
data.e3.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                        alpha1=0.05,alpha2=0.05)
summary(data.e3.2,only=TRUE)

## -----------------------------------------------------------------------------
# multicategorical predictor
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
# continuous y
lu<--0.5363+ifelse(pred=="1",0.3,0)+ifelse(pred=="2",0.7,0)+0.801*m[,1]+0.518*m[,2]+
    1.402*m[,11]+0.773*m[,12]+ifelse(m[,16]=="2",2.15,0)+ifelse(m[,16]=="1",0.201,0)
y<-rnorm(n,lu,1)

data.m.e1<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),testtype=1)
summary(data.m.e1,only=TRUE)

## -----------------------------------------------------------------------------
# multivariate predictor
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
# binary y
lu<--0.6852+0.3*pred[,1]+0.7*pred[,2]+0.801*m[,1]+0.518*m[,2]+1.402*m[,11]+0.773*m[,12]+
     ifelse(m[,16]=="2",2.15,0)+ifelse(m[,16]=="1",0.201,0)
y<-rbinom(n,1,exp(lu)/(1+exp(lu)))

data.m.c.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                          testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.m.c.e2,only=TRUE)

## -----------------------------------------------------------------------------
# multivariate predictor
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
#multivariate responses
y<-cbind(rnorm(n,lu,1),rbinom(n,1,exp(lu)/(1+exp(lu))))
colnames(y)=c("y1","y2")

data.m.m.c.e2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                            testtype=1,alpha1=0.05,alpha2=0.05)
summary(data.m.m.c.e2,only=TRUE)

data.m.m.c.e2.2<-data.org.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                              alpha1=0.05,alpha2=0.05)
summary(data.m.m.c.e2.2,only=TRUE)

## -----------------------------------------------------------------------------
med.e1<-med.big(data.e1)

## -----------------------------------------------------------------------------
med.e1

## ---- include = FALSE---------------------------------------------------------
med.e3<-med.big(data.e3)
med.e3

## -----------------------------------------------------------------------------
med.m.m.c.e2.2<-med.big(data.m.m.c.e2.2)
med.m.m.c.e2.2

## ----include=F----------------------------------------------------------------
#print.med.big can be used to print the estimation of mediation effects from an med.big object with a user-defined new set of predictors.
pred.new=cbind(runif(10,-1,1),rnorm(10))
colnames(pred.new)=c("x1","x2")
print(med.m.m.c.e2.2,pred.new=pred.new)

## ---- fig.show='hold', fig.height=5, fig.width=7------------------------------
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
y<-rnorm(n,lu,1)

mma.e1<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),
                alpha=1,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.e1)

## ---- fig.show='hold', fig.height=5, fig.width=7------------------------------
mma.e1.2<-mma.big(data=data.e1.2,alpha1=0.05,alpha2=0.05,n2=3) 
summary(mma.e1.2,RE=TRUE)

## ---- fig.height=7, fig.width=5-----------------------------------------------
#plot(mma.e1.2,vari="m16")
plot(mma.e1.2,vari="m11")

## ---- fig.show='hold', fig.height=5, fig.width=7------------------------------
lambda=1/500
survt=-log(runif(n))/lambda/exp(lu)
st=round(runif(n,1,500),0)
time=ifelse(st+survt>600,600,st+survt)-st
cen=ifelse(st+survt>600,0,1)
y=Surv(time,cen)

mma.e3<-mma.big(x=m,y=data.frame(y),mediator=1:ncol(m),pred=data.frame(pred),alpha=1,alpha1=0.05,
                alpha2=0.05,n2=3)  #use only the test results.
mma.e3.2<-mma.big(data=data.e3.2,alpha1=0.05,alpha2=0.05,n2=3)  #use only the test results.
summary(mma.e3)
summary(mma.e3.2,RE=TRUE,quant=FALSE)

## ---- fig.height=7, fig.width=5-----------------------------------------------
plot(mma.e3.2,vari="m16")
plot(mma.e3.2,vari="m11")


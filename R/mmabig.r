#to organize data
data.org.big<-function(x,y,pred,mediator=NULL,contmed=NULL,binmed=NULL,binref=NULL,catmed=NULL,
                       catref=NULL,jointm=NULL, 
                       family1=as.list(rep(NA,ncol(data.frame(y)))),
                       predref=NULL,alpha=1,alpha1=0.01,alpha2=0.01,testtype=2, w=NULL,lambda=exp(seq(log(0.001), log(5), length.out=15)))
{cattobin<-function(x,cat1,cat2=rep(1,length(cat1)))
{ ad1<-function(vec)
{vec1<-vec[-1]
 vec1[vec[1]]<-1
 vec1
}
 dim1<-dim(x)
 catm<-list(n=length(cat1))
 g<-dim1[2]-length(cat1)
 ntemp<-names(x)[cat1]
 j<-1
 for (i in cat1)
 {a<-factor(x[,i])
  d<-rep(0,dim1[1])
  b<-sort(unique(a[a!=cat2[j]]))
  l<-1
  for (k in b)
  {d[a==k]<-l
   l<-l+1}
  d[a==cat2[j]]<-l
  f<-matrix(0,dim1[1],l-1) 
  colnames(f)<-paste(ntemp[j],b,sep="") #changed for error info
  hi<-d[d!=l & !is.na(d)]
  f[d!=l & !is.na(d),]<-t(apply(cbind(hi,f[d!=l & !is.na(d),]),1,ad1))
  f[is.na(d),]<-NA
  x<-cbind(x,f)
  catm<-append(catm,list((g+1):(g+l-1)))
  g<-g+length(b)
  j<-j+1
 }
 x<-x[,-cat1]
 list(x=x,catm=catm)
}

anymissing<-function(vec) #return T if there is any missing in the vec
{if(sum(is.na(vec))>0)
  return(F)
  else return(T)
}

colnum<-function(vec,cutx) 
{z<-vec
 for (i in 1:length(vec))
   z[i]<-vec[i]-sum(vec[i]>cutx)
 z}

#y2 is the response dataframe. Family changes to the family for glmnet but not glm
#did not use the multivariate responses function
y2<-data.frame(y)                     
ny<-ncol(y2)
y_type<-rep(4,ny)                     #1 is continuous, 2 is binary, 3 is multi-categorical, 4 is survival
if(is.null(family1))
  family1=as.list(rep(NA,ncol(data.frame(y))))
for (i in 1:ny)
{if(class(y2[,i])!="Surv"){
  if(nlevels(droplevels(as.factor(y2[,i])))==2)   
  {y_type[i]<-2
   if(is.na(family1[[i]]))
    family1[[i]] = "binomial"        #family1 changed according to glmnet
  }
 else if(is.character(y2[,i]) | is.factor(y2[,i]))
 {y_type[i]<-3
  y2[,i]<-droplevels(y2[,i])  #drop those empty levles
  if(is.na(family1[[i]]))
    family1[[i]] = "multinomial"        
 }
 else
 {y_type[i]=1
  if(is.na(family1[[i]]))
   family1[[i]] = "gaussian"
 }}
else
{if(is.na(family1[[i]]))
  family1[[i]] = "cox"
}
}

if(sum(y_type==3)>0) #transfer the multicategorical type response
#the last level is used as the reference (note:by cat2bin, 0 for all columns) --need to check with glmnet which use the 1 as the target group(????)
{temp1<-(1:length(y_type))[y_type==3]
 temp2<-cattobin(y2,temp1,tail(nlevels(as.factor(y2[,temp1])),n=1))
 y2<-data.frame(temp2$x)
 y_type<-y_type[-temp1]
 y_type<-c(y_type,rep(2,ncol(y2)-length(y_type)))
 family1[[temp1]]<-NULL
 family1<-append(family1,rep(list("binomial"),ncol(y2)-length(family1)))
}

xnames<-colnames(x)

#predictors have to be one categorical or all continuous, pred is the exposure vector/matrix
pred1<-data.frame(pred)
binpred=T
for (i in 1:ncol(pred1))
 if(nlevels(as.factor(pred1[,i]))==2)
 {if(!is.null(predref))
      pred1[,i]<-as.factor(ifelse(pred1[,i]==predref,0,1))
  else
      {pred1[,i]<-as.factor(pred1[,i])
       pred1[,i]<-as.factor(ifelse(pred1[,i]==levels(pred1[,i])[1],0,1))}
 }
 else if(is.character(pred1[,i]) | is.factor(pred1[,i]))
 {if(!is.null(predref))
   pred1<-cattobin(data.frame(pred1[,i]),1,predref)$x
  else
   pred1<-cattobin(data.frame(pred1[,i]),1,levels(as.factor(pred1[,i]))[1])$x
 }
else
  binpred=F
pred<-data.frame(pred1)

if(!is.null(jointm))
 jointm<-list(n=1,jointm=jointm)

if(is.character(contmed))
  contmed<-unlist(sapply(contmed,grep,xnames))
if(is.character(binmed))
  binmed<-unlist(sapply(binmed,grep,xnames))
if(is.character(catmed))
  catmed<-unlist(sapply(catmed,grep,xnames))
if(!is.null(jointm))
  for (i in 2:length(jointm))
   if(is.character(jointm[[i]]))
     jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))

if(!is.null(binmed))  #revise: binref does not have to be null or full, & is.null(binref)
 {j<-1
  for(i in binmed)
  {x[,i]=as.factor(x[,i])
   if(is.null(binref))
   binref[i]=levels(x[,i])[1]
   else if(is.na(binref[i]))
     binref[i]=levels(x[,i])[1]
   j<-j+1}}

if(!is.null(catmed))   
 {j<-1  
  for(i in catmed)
  {x[,i]=as.factor(x[,i])
    if(is.null(catref))
      catref[i]=levels(x[,i])[1]
    else if(is.na(catref[j]))
      catref[i]=levels(x[,i])[1]
   j<-j+1}}

if(!is.null(mediator))   #for all potential mediators
 {if(is.character(mediator))
   mediator<-unlist(sapply(mediator,grep,xnames))
  for (i in 1:length(mediator))
    {if(is.character(x[,mediator[i]]))
      x[,mediator[i]]<-as.factor(x[,mediator[i]])
     if(is.factor(x[,mediator[i]]) | nlevels(as.factor(x[,mediator[i]]))==2)
       if(nlevels(as.factor(x[,mediator[i]]))==2)
       {if(sum(binmed==mediator[i])==0)
        {x[,mediator[i]]<-as.factor(x[,mediator[i]])
         binmed<-c(binmed,mediator[i])
         binref<-c(binref,(levels(x[,mediator[i]]))[1])
         }}
      else
       {if(sum(catmed==mediator[i])==0)
        {catmed<-c(catmed,mediator[i])
         catref<-c(catref,levels(x[,mediator[i]])[1])}}
     else if(sum(contmed==mediator[i])==0 & sum(catmed==mediator[i])==0)
        contmed<-c(contmed,mediator[i])}
 }

if(!is.null(jointm))  #mediators that should be jointly considered are forced in as mediators
{joint<-NULL 
 for (i in 2:(jointm[[1]]+1))
   joint=c(joint,jointm[[i]])
 joint1<-unique(joint)
 
 if(!is.null(contmed))
 {cont1<-rep(F,length(contmed))
  for (i in 1:length(contmed))
    cont1[i]<-ifelse(sum(contmed[i]==joint1)>0,T,F)
 }
 if(!is.null(binmed))
 {bin1<-rep(F,length(binmed))
  for (i in 1:length(binmed))
    bin1[i]<-ifelse(sum(binmed[i]==joint1)>0,T,F)
 }
 if(!is.null(catmed))
 {cat1<-rep(F,length(catmed))
  for (i in 1:length(catmed))
    cat1[i]<-ifelse(sum(catmed[i]==joint1)>0,T,F)
 }
}
else
{if(!is.null(contmed))
  cont1<-rep(F,length(contmed))
 
 if(!is.null(binmed))
   bin1<-rep(F,length(binmed))
 
 if(!is.null(catmed))
   cat1<-rep(F,length(catmed))
}

if(!is.null(binmed))
{j<-1
 for (i in binmed)
 {x[,i]<-ifelse(x[,i]==binref[j],0,1)
  j<-j+1}
}

if(!is.null(catmed))
{tempx<-cattobin(x,catmed,catref)
 newx<-tempx$x
 catnewx<-tempx$catm}
else newx<-x

fullmodel1<-NULL
x.all<-cbind(x,pred)
for (j in 1:length(y_type))
{y.all<-y2[,j]
 nonmissing<-apply(cbind(y.all,x.all),1,anymissing)
 y.all<-y.all[nonmissing]
 x.all1<-x.all[nonmissing,]
 if(is.null(w) & y_type[j]!=4)
   w1<-rep(1,length(y.all))
 else if(is.null(w))
   w1<-rep(1,nrow(y.all))
 else
   w1<-w[nonmissing]
 x.temp<-model.matrix(y.all~.,x.all1)
 fit.all<-suppressWarnings(cv.glmnet(x.temp,y.all,weights=w1,family=family1[[j]],alpha=alpha,lambda=lambda))
 if(y_type[j]!=4)
   fullmodel1<-cbind(fullmodel1,coef(fit.all)[-(1:2),])
 else     fullmodel1<-cbind(fullmodel1,coef(fit.all)[-1,]) # survival outcome does not fit an intercept
}
xname<-names(x.all)
xnames3<-rownames(fullmodel1)
if(testtype==2)
  {P1<-matrix(NA,length(xname),ncol(y2))  #the type III for predictor and one mediator only model
   rownames(P1)<-xname}
else
  P1<-fullmodel1
  
colnames(P1)<-colnames(y2)
prednames<-colnames(pred)

covr.cont<-rep(F,length(contmed))
covr.bin<-rep(F,length(binmed))
covr.cat<-rep(F,length(catmed))

for (j in 1:ncol(y2))  
{if(testtype==2)        #univariate comparison
 {if(y_type[j]==4 & is.null(w))
   temp.fullmodel1<-coxph(y2[,j]~.,data=pred)
  else if (y_type[j]==4)
   temp.fullmodel1<-coxph(y2[,j]~.,data=pred,weights=w)
  else
   temp.fullmodel1<-glm(y2[,j]~.,data=pred,family=family1[[j]],weights=w) 
  temp.type3<-summary(temp.fullmodel1)$coefficients#Anova(temp.fullmodel1,type="III")
  for (k in 1:ncol(pred))
    if(y_type[j]==4)
      P1[grep(prednames[k],xname),j]<-temp.type3[grep(prednames[1],rownames(temp.type3)),5]
    else
      P1[grep(prednames[k],xname),j]<-temp.type3[grep(prednames[1],rownames(temp.type3)),4]
}

 if(!is.null(contmed))
  for (i in 1:length(contmed))
   if(testtype==1)
    covr.cont[i]<-ifelse(sum(fullmodel1[grep(xname[contmed[i]],xnames3),j])!=0,T,covr.cont[i])
   else if(testtype==2)
   {temp.data<-cbind(x[,contmed[i]],pred)
    names(temp.data)<-c(xname[contmed[i]],names(pred))
    if(y_type[j]==4)
     temp.fullmodel1<-coxph(y2[,j]~.,data=data.frame(temp.data),weights=w)
    else
     temp.fullmodel1<-glm(y2[,j]~.,data=data.frame(temp.data),family=family1[[j]],weights=w) 
    temp.type3<-summary(temp.fullmodel1)$coefficients
    if(y_type[j]==4)
      temp.p1<-temp.type3[rownames(temp.type3)==xname[contmed[i]],5]
    else
      temp.p1<-temp.type3[rownames(temp.type3)==xname[contmed[i]],4]
    P1[grep(xname[contmed[i]],xname),j]<-temp.p1
    covr.cont[i]<-ifelse(temp.p1<alpha1,T,covr.cont[i])
   }

 if(!is.null(binmed))
  for (i in 1:length(binmed))
   if(testtype==1)
    covr.bin[i]<-ifelse(fullmodel1[xnames3==xname[binmed[i]],j]!=0,T,covr.bin[i])
   else if(testtype==2)
   {temp.data<-cbind(x[,binmed[i]],pred)
    names(temp.data)<-c(xname[binmed[i]],names(pred))
    if(y_type[j]==4)
     temp.fullmodel1<-coxph(y2[,j]~.,data=data.frame(temp.data),weights=w)
    else
     temp.fullmodel1<-glm(y2[,j]~.,data=data.frame(temp.data),family=family1[[j]],weights=w) 
   temp.type3<-summary(temp.fullmodel1)$coefficients
   if(y_type[j]==4)
     temp.p1<-temp.type3[rownames(temp.type3)==xname[binmed[i]],5]
   else
     temp.p1<-temp.type3[rownames(temp.type3)==xname[binmed[i]],4]
   P1[grep(xname[binmed[i]],xname),j]<-temp.p1
   covr.bin[i]<-ifelse(temp.p1<alpha1,T,covr.bin[i])
  }
 
 if(!is.null(catmed))
  for (i in 1:length(catmed))
   if(testtype==1)
     covr.cat[i]<-ifelse(sum(fullmodel1[grep(xname[catmed[i]],xnames3),j])!=0,T,covr.cat[i]) 
   else if(testtype==2)
    {temp.data<-cbind(x[,catmed[i]],pred)
     names(temp.data)<-c(xname[catmed[i]],names(pred))
     if(y_type[j]==4)
      temp.fullmodel1<-coxph(y2[,j]~.,data=data.frame(temp.data),weights=w)
     else
      temp.fullmodel1<-glm(y2[,j]~.,data=data.frame(temp.data),family=family1[[j]],weights=w) #use type three error to test the full model
    temp.type3<-summary(temp.fullmodel1)$coefficients
    if(y_type[j]==4)
      temp.p1<-temp.type3[grep(xname[catmed[i]],rownames(temp.type3)),5]
    else
      temp.p1<-temp.type3[grep(xname[catmed[i]],rownames(temp.type3)),4]
    P1[grep(xname[catmed[i]],xname),j]<-min(temp.p1)
    covr.cat[i]<-ifelse(min(temp.p1)<alpha1,T,covr.cat[i])
   }
} 

if(!is.null(contmed))
 {covr.cont<-ifelse(covr.cont|cont1,T,F)
  cont2<-cont1[covr.cont]
  contmed1<-contmed[covr.cont]}
if(!is.null(binmed))
 {covr.bin<-ifelse(covr.bin+bin1>0,T,F) 
  bin2<-bin1[covr.bin]
  binmed1<-binmed[covr.bin]}
if(!is.null(catmed))
 {covr.cat<-ifelse(covr.cat+cat1>0,T,F)
  cat2<-cat1[covr.cat]
  catmed1<-catmed[covr.cat]
  catref1<-catref[covr.cat]}

a1<-c(contmed,binmed,catmed)
a2<-c(covr.cont,covr.bin,covr.cat)

cutx<-a1[!a2]

if (sum(a2)==0)
  return ("no mediators found")
else if(length(cutx)==0)
{newx1<-x
 contm1<-contmed
 binm1<-binmed
 catm1<-catmed
 catref1<-catref
}
else {newx1<-x[,-cutx]
      if(sum(covr.cont)==0)
       contm1<-NULL  
      else 
       contm1<-colnum(contmed[covr.cont],cutx)
      if(sum(covr.bin)==0)
        binm1<-NULL
      else 
        binm1<-colnum(binmed[covr.bin],cutx)
      if(sum(covr.cat)==0)
      {catm1<-NULL
       catref1<-NULL}
      else 
      {catm1<-colnum(catmed[covr.cat],cutx)
       catref1<-catref[covr.cat]}
      if(!is.null(jointm))
        for(i in 2:(jointm[[1]]+1))
          jointm[[i]]<-colnum(jointm[[i]],cutx)
}

#delete nonmediators
rela_var<-NULL               #to store the relationships between mediators and predictor   
rela_p<-NULL
name_newx<-names(newx1)
if (binpred)             #for binary (x)
{contm2<-contm1
 if(length(contm1)>0)
  {med.cont<-rep(F,length(contm1))
   for (i in 1:length(contm1))
   {tempmodel<-summary(glm(newx1[,contm1[i]]~.,weights=w,data=pred)) #allowing multivariate predictors
    med.cont[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
    rela_var<-c(rela_var,name_newx[contm1[i]])
    rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
   }
  med.cont<-ifelse(med.cont+cont2>0,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
    {tempmodel<-summary(glm(newx1[,binm1[i]]~.,weights=w,family="binomial",data=pred)) #allowing multivariate predictors
     med.bin[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
     rela_var<-c(rela_var,name_newx[binm1[i]])
     rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
    }
  med.bin<-ifelse(med.bin+bin2>0,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))  
   {temp.p<-NULL                                 #allowing multivariate predictors
    for (j in 1:ncol(pred))          
     temp.p<-c(temp.p,min(summary(glm(pred[,j]~newx1[,catm1[i]],weights=w,family="binomial"))$coef[-1,4]))
    med.cat[i]<-ifelse(min(temp.p)<alpha2,T,F)
    rela_var<-c(rela_var,name_newx[catm1[i]])
    rela_p<-rbind(rela_p,temp.p)
   }
  med.cat<-ifelse(med.cat+cat2>0,T,F)
  cat3<-cat2[med.cat]
  catm2<-catm1[med.cat]
  catref2<-catref1[med.cat]}
}
else
{contm2<-contm1
 if(length(contm1)>0)
 {med.cont<-rep(F,length(contm1))
  for (i in 1:length(contm1))
  {tempmodel<-summary(glm(newx1[,contm1[i]]~.,weights=w,data=pred))
   med.cont[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[contm1[i]])
   rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
  }
 med.cont<-ifelse(med.cont|cont2,T,F)
  contm2<-contm1[med.cont]}
 
 binm2<-binm1
 if(length(binm1)>0) 
 {med.bin<-rep(F,length(binm1))
  for (i in 1:length(binm1))   
  {tempmodel<-summary(glm(newx1[,binm1[i]]~.,weights=w,family="binomial",data=pred))
   med.bin[i]<-ifelse(min(tempmodel$coef[-1,4])<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[binm1[i]])
   rela_p<-rbind(rela_p,tempmodel$coef[-1,4])
  }    
  med.bin<-ifelse(med.bin|bin2,T,F)
  binm2<-binm1[med.bin]}
 
 catm2<-catm1
 if(length(catm1)>0) 
 {med.cat<-rep(F,length(catm1))
  for (i in 1:length(catm1))   
  {temp.p<-NULL
   for(j in 1:ncol(pred))
    temp.p<-c(temp.p,min(summary(glm(pred[,j]~newx1[,catm1[i]],weights=w))$coef[-1,4]))
   med.cat[i]<-ifelse(min(temp.p)<alpha2,T,F)
   rela_var<-c(rela_var,name_newx[catm1[i]])
   rela_p<-rbind(rela_p,temp.p)
  }
  med.cat<-(med.cat|cat2)
  cat3<-cat2[med.cat]
  catm2<-catm1[med.cat]
  catref2<-catref1[med.cat]}
}

if(length(catm2)==0)
  catm2<-NULL
if(length(binm2)==0)
  binm2<-NULL
if(length(contm2)==0)
  contm2<-NULL

newx2<-newx1

catm<-NULL
 if(!is.null(catm2))
 {tempx<-cattobin(newx1,catm2,catref2)
  newx2<-tempx$x
  catm<-tempx$catm
  if(!is.null(binm2))
    binm2<-colnum(binm2,catm2)
  if(!is.null(contm2))
    contm2<-colnum(contm2,catm2)
  if (!is.null(jointm))
  {if (sum(cat3)==0)
    for(i in 2:(jointm[[1]]+1))
      jointm[[i]]<-colnum(jointm[[i]],catm2)
   else
     for(i in 2:(jointm[[1]]+1))
     {a<-jointm[[i]]
      b<-NULL
      for (j in a)
        if(sum(j==catm2)==0)
          b<-c(b,colnum(j,catm2))
      else for (k in 1:length(catm2))
        if (j==catm2[k])
          b<-c(b,catm[[k+1]])
      jointm[[i]]<-b}
  }
 }
 rownames(rela_p)<-rela_var
 results<-list(x=newx2,dirx=pred,contm=contm2,binm=binm2,catm=catm, jointm=jointm, y=y2, y_type=y_type,
               fullmodel=fullmodel1,rela=rela_p,binpred=binpred,family1=family1,
               testtype=testtype,P1=P1,w=w,testtype=testtype)
 class(results)<-"med_iden.big"
 return(results)
#}
}

summary.med_iden.big<-function(object,...,only=F)
{var.name<-colnames(object$x)
if(is.list(object$catm))                  #revised to show catm when it is a list (when x is continuous)
{t.catm<-NULL
for (i in 2:length(object$catm))
  t.catm<-c(t.catm,object$catm[[i]])}
else
{t.catm<-object$catm}
mediator<-var.name[c(object$contm,object$binm,t.catm)]
covariates<-var.name[-c(object$contm,object$binm,t.catm)]
tests<-NULL
tests<-object$P1 
temp<-rownames(tests)
rownames(tests)<-temp
temp.name<-rownames(object$rela)
temp2<-matrix(NA,length(temp),ncol(object$rela))
for (i in 1:nrow(object$rela))
  temp2[grep(temp.name[i],temp),]<-object$rela[i,]
tests<-cbind(tests,temp2)
dimnames(tests)[[2]]<-c(paste("P-Value 1",colnames(object$y),sep="."), paste("P-Value 2", colnames(object$dirx),sep="."))
result<-list(mediator=mediator, covariates=covariates,tests=tests, results=object,only=only,testtype=object$testtype)
class(result)<-"summary.med_iden.big"
result
}

print.summary.med_iden.big<-function(x,...)  #version 6: changed typos in the function--tests->x$tests
{cat("Identified as mediators: \n")
  print(x$mediator)
  cat("Selected as covariates: \n")
  print(x$covariates)
  cat("Tests: \n")
  temp<-rownames(x$tests)
  tests.1<-NULL
  for (z in 1:length(temp))
    if(length(grep(temp[z],x$mediator))>0)
      dimnames(x$tests)[[1]][z]<-paste(temp[z],"*")
  if (!is.null(x$results$jointm))
  {tt<-NULL
  for (i in 1:(length(x$results$jointm)-1)) 
    tt<-c(tt,x$results$jointm[[i+1]])
  tt<-unique(tt)
  for (z in 1:length(temp))
    if(length(grep(temp[z],names(x$results$x)[tt]))>0)
      dimnames(x$tests)[[1]][z]<-paste(temp[z],"-")}
  if(x$testtype==2)
    dimnames(x$tests)[[2]]<-c(paste("P-Value 1",colnames(x$results$y),sep="."), paste("P-Value 2", colnames(x$results$dirx),sep="."))
  else
    dimnames(x$tests)[[2]]<-c(paste("Coefficients",colnames(x$results$y),sep="."), paste("P-Value 2", colnames(x$results$dirx),sep="."))
  if(!x$only)
    print(round(x$tests,3))
  else
  {temp.name<-NULL
  temp1<-rownames(x$tests)
  for (z in 1:length(temp))
    if(length(grep(temp[z],x$mediator))>0)
    {tests.1<-rbind(tests.1,x$tests[z,])
    temp.name<-c(temp.name,temp1[z])}
  else if(length(grep(temp[z],x$covariates))>0) 
  {tests.1<-rbind(tests.1,x$tests[z,])
  temp.name<-c(temp.name,temp1[z])}
  rownames(tests.1)<-temp.name
  print(round(tests.1,3))
  } 
  if(x$testtype==2)
    cat("----\n *:mediator,-:joint mediator\n P-Value 1:univariate relationship test with the outcome;P-Value 2:Tests of relationship with the Predictor\n")
  else
    cat("----\n *:mediator,-:joint mediator\n Coefficients: estimated coefficients; P-Value 2:Tests of relationship with the Predictor\n")
  
}


med.big<-function(data, x=data$x, y=data$y, dirx=data$dirx, binm=data$binm, contm = data$contm, 
                        catm = data$catm, jointm = data$jointm, allm = c(contm, binm, unique(unlist(catm)[-1])), 
                        margin=1,df=1,family1=data$family1,refy=rep(NA,ncol(y)),binpred=data$binpred,type=NULL, 
                        w=NULL,alpha=1,lambda=exp(seq(log(0.001), log(5), length.out=15)))
  {if (is.null(c(binm,contm,catm)))
    stop("Error: no potential mediator is specified")
   
  if(is.null(data)){
    surv=rep(F,ncol(y))
    biny=rep(F,ncol(y))
    for(j in 1:ncol(y)) {
      if(class(y[,j])=="Surv"){
        surv[j]=T
        if(is.null(type))
          type="response"
        else if (is.null(type))
          type="risk"
      }
      else if(is.character(y[,j]) | is.factor(y[,j]) | nlevels(as.factor(y[,j]))==2)
      {biny[j]=T
      if(is.na(family1[[j]]))
        family1[[j]] = binomial("logit")
      if(!is.na(refy[j]))
        y[,j]<-ifelse(y[,j]==refy[j],0,1)
      else
        y[,j]<-ifelse(as.factor(y[,j])==levels(as.factor(y[,j]))[1],0,1)
      }
      else { 
        if(is.na(family1[[j]]))
          family1[[j]] = gaussian(link = "identity")
      }
    }
  }
  else
  {biny=data$y_type==2
  surv=data$y_type==4
  if(sum(surv)>0 & is.null(type))
    type="response"
  for (j in 1:ncol(y))
   if(biny[j] & !is.na(refy[j]))
     y[,j]<-ifelse(y[,j]==refy[j],0,1)
   else if(biny[j])
     y[,j]<-ifelse(as.factor(y[,j])==levels(as.factor(y[,j]))[1],0,1)  #the first level is the reference group
  }
  
    xnames<-colnames(x)
    pred_names<-colnames(dirx)
    ynames<-colnames(y)
    if(is.character(contm))
      contm<-unlist(sapply(contm,grep,xnames))
    if(is.character(binm))
      binm<-unlist(sapply(binm,grep,xnames))
    if(!is.null(catm))
      for (i in 2:length(catm))
        if(is.character(catm[[i]]))
          catm[[i]]<-unlist(sapply(catm[[i]],grep,xnames))
    if(!is.null(jointm))
      for (i in 2:length(jointm))
        if(is.character(jointm[[i]]))
          jointm[[i]]<-unlist(sapply(jointm[[i]],grep,xnames))
    
    anymissing<-function(vec) #return T if there is any missing in the vec
    {if(sum(is.na(vec))>0)
      return(F)
     else return(T)
    }
    
    exp.m.given.x<-function(x,dirx,binm,contm=NULL,df,w) #give the model and residual of m given x
    {models<-NULL
     res<-NULL
     prednames=names(dirx)
     if (ncol(dirx)==1)
       expr=paste("ns(",prednames,",df=",df,")",sep="")
     else
       {expr=paste("ns(",prednames[1],",df=",df, ")",sep="")
        for (i in 2:ncol(dirx))
          expr=paste(expr,paste("ns(",prednames[i],",df=",df,")",sep=""),sep="+")
       }
    j<-1
    forml=paste("y~",expr,sep="")
    if(!is.null(binm))
    {for(i in binm)
     {y=x[,i]
      if(is.null(w))
        models[[j]]<-glm(forml,data=data.frame(y,dirx),family=binomial(link = "logit"))
      else
        models[[j]]<-glm(forml,data=data.frame(y,dirx,w),family=binomial(link = "logit"),weights=w)
      #res<-cbind(res,x[,i]-predict(models[[j]],type = "response",newdata=data.frame(z=z)))
      j<-j+1}
    }
    for (i in contm)
    {y=x[,i]
    if(is.null(w))
      models[[j]]<-glm(forml,data=data.frame(y,dirx),family=gaussian(link="identity"))
    else
      models[[j]]<-glm(forml,data=data.frame(y,dirx,w),family=gaussian(link="identity"),weights=w)
    #res<-cbind(res,models[[j]]$res)
     j<-j+1
    }
    list(models=models)
    }
    
    
    dM<-function(distmgivenx,dirx,binm1,contm,df,margin,m.names)  
    {means<-NULL
    if(!is.null(binm1))
      for (i in 1:length(binm1))
        means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(dirx)))
    if(!is.null(contm))
      for (i in (length(binm1)+1):length(c(binm1,contm)))
        means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(dirx)))
    
    dm<-vector("list", ncol(dirx))
    names(dm)<-colnames(dirx)
    for (j in 1:ncol(dirx)){  
     dirx1<-dirx
     dirx1[,j]=dirx[,j]+margin
     
     means1<-NULL
     if(!is.null(binm1))
       for (i in 1:length(binm1))
         means1<-cbind(means1,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(dirx1)))
     if(!is.null(contm))
       for (i in (length(binm1)+1):length(c(binm1,contm)))
         means1<-cbind(means1,predict(distmgivenx$models[[i]],newdata=data.frame(dirx1)))
     dm[[j]]<-means1-means
     colnames(dm[[j]])=m.names
    } 
  dm
    }
    
    dM.bin<-function(dirx,x,binm1,contm,m.names,w)  
    {dm<-NULL
     ref=apply(dirx!=0,1,sum)==0  #the reference group are all zeros. 
     if(is.null(w))
       means=apply(x[ref,c(binm1,contm)],2,mean)
     else
       means=apply(x[ref,c(binm1,contm)],2,weighted.mean,w[ref])

     for (j in 1:ncol(dirx)){  
       if(is.null(w))
         means1<-apply(x[dirx[,j]==1,c(binm1,contm)],2,mean)
       else
         means1<-apply(x[dirx[,j]==1,c(binm1,contm)],2,weighted.mean,w[dirx[,j]==1])
          dm<-rbind(dm,means1-means)
        } 
     colnames(dm)=m.names
     rownames(dm)=colnames(dirx)
     dm
    }
    
    match1<-function(vec1,vec2)
  {vec<-rep(F,length(vec1))
    for (i in 1:length(vec2))
     vec= vec | vec1==vec2[i]
    }
    
    if(is.null(catm))
      multi=jointm
    else if(is.null(jointm))
      multi=catm
    else {temp1<-catm
    temp2<-jointm
    temp1[[1]]=catm[[1]]+jointm[[1]]
    temp2[[1]]<-NULL
    multi=append(temp1,temp2)} 
    listm=list(single=c(contm,binm),multi=multi)
    
    if (is.null(multi))                      #allm list all mediators
    {tempm<-multi
    tempm[[1]]<-NULL}
    else  tempm<-NULL
    allm<-unique(c(contm,binm,unlist(tempm)))
    
    nonmissing<-apply(cbind(y,x,dirx),1,anymissing)
    x<-x[nonmissing,]
    y<-data.frame(y[nonmissing,])
    colnames(y)<-ynames
    pred<-data.frame(dirx[nonmissing,])
    colnames(pred)<-pred_names
    w<-w[nonmissing]

    binm1<-binm
    if(!is.null(catm))
    {for (i in 2:(catm$n+1))
      binm1<-c(binm1,catm[[i]])}
    
    n1<-dim(x)[1]
    m.names<-xnames[c(binm1,contm)]
    if(binpred){
      dm<-dM.bin(pred,x,binm1, contm, m.names,w)
      diag.temp<-abs(dm)
#      diag.temp1<-diag.temp
      distmgivenx=NULL
    }
    else{
     distmgivenx<-exp.m.given.x(x,pred,binm1,contm,df,w)
     dm<-dM(distmgivenx,pred,binm1,contm,df,margin,m.names) #draw ms conditional on x
     diag.temp<-NULL
#     diag.temp1<-NULL
     for (k in 1:ncol(pred))
       if(is.null(w))
         {diag.temp<-rbind(diag.temp, apply(dm[[k]],2,mean))
#         diag.temp1<-rbind(diag.temp1, apply(abs(dm[[k]]),2,mean))
          }
       else
         {diag.temp<-rbind(diag.temp, apply(dm[[k]],2,weighted.mean,w))
#         diag.temp1<-rbind(diag.temp1, apply(abs(dm[[k]]),2,weighted.mean,w))
         }
#     diag.temp=abs(diag.temp)
    }
    diag1<-rep(1,ncol(x))
    deltaM<-apply(abs(diag.temp),2,max)
    diag1[c(binm1,contm)]=1/deltaM
    x.star<-data.matrix(x)%*%diag(diag1)
    x2<-cbind(x.star,pred)
    colnames(x2)<-c(xnames,pred_names)
    full.model<-vector("list",ncol(y))  #D items, with dth item the model for the dth response
    coef.model<-NULL #D rows, with each row the fitted coefficient
    coef.m<-NULL  #D rows, with each row the fitted coefficient for mediators
    denm<-NULL # D*K, each row is the direct effect for the K predictors

    for(d in 1:ncol(y)){
      if(is.null(w))
         full.model[[d]]<-suppressWarnings(cv.glmnet(x=data.matrix(x2), y=y[,d],standardize=F,family=data$family1[[d]],alpha=alpha,lambda=lambda))
      else
         full.model[[d]]<-suppressWarnings(cv.glmnet(x=data.matrix(x2), y=y[,d],standardize=F,family=data$family1[[d]],weights=w,alpha=alpha,lambda=lambda))
      if(data$y_type[d]!=4)
          coef.temp<-coef(full.model[[d]])[-1,]
       else  coef.temp<-coef(full.model[[d]])[1:length(coef(full.model[[d]])),]  #survival outcome does not fit an intercept
       coef.model<-cbind(coef.model, coef.temp) 
       coef.m<-cbind(coef.m,coef.temp[c(binm1,contm)])

      #get the direct effect
       denm<-rbind(denm,tail(coef.model[,d],ncol(pred))/margin)
    }
    
    colnames(denm)<-pred_names
    rownames(denm)<-colnames(y)
    rownames(coef.m)=m.names
    colnames(coef.m)=colnames(y)
    rownames(coef.model)=c(xnames,pred_names)
    colnames(coef.model)=colnames(y)
    names(full.model)<-colnames(y)

    a<-list(denm=denm,model=list(full.model=full.model,family=family1,surv=surv,biny=biny),
            coef.m=coef.m,coef.model=coef.model,deltaM=deltaM,dm=dm,data=data,margin=margin,
            diag.temp=diag.temp,distmgivenx=distmgivenx,binm1=binm1,contm=contm,df=df,
            margin=margin,m.names=m.names,data=data)
    class(a)<-"med.big"
    return(a)
  }
  
print.med.big<-function(x,...,pred.new=NULL,w.new=NULL,print=T)
{dM<-function(distmgivenx,dirx,binm1,contm,df,margin,m.names)  
{means<-NULL
if(df>1)
{z<-NULL
for(i in 1:ncol(dirx))
  z<-cbind(z,ns(dirx[,i],df=df))}
else
  z<-dirx

if(!is.null(binm1))
  for (i in 1:length(binm1))
    means<-cbind(means,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z)))
  if(!is.null(contm))
    for (i in (length(binm1)+1):length(c(binm1,contm)))
      means<-cbind(means,predict(distmgivenx$models[[i]],newdata=data.frame(z)))
    
    dm<-vector("list", ncol(dirx))
    names(dm)<-colnames(dirx)
    for (j in 1:ncol(dirx)){  
      dirx1<-dirx
      dirx1[,j]=dirx[,j]+margin
      
      if(df>1)
      {z<-NULL
      for(i in 1:ncol(dirx))
        z<-cbind(z,ns(dirx1[,i],df=df))}
      else
        z<-dirx1
      
      means1<-NULL
      if(!is.null(binm1))
        for (i in 1:length(binm1))
          means1<-cbind(means1,predict(distmgivenx$models[[i]],type = "response",newdata=data.frame(z)))
      if(!is.null(contm))
        for (i in (length(binm1)+1):length(c(binm1,contm)))
          means1<-cbind(means1,predict(distmgivenx$models[[i]],newdata=data.frame(z)))
      dm[[j]]<-means1-means
      colnames(dm[[j]])=m.names
    } 
    dm
}

 results<-vector("list",ncol(x$data$y))
 names(results)<-colnames(x$data$y)
 if(is.null(pred.new))
   {diag.temp=x$diag.temp
    dm<-x$dm}
 else
  {dm<-dM(x$distmgivenx,pred.new,x$binm1,x$contm,x$df,x$margin,x$m.names) #draw ms conditional on x
   diag.temp<-NULL
   for (k in 1:ncol(pred.new))
    if(is.null(w.new))
     diag.temp<-rbind(diag.temp, apply(dm[[k]],2,mean))
    else
     diag.temp<-rbind(diag.temp, apply(dm[[k]],2,weighted.mean,w.new))
  }

 if(x$data$binpred)
   {dm=vector("list",ncol(x$data$dirx))
    for (k in 1:ncol(x$data$dirx))
     dm[[k]]=x$dm[k,]
   }
 
 for(d in 1:ncol(x$data$y))
 {if(x$data$binpred)
  {results[[d]]<-x$denm[d,]
   temp.coef<-x$coef.m[,d]/x$deltaM
   if(length(temp.coef)==1)
     results[[d]]<-rbind(results[[d]],t(x$dm)*temp.coef)
   else
     results[[d]]<-rbind(results[[d]],t(x$dm%*%diag(temp.coef)))
  }
  else
  {results[[d]]<-x$denm[d,]/x$margin
   temp.coef<-x$coef.m[,d]/x$deltaM/x$margin
   if(length(temp.coef)==1)
     results[[d]]<-rbind(results[[d]],t(diag.temp)*temp.coef)
   else
     results[[d]]<-rbind(results[[d]],t(diag.temp%*%diag(temp.coef)))
  }
 results[[d]]<-rbind(results[[d]],apply(results[[d]],2,sum))
 colnames(results[[d]])=colnames(x$data$dirx)
 rownames(results[[d]])=c("de",rownames(x$coef.m),"te")
 }
 re=list(results=results,dm=dm)
 if(print)
  print(re$results)
 return(re)
}




mma.big<-function(data=NULL,x=data$x, y=data$y,pred=data$dirx, mediator=NULL, binm=data$binm,contm=data$contm,catm=data$catm,
                  jointm=data$jointm,margin=1,df=1,binref=NULL,catref=NULL,predref=NULL,alpha=1,alpha1=0.01,alpha2=0.01,
                  family1=data$family1,n2=50,w=rep(1,nrow(x)),
                  refy=NULL,pred.new=NULL,binpred=data$binpred,type=NULL,w.new=NULL,lambda=exp(seq(log(0.001), log(5), length.out=15)))
{if(is.null(data)){
  mediator=unique(c(contm,mediator,binm,catm,jointm))
  data=data.org.big(x,y,pred,mediator,contmed=contm,binmed=binm,binref=binref,catmed=catm,
                    catref=catref,jointm=mediator,family1=family1,predref=predref,
                    alpha=alpha,alpha1=alpha1,alpha2=alpha2, testtype=1, w=w)
  x=data$x
  y=data$y
  pred=data$dirx
  binm=data$binm
  contm=data$contm
  catm=data$catm
  jointm=data$jointm
  family1=data$family1
  binpred=data$binpred
}
biny=data$y_type==2
surv=data$y_type==4
if (sum(surv)>0 & is.null(type))
  type="risk"

if (is.null(c(binm,contm,catm)))
  stop("Error: no potential mediator is specified")

xnames<-colnames(x)
pred_names<-colnames(pred)

if(is.null(catm))
{multi=jointm
name1<-NULL                       #added in the new program
if (!is.null(multi))              #added in the new program, in case that multi is NULL
  name1<-paste("j",1:multi[[1]],sep="")}
else if(is.null(jointm))
{multi=catm
name1<-NULL
for (i in 2:(catm[[1]]+1))
  name1<-c(name1,colnames(x)[multi[[i]][1]])}
else {temp1<-catm
temp2<-jointm
temp1[[1]]=catm[[1]]+jointm[[1]]
temp2[[1]]<-NULL
multi=append(temp1,temp2)
name1<-NULL
for (i in 2:(catm[[1]]+1))
  name1<-c(name1,colnames(x)[multi[[i]][1]])
name1<-c(name1,paste("j",1:jointm[[1]],sep=""))} 
listm=list(single=c(contm,binm),multi=multi)

D=ncol(y)
K=ncol(pred)
binm1<-binm
if(!is.null(catm))
{for (i in 2:(catm$n+1))
  binm1<-c(binm1,catm[[i]])}
J=length(c(binm1,contm))
dm<-vector("list",K)
coef<-vector("list",D)
deltaM<-matrix(0,n2,J)
colnames(deltaM)=xnames[c(binm1,contm)]
names(dm)=colnames(pred)
names(coef)=colnames(y)
bootresults=coef

results<-med.big(data=data,margin=margin,df=df,type=type,w=w,alpha=alpha,lambda=lambda)

for (i in 1:n2)
{boots<-sample(1:nrow(x),replace=T, prob=w)
 data1=data
 data1$x<-data.frame(x[boots,])
 data1$y<-data.frame(y[boots,])
 data1$dirx<-data.frame(pred[boots,])
 if(is.null(w))
   w1=NULL
 else
   w1=w[boots]
 temp<-med.big(data=data1,margin=margin,df=df,type=type,w=w,alpha=alpha,lambda=lambda) #added to the new codel, change the seed to make different results
 temp1<-print(temp,pred.new=pred.new,w.new=w.new,print=F)
 for (k in 1:K)
  dm[[k]]<-rbind(dm[[k]],temp1$dm[[k]])
 for (d in 1:D)
 {coef[[d]]=rbind(coef[[d]],t(temp$coef.m[,d]))
  deltaM[i,]=temp$deltaM
  bootresults[[d]]=cbind(bootresults[[d]],temp1$results[[d]])}
print(i)
}

a<-list(dm=dm,coef=coef,deltaM=deltaM,bootresults=bootresults,results=results,data=data,pred.new=pred.new,w.new=w.new,margin=margin)
class(a)<-"mma.big"
return(a)
}


print.mma.big<-function(x,...)
{print(x$result)
}


summary.mma.big<-function(object,...,alpha=0.05,plot=TRUE,RE=FALSE,quant=T,ball.use=T)
{sqr.dist<-function(vec1,vec2)
 {mean((vec1-vec2)^2,na.rm=T)}

 bound.ball<-function(mat1,mat2,alpha)
 {n1<-ncol(mat1)
  n2=nrow(mat2)
  upbd=NULL
  lwbd=NULL
  for(i in 1:n1)
  {temp.t<-i%%n1
   temp.z<-(1:ncol(mat2))%%n1==temp.t
   temp.m<-mat2[,temp.z]
   temp.d<-apply(temp.m[-n2,],2,sqr.dist,mat1[-n2,i])
   temp.ball<-temp.d<quantile(temp.d,1-alpha,na.rm=T) #temp.rank<=length(temp.dist)*(1-alpha)
   upbd<-c(upbd,apply(as.matrix(temp.m[,temp.ball]),1,max,na.rm=T))
   lwbd<-c(lwbd,apply(as.matrix(temp.m[,temp.ball]),1,min,na.rm=T))
  }

  return(cbind(upbd,lwbd)) 
 }

a1<-alpha/2
a2<-1-a1
b1<-qnorm(a1)
b2<-qnorm(a2)
x<-object
ny<-ncol(x$data$y)
nx<-ncol(x$data$dirx)
temp1<-x$bootresults
temp.0<-print(x$results,print=F)
temp2<-temp1
for(l in 1:ny)
  temp2[[l]]=temp2[[l]]%*%diag(1/temp2[[l]][nrow(temp2[[l]]),])
temp<-vector("list",ny)
names(temp)<-colnames(x$data$y)
temp.4<-vector("list",ny)                #save relative effects
names(temp.4)<-colnames(x$data$y)

for(l in 1:ny){
  temp.1<-NULL
  for (j in 1:nrow(temp1[[l]]))
    temp.1<-rbind(temp.1,matrix(temp1[[l]][j,],nrow=nx))
  temp[[l]]<-rbind(est=as.vector(temp.0$results[[l]]),
                   mean=as.vector(t(matrix(apply(temp.1,1,mean,na.rm=T),nrow=nx))),
                   sd=as.vector(t(matrix(apply(temp.1,1,sd,na.rm=T),nrow=nx))),
                   upbd_q=as.vector(t(matrix(apply(temp.1,1,quantile,a2,na.rm=T),nrow=nx))), 
                   lwbd_q=as.vector(t(matrix(apply(temp.1,1,quantile,a1,na.rm=T),nrow=nx))))
  temp[[l]]=rbind(temp[[l]],upbd=temp[[l]][1,]+b2*temp[[l]][3,])
  temp[[l]]=rbind(temp[[l]],lwbd=temp[[l]][1,]+b1*temp[[l]][3,])
  temp.3=temp.0$results[[l]]
  temp.2=bound.ball(temp.3,temp1[[l]],alpha)
  temp[[l]]=rbind(temp[[l]],upbd_b=temp.2[,1])
  temp[[l]]=rbind(temp[[l]],lwbd_b=temp.2[,2])
  colnames(temp[[l]])=paste(rep(colnames(x$data$dirx),each=nrow(temp1[[l]])),rownames(temp1[[l]]),sep=".")
}

for(l in 1:ny){
  temp.1<-NULL
  for (j in 1:nrow(temp2[[l]]))
    temp.1<-rbind(temp.1,matrix(temp2[[l]][j,],nrow=nx))
  if(nx>1)
    est=as.vector(temp.0$results[[l]]%*%diag(1/temp.0$results[[l]][nrow(temp.0$results[[l]]),]))
  else
    est=temp.0$results[[l]][,1]/temp.0$results[[l]][nrow(temp.0$results[[l]]),1]
  temp.4[[l]]<-rbind(est=est,
                   mean=as.vector(t(matrix(apply(temp.1,1,mean,na.rm=T),nrow=nx))),
                   sd=as.vector(t(matrix(apply(temp.1,1,sd,na.rm=T),nrow=nx))),
                   upbd_q=as.vector(t(matrix(apply(temp.1,1,quantile,a2,na.rm=T),nrow=nx))), 
                   lwbd_q=as.vector(t(matrix(apply(temp.1,1,quantile,a1,na.rm=T),nrow=nx))))
  temp.4[[l]]=rbind(temp.4[[l]],upbd=temp.4[[l]][1,]+b2*temp.4[[l]][3,])
  temp.4[[l]]=rbind(temp.4[[l]],lwbd=temp.4[[l]][1,]+b1*temp.4[[l]][3,])
  colnames(temp.4[[l]])=paste(rep(colnames(x$data$dirx),each=nrow(temp2[[l]])),rownames(temp2[[l]]),sep=".")
}

results=list(result=temp,re=temp.4,alpha=alpha,plot=plot,obj=x,RE=RE,quant=quant,nx=nx,ny=ny,ball.use=ball.use)
class(results)<-"summary.mma.big"
results
}

print.summary.mma.big<-function(x,...,digit=3)
{temp=x$result
 temp.4=x$re
 RE=x$RE
 plot=x$plot
 quant=x$quant
 nx=x$nx
 ny=x$ny
 ball.use=x$ball.use
 results<-x
if(RE)
  {cat("The relative effects:\n")
   temp.name<-names(results$re)
   for (l in 1:length(results$re))
   {cat(paste("\nFor the response variable",temp.name[l],"\n",sep=","))
    print(round(results$re[[l]],digit))}
  }
  else  
   {cat("\nThe mediaiton effects:\n")
    temp.name<-names(results$result)
   for (l in 1:length(results$result))
   {cat(paste("\nFor the response variable",temp.name[l],"\n",sep=","))
    print(round(results$result[[l]],digit))}
   }

J<-ncol(temp.4[[1]])/nx
pred.names<-colnames(x$obj$data$dirx)
if(plot)
  if(RE)
    for (m in 1:ny)
     for (l in 1:nx)
      {re<-temp.4[[m]][1,((l-1)*J+1):((l*J)-1)]
       if(quant)
       {re<-temp.4[[m]][2,((l-1)*J+1):((l*J)-1)]
        upper<-temp.4[[m]][4,((l-1)*J+1):((l*J)-1)]
        lower<-temp.4[[m]][5,((l-1)*J+1):((l*J)-1)]}
       else
       {upper<-temp.4[[m]][6,((l-1)*J+1):((l*J)-1)]
        lower<-temp.4[[m]][7,((l-1)*J+1):((l*J)-1)]}
      d<-order(re)
      name1<-colnames(temp.4[[m]][,((l-1)*J+1):((l*J)-1)])
      par(mfrow=c(1,1),mar=c(1,6,1,1),oma=c(3,2,2,4))
      bp <- barplot2(re[d], horiz = TRUE, main=paste("Relative Effects on y",m," on Predictor ",pred.names[l],sep=""), 
                     names.arg=name1[d],plot.ci = TRUE, ci.u = upper[d], ci.l = lower[d],
                     cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(c(upper,lower)),
                     col = rainbow(length(d), start = 3/6, end = 4/6))
      }
else
  for (m in 1:ny)
    for (l in 1:nx)
    {re<-temp[[m]][1,((l-1)*J+1):((l*J)-1)]
     est1=temp[[m]][1,l*J]
    if(ball.use)
    {upper<-temp[[m]][8,((l-1)*J+1):((l*J)-1)]
     lower<-temp[[m]][9,((l-1)*J+1):((l*J)-1)]
     upper1=temp[[m]][8,l*J]
     lower1=temp[[m]][9,l*J]
    }
    if(quant)
    {re<-temp[[m]][2,((l-1)*J+1):((l*J)-1)]
     upper<-temp[[m]][4,((l-1)*J+1):((l*J)-1)]
     lower<-temp[[m]][5,((l-1)*J+1):((l*J)-1)]
     est1=temp[[m]][2,l*J]
     upper1=temp[[m]][4,l*J]
     lower1=temp[[m]][5,l*J]}
    else
    {upper<-temp[[m]][6,((l-1)*J+1):((l*J)-1)]
     lower<-temp[[m]][7,((l-1)*J+1):((l*J)-1)]
     upper1=temp[[m]][6,l*J]
     lower1=temp[[m]][7,l*J]}
    d<-order(re)
    name1<-colnames(temp[[m]][,((l-1)*J+1):((l*J)-1)])
    par(mfrow=c(1,1),mar=c(1,6,1,1),oma=c(3,2,2,4))
    bp <- barplot2(c(re[d],est1), horiz = TRUE, main=paste("Relative Effects on y",m," on Predictor ",pred.names[l],sep=""), 
                   names.arg=c(name1[d],"te"),plot.ci = TRUE, ci.u = c(upper[d],upper1), ci.l = c(lower[d],lower1),
                   cex.names=0.9,beside=FALSE,cex.axis=0.9,las=1,xlim=range(c(upper,lower,upper1,lower1)),
                   col = rainbow(length(d)+1, start = 3/6, end = 4/6))
    }
}


plot.mma.big<-function(x,vari,...)
{D<-ncol(x$data$y)
 n<-nrow(x$data$y)
 K<-ncol(x$data$dirx)
 ynames=colnames(x$data$y)
 prednames=colnames(x$data$dirx)
 mnames<-colnames(x$coef[[1]])
 if(is.character(vari))
   vari=grep(vari,mnames)
 a=length(vari)
 if(a==0)
   return ("Select a mediator for vari")
 deltam=x$deltaM[,vari]
 coef3<-NULL
 for (d in 1:D)
 {coef1=x$coef[[d]][,vari]
  coef2=coef1/deltam
  if(d==1)
    par(fig=c((d-1)/D,d/D,0.5,1))
  else
    par(fig=c((d-1)/D,d/D,0.5,1),new=T)
  boxplot(coef2,names=mnames[vari],main=ynames[d],ylab=paste("coef for", mnames[vari])) 
  coef3<-cbind(coef3,coef2)
 }
 dm=vector("list",K)  
 for(k in 1:K)
  if(x$data$binpred)
  {dm1<-x$dm[[k]][,vari]
   par(fig=c((k-1)/K,k/K,0,0.47),new=T)
    boxplot(dm1,names=mnames[vari],main=prednames[k],ylab="delta M") 
  }
  else
  {if(is.null(x$pred.new))
    pred=x$data$dirx[,k]
   else 
    pred=x$pred.new[,k]
   z=order(pred)
   min1=NULL
   max1=NULL
   mean1=NULL
   for (i in 1:a)
    {dm1=matrix(x$dm[[k]][,vari[i]]/x$margin,nrow=n)
     min1=cbind(min1,apply(dm1,1,quantile,0.025,na.rm=T))
     max1=cbind(max1,apply(dm1,1,quantile,0.975,na.rm=T))
     mean1=cbind(mean1,apply(dm1,1,mean,na.rm=T))
     dm[[k]]<-rbind(dm[[k]],dm1)
    }
   range1=c(min(min1),max(max1))
   par(fig=c((k-1)/K,k/K,0,0.5),new=T)
   plot(stats::lowess(data.frame(pred[z],mean1[z,1])),type="l",col=1,xlab=prednames[k],ylab="delta M",ylim=range1)
   lines(stats::lowess(data.frame(pred[z],min1[z,1])),lty=2,col=1)
   lines(stats::lowess(data.frame(pred[z],max1[z,1])),lty=2,col=1)
   if(length(vari)>1)
     for (i in 2:length(vari))
     { lines(stats::lowess(data.frame(pred[z],mean1[z,i])),lty=1,col=i)
       lines(stats::lowess(data.frame(pred[z],min1[z,i])),lty=2,col=i)
       lines(stats::lowess(data.frame(pred[z],max1[z,i])),lty=2,col=i)
     }
  }

me<-vector("list",D)  #store the mediation effects
names(me)=ynames
for (d in 1:D){
  me[[d]]=vector("list",K)
  names(me[[d]])=prednames}

if(!x$data$binpred)
 for(d in 1:D)
   for (k in 1:K)
   {if(is.null(x$pred.new))
     pred=x$data$dirx[,k]
    else 
     pred=x$pred.new[,k]
    z=order(pred)
    min2=NULL
    max2=NULL
    mean2=NULL
    for (i in 1:a)
    {dm2=dm[[k]][((i-1)*n+1):(i*n),]%*%diag(coef3[,(d-1)*a+i])
     min2=cbind(min2,apply(dm2,1,quantile,0.025,na.rm=T))
     max2=cbind(max2,apply(dm2,1,quantile,0.975,na.rm=T))
     mean2=cbind(mean2,apply(dm2,1,mean,na.rm=T))
     me[[d]][[k]]<-rbind(me[[d]][[k]],dm2)
   }
   range1=c(min(min2),max(max2))
   if(d==1 & k==1)
     par(fig=c(0,1/K,(D-1)/D,1))
   else
     par(fig=c((k-1)/K,k/K,(D-d)/D,(D-d+1)/D),new=T)
   plot(stats::lowess(data.frame(pred[z],mean2[z,1])),type="l",col=1,xlab=prednames[k],
        ylab=paste("ie",ynames[d],sep=":"),ylim=range1)
   lines(stats::lowess(data.frame(pred[z],min2[z,1])),lty=2,col=1)
   lines(stats::lowess(data.frame(pred[z],max2[z,1])),lty=2,col=1)
   if(length(vari)>1)
     for (i in 2:length(vari))
     { lines(stats::lowess(data.frame(pred[z],mean2[z,i])),lty=1,col=i)
       lines(stats::lowess(data.frame(pred[z],min2[z,i])),lty=2,col=i)
       lines(stats::lowess(data.frame(pred[z],max2[z,i])),lty=2,col=i)
     }
   }
 
invisible(me)
}


joint.effect<-function(object,vari,alpha=0.05)
{
a1<-alpha/2
a2<-1-a1
b1<-qnorm(a1)
b2<-qnorm(a2)
x<-object
ny<-ncol(x$data$y)
nx<-ncol(x$data$dirx)
temp1<-x$bootresults
temp.0<-print(x$results,print=F)
temp2<-temp1
for(l in 1:ny)
  temp2[[l]]=temp2[[l]]%*%diag(1/temp2[[l]][nrow(temp2[[l]]),])
temp<-vector("list",ny)
names(temp)<-colnames(x$data$y)
temp.4<-vector("list",ny)                #save relative effects
names(temp.4)<-colnames(x$data$y)

for(l in 1:ny){
  temp.namea=rownames(temp1[[l]])
  temp.a=unique(unlist(lapply(vari,grep,temp.namea)))
  if (length(temp.a)==1)
    return("Please select multiple mediators.")
  else
    temp.1<-apply(temp1[[l]][temp.a,],2,sum)
  temp.1<-matrix(temp.1,nrow=nx)
  temp.nameb=rownames(temp.0$results[[l]])
  temp.b=unique(unlist(lapply(vari,grep,temp.nameb)))
  temp[[l]]<-rbind(est=apply(as.matrix(temp.0$results[[l]][temp.b,]),2,sum),
                   mean=apply(temp.1,1,mean,na.rm=T),
                   sd=apply(temp.1,1,sd,na.rm=T),
                   upbd_q=apply(temp.1,1,quantile,a2,na.rm=T), 
                   lwbd_q=apply(temp.1,1,quantile,a1,na.rm=T))
  temp[[l]]=rbind(temp[[l]],upbd=temp[[l]][1,]+b2*temp[[l]][3,])
  temp[[l]]=rbind(temp[[l]],lwbd=temp[[l]][1,]+b1*temp[[l]][3,])
  colnames(temp[[l]])=colnames(x$data$dirx)
#  temp.3=temp.0$results[[l]]
#  temp.2=bound.ball(temp.3,temp1[[l]],alpha)
#  temp[[l]]=rbind(temp[[l]],upbd_b=temp.2[,1])
#  temp[[l]]=rbind(temp[[l]],lwbd_b=temp.2[,2])
#  colnames(temp[[l]])=paste(rep(colnames(x$data$dirx),each=nrow(temp1[[l]])),rownames(temp1[[l]]),sep=".")
}

for(l in 1:ny){
  if (length(temp.a)==1)
    temp.1=temp2[[l]][temp.a,]
  else
    temp.1<-apply(temp2[[l]][temp.a,],2,sum)
  temp.1<-matrix(temp.1,nrow=nx)
  if(nx>1)
    est=apply((temp.0$results[[l]]%*%diag(1/temp.0$results[[l]][nrow(temp.0$results[[l]]),]))[temp.b,],2,sum)
  else
    est=sum(temp.0$results[[l]][temp.b,1]/temp.0$results[[l]][nrow(temp.0$results[[l]]),1])
  temp.4[[l]]<-rbind(est=est,
                     mean=apply(temp.1,1,mean,na.rm=T),
                     sd=apply(temp.1,1,sd,na.rm=T),
                     upbd_q=apply(temp.1,1,quantile,a2,na.rm=T), 
                     lwbd_q=apply(temp.1,1,quantile,a1,na.rm=T))
  temp.4[[l]]=rbind(temp.4[[l]],upbd=temp.4[[l]][1,]+b2*temp.4[[l]][3,])
  temp.4[[l]]=rbind(temp.4[[l]],lwbd=temp.4[[l]][1,]+b1*temp.4[[l]][3,])
  colnames(temp.4[[l]])=colnames(x$data$dirx)
}
#cat("The joint mediation effect and relative effect of")
list(variables=temp.namea[temp.a],effect=temp,relative.effect=temp.4)
}











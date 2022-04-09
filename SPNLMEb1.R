#author Maritza Marquez
#March 2022

#######################################################################
   ### SAEMIX code - Pspline (beta-HCG hormone concentration) ###
     #################### Simulated Data ##################### 
#######################################################################
 
library(saemix)
library(sitar)
library(ggplot2)
require(gridExtra)
par(mar=rep(2,4))
library(colourpicker)

data <- read.csv("data/datasimul.csv", header=T)
data$X <- NULL
 

############################# 
     ### NLME Model ###
#############################
 
#Longitudinal data structure used by SAEMIX
data.nlme <- saemixData(name.data=data,header=TRUE,sep=" ",na=NA,
                    name.group=c("id"),name.predictors=c("time"),
                    name.response=c("y"),
                    units=list(x="days",y="log(beta-HCG)"),name.X="time")
 

#Longitudinal data model
model.log<-function(psi,id,xidep) { 
   tim<-xidep[,1] 
   
   a<-psi[id,1]
   b<-psi[id,2]
   c<-psi[id,3]
   
  ypred<-a/(1+exp(-(tim-b)/c))
   
  return(ypred)
 }
 
#Saemix model 
model.nlme<-saemixModel(model=model.log,description="Logitudinal Model",
                           psi0=matrix(c(4,15,7,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("a","b","c"))),
                           transform.par=c(0,0,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
                           covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
                                                   byrow=TRUE))
 
#RUNS
K1 = 500
K2 = 100
iterations = 1:(K1+K2+1)
end = K1+K2
 
options.nlme<-list(seed=139546,map=FALSE,fim=FALSE,ll.is=TRUE,
                  nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),
                  displayProgress = TRUE,save.graphs = FALSE,nbiter.burn = 0)
 
nlme.fit <- saemix(model.nlme,data.nlme,options.nlme) 
 
plot(nlme.fit, plot.type="individual.fit")
 
 
#Plot NLME
library(lattice)
 
xyplot(y + fitted(nlme.fit) ~ time | id, data = data,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130), which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l"), distribute.type = TRUE, col.line=c("blue","#DB3B1B"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME", sub="id - abnormal and normal group",
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME"),
                        col=c("blue","#DB3B1B")))


############################## 
   ### Bsplines Models ###
#############################

library(splines)
 
#B-spline bases
bspline <- function(x,xl,xr,ndx,bdeg){
   dx <- (xr-xl)/ndx
   knots <- seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
   B <- spline.des(knots,x,bdeg+1,0*x,outer.ok=TRUE)$design
   output <- list(knots=knots,B=B)
   return(output)
 }
 
#Mixed model for B-splines
mixed.model<-function(x,xl,xr,ndx,bdeg,pord,type="Eilers"){
   Bbasis=bspline(x,xl,xr,ndx,bdeg)
   B <- Bbasis$B
   m=ncol(B)
   D=diff(diag(m),differences=pord)
   
   if(type=="Eilers"){
     Z <- B%*%t(D)%*%solve(D%*%t(D))
   }else if(type=="SVD"){ print("SVD method")
     P.svd=svd(t(D)%*%D)
     U=(P.svd$u)[,1:(m-pord)]
     d=(P.svd$d)[1:(m-pord)]
     Delta=diag(1/sqrt(d))
     Z=B%*%U%*%Delta
   }
   X=NULL
   for(i in 0:(pord-1)){
     X=cbind(X,x^i)
   }
   output <- list(X=X,Z=Z)
   return(output)
 }


#Bspline model with 5 knots
ndx = 5 # knots numbers
bdeg  = 3 # polynomial degree 
pord = 2 # penalty order
type = "SVD"
 
x <- data$time
mm <-mixed.model(x,min(x)-0.01,max(x)+0.01,ndx,bdeg,pord,type=type)
 
X <- mm$X
Z <- mm$Z
 
aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(data)[2]
 
RES <- cbind(X,Z)
 
data.full.5k<-cbind(data,Z)
colnames(data.full.5k)[(aD+1):(aD+aZ)] <- paste0("Z",1:aZ)
 
data.bs.5k<- saemixData(name.data = data.full.5k, header = TRUE, 
                                 sep = " ", na = NA, 
                                 name.predictors = c("time",paste0("Z",1:aZ)), 
                                 name.response = c("y"), 
                                 units = list(x = "", y = "", 
                                              covariates = c("-", "-")), 
                                 name.X = c("time"))
 
#Data 5-knots bspline model to run in monolix  
write.csv(data.full.5k, file = "databs5k.csv", row.names = FALSE)
 
plot(data.bs.5k, plot.type="data")
 
modelbsplines.5k <- function(psi,id,x) { 
  RES <- x
   
  b1<- psi[id, 1]
  b2<- psi[id, 2]
  b3<- psi[id, 3]
  b4<- psi[id, 4]
  b5<- psi[id, 5]
  b6<- psi[id, 6]
   
  ypred <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + RES[,7]*b6
   
  return(ypred)
 }
 
model.bs.5k<-saemixModel(model = modelbsplines.5k, 
                                description = "Bsplines model 5 knots", 
                                #type='structural',
                                psi0 = matrix(c(rep(0,6)), ncol = 6, byrow = TRUE,
                                              dimnames = list(NULL, c("b1","b2","b3","b4","b5","b6"))),
                                transform.par = c(rep(0,6)),
                                fixed.estim=c(rep(1,6)),
                                covariance.model = cbind(rbind(diag(6))))
 
options.bs<-list(seed=39546,map=FALSE,fim=FALSE,ll.is=TRUE,nb.chains = 5,
                    nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),
                    displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
 

bs.fit.5k <- saemix(model.bs.5k, data.bs.5k, options.bs)
 
plot(bs.fit.5k,plot.type="individual.fit")
 

#Plot Bspline 5 knots
library(lattice)
 
xyplot(y + fitted(nlme.fit) + fitted(bs.fit.5k) ~ time | id, data = data.full.5k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#E6CE17"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Bspline 5 knots (SVD)", sub="id - abnormal and normal group",
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME", "Bspline 5 knots"),
                col=c("blue","#DB3B1B","#E6CE17")))
 
#Bspline model with 8 knots
ndx = 8 
bdeg  = 3
pord = 2
type = "SVD"
 
x <- data$time
mm <- mixed.model(x,min(x)-0.01,max(x)+0.01,ndx,bdeg,pord,type=type)
 
X <- mm$X
Z <- mm$Z
 
aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(data)[2]
 
RES <- cbind(X,Z)
data.full.8k <-cbind(data,Z)
colnames(data.full.8k)[(aD+1):(aD+aZ)] <- paste0("Z",1:9)
 
#Data 8-knots bspline model to run in monolix   
write.csv(data.full.8k, file = "databs8k.csv", row.names = FALSE)
 
data.bs.8k <- saemixData(name.data = data.full.8k, header = TRUE, 
                                 name.group=c("id"),
                                 sep = " ", na = NA, 
                                 name.predictors = c("time",paste0("Z",1:aZ)), 
                                 name.response = c("y"), 
                                 units = list(x = "days", y = "log(beta-HCG)", 
                                              covariates = c("-", "-")), 
                                 name.X = c("time"))

plot(data.bs.8k, plot.type="data")
 

modelbsplines.8k <- function(psi,id,x) { 
   RES <- x
   
   b1<- psi[id, 1]
   b2<- psi[id, 2]
   b3<- psi[id, 3]
   b4<- psi[id, 4]
   b5<- psi[id, 5]
   b6<- psi[id, 6]
   b7<- psi[id, 7]
   b8<- psi[id, 8]
   b9<- psi[id, 9]
   
   ypred <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + RES[,10]*b9
   
   return(ypred)
 }
 
 
model.bs.8k<-saemixModel(model = modelbsplines.8k, 
                                description = "Bsplines model 8 knots", 
                                #type="structural",
                                psi0 = matrix(c(rep(0,9)), ncol = 9, byrow = TRUE, dimnames = list(NULL, c("b1","b2","b3","b4","b5","b6","b7","b8","b9"))),
                                transform.par = c(rep(0,9)),
                                fixed.estim=c(rep(1,9)),
                                covariance.model = cbind(rbind(diag(9))))
 
 
options.bs<-list(seed=39546,map=FALSE,fim=FALSE,ll.is=TRUE,nb.chains = 5,
                    nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),
                    displayProgress=TRUE,save.graphs=FALSE,nbiter.burn = 0)
 
 
bs.fit.8k <- saemix(model.bs.8k, data.bs.8k, options.bs)
 
plot(bs.fit.8k,plot.type="individual.fit")
 
 
#Plot Bsplines 8 knots
 
library(lattice)
 
xyplot(y + fitted(nlme.fit) + fitted(bs.fit.8k) ~ time | id, data = data.full.8k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#1AD98D"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Bspline 8 knots (SVD)", sub="id - abnormal and normal group",                               
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME", "Bspline 8 knots"),
                        col=c("blue","#DB3B1B","#1AD98D")))
       

#Bspline model with 10 knots 
ndx = 10
bdeg  = 3 
pord = 2
type = "SVD"
 
x <- data$time
mm <- mixed.model(x,min(x)-0.01,max(x)+0.01,ndx,bdeg,pord,type=type)
 
X <- mm$X
Z <- mm$Z
 
aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(data)[2]
 
RES <- cbind(X,Z)
data.full.10k <-cbind(data,Z)
colnames(data.full.10k)[(aD+1):(aD+aZ)] <- paste0("Z",1:aZ)
 
#Data 10-knots bspline model to run in monolix 
write.csv(data.full.10k, file = "databs10k.csv", row.names = FALSE)
 
data.bs.10k <- saemixData(name.data = data.full.10k, header = TRUE, 
                                 name.group=c("id"),
                                 sep = " ", na = NA, 
                                 name.predictors = c("time",paste0("Z",1:aZ)), 
                                 name.response = c("y"), 
                                 units = list(x = "days", y = "log(beta-HCG)", 
                                              covariates = c("-", "-")), 
                                 name.X = c("time"))
 
plot(data.bs.10k, plot.type="data")
 
 
modelbsplines.10k <- function(psi,id,x) { 
   RES <- x
   
   b1<- psi[id, 1]
   b2<- psi[id, 2]
   b3<- psi[id, 3]
   b4<- psi[id, 4]
   b5<- psi[id, 5]
   b6<- psi[id, 6]
   b7<- psi[id, 7]
   b8<- psi[id, 8]
   b9<- psi[id, 9]
   b10<-psi[id,10]
   b11<-psi[id,11]
   
   ypred <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + RES[,10]*b9 + RES[,11]*b10 + RES[,12]*b11
   
   return(ypred)
 }
 
model.bs.10k<-saemixModel(model = modelbsplines.10k, 
                                description = "Bspline model 10 knots", 
                                #type="structural",
                                psi0 = matrix(c(rep(0,11)), ncol = 11, byrow = TRUE, dimnames = list(NULL, c("b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11"))),
                                transform.par = c(rep(0,11)),
                                fixed.estim=c(rep(1,11)),
                                covariance.model = cbind(rbind(diag(11))))
 
 
options.bs<-list(seed=39546,map=FALSE,fim=FALSE,ll.is=TRUE,nb.chains = 5,
                    nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),
                    displayProgress=TRUE,save.graphs=FALSE,nbiter.burn = 0)
 
 
bs.fit.10k <- saemix(model.bs.10k, data.bs.10k,options.bs)
 
plot(bs.fit.10k,plot.type="individual.fit")
 
 
#Plot Bsplines 10 knots
library(lattice)
 
xyplot(y + fitted(nlme.fit) + fitted(bs.fit.10k) ~ time | id, data = data.full.10k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#3E7008"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Bspline 10 knots (SVD)", sub="id - abnormal and normal group",
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME", "Bspline 10 knots"),
                       col=c("blue","#DB3B1B","#3E7008")))
                               
 
### Bsplines Plots

xyplot(y + fitted(nlme.fit) + fitted(bs.fit.5k) + fitted(bs.fit.8k) + fitted(bs.fit.10k) ~ time | id, data = data.full.10k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#E6CE17","#1AD98D","#3E7008"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)",main="NLME - Bsplines (SVD)", sub="id - abnormal and normal group",
        auto.key =list(points=FALSE,lines = FALSE, space = "right",text=c("Observations","NLME","Bspline 5 knots","Bspline 8 knots","Bspline 10 knots"),
                       col=c("blue","#DB3B1B","#E6CE17","#1AD98D","#3E7008")))

 

################################### 
  ### Semiparametric Models ###
      ## NLME + Bspline ##
###################################

# Semiparametric 5 knots
 
modelsp.5k <-function(psi, id, x){ 
   RES <- x
   
   tim<-RES[,1]  
   a<-psi[id,1]
   b<-psi[id,2]
   c<-psi[id,3]
   
   ypred<-a/(1+exp(-(tim-b)/c))
   
   return(ypred)
   
   ##
   
   b1<- psi[id, 4]
   b2<- psi[id, 5]
   b3<- psi[id, 6]
   b4<- psi[id, 7]
   b5<- psi[id, 8]
   b6<- psi[id, 9]
   
   
   yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + RES[,7]*b6
   resp <- ypred + yspline 
   
   return(resp)
 }
 
ndx = 5
bdeg  = 3
pord = 2
type = "SVD"
 
x <- data$time
mm <-mixed.model(x,min(x)-0.01,max(x)+0.01,ndx,bdeg,pord,type=type)
 
X <- mm$X
Z <- mm$Z
 
aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(data)[2]
 
RES <- cbind(X,Z)
datasp.full.5k <-cbind(data,Z)
colnames(datasp.full.5k)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))
 
data.sp.5k <- saemixData(name.data = datasp.full.5k, header = TRUE, 
                                 sep = " ", na = NA, 
                                 name.predictors = c("time",paste0("Z",1:aZ)), 
                                 name.response = c("y"), 
                                 units = list(x = "days", y = "log(beta-HCG)", 
                                              covariates = c("-")), 
                                 name.X = c("time"))
 
#Data 5-knots semiparametric model to run in monolix 
write.csv(datasp.full.5k, file = "datasp5k.csv", row.names = FALSE)
 
plot(data.sp.5k, plot.type="data")
 
model.sp.5k <-saemixModel(model= modelsp.5k, 
                                      #type='structural',
                                      description = "Semiparametric model 5 knots", 
                                      psi0 = matrix(c(coef(nlme.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                                    dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                                      transform.par = c(rep(0,3),rep(0,aZ)),
                                      fixed.estim=c(rep(1,3+aZ)),
                                      covariance.model = cbind(rbind(diag(9))))
 
 
options.sp<-list(seed=39546,map=FALSE,fim=FALSE,ll.is=TRUE,nb.chains = 5,
                    nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0.9,
                    displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
 
 
sp.fit.5k <- saemix(model.sp.5k, data.sp.5k, options.sp)
 
plot(sp.fit.5k,plot.type="individual.fit")
 
 
#Plot semiparametric model with 5 knots
library(lattice)
 
xyplot(y + fitted(nlme.fit) + fitted(sp.fit.5k) ~ time | id, data = data.full.5k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#E6CE17"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Semiparametric (SVD)", sub="id - abnormal and normal group",
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME", "Semiparametric 5 knots"),
                       col=c("blue","#DB3B1B","#E6CE17")))

 
### Semiparametric model with 8 knots
modelsp.8k <-function(psi, id, x){ 
   RES <- x
   
   tim<-RES[,1]
   a<-psi[id,1]
   b<-psi[id,2]
   c<-psi[id,3]
   ypred<-a/(1+exp(-(tim-b)/c))
   
   return(ypred)
   
   ##
   
   b1<- psi[id, 4]
   b2<- psi[id, 5]
   b3<- psi[id, 6]
   b4<- psi[id, 7]
   b5<- psi[id, 8]
   b6<- psi[id, 9]
   b7<- psi[id, 10]
   b8<- psi[id, 11]
   b9<- psi[id, 12]
   
   yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + RES[,10]*b9
   resp <- ypred + yspline
   
   return(resp)
 }
 
ndx = 8
bdeg  = 3
pord = 2
type = "SVD"
 
x <- data$time
mm <-mixed.model(x,min(x)-0.01,max(x)+0.01,ndx,bdeg,pord,type=type)
 
X <- mm$X
Z <- mm$Z
 
aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(data)[2]
 
RES <- cbind(X,Z)
datasp.full.8k <-cbind(data,Z)
colnames(datasp.full.8k)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))
 
#Data 5-knots semiparametric model to run in monolix 
write.csv(datasp.full.8k, file = "datasp8k.csv", row.names = FALSE)
 
 
data.sp.8k <- saemixData(name.data = data.full.8k, header = TRUE, 
                                 sep = " ", na = NA, 
                                 name.predictors = c("time",paste0("Z",1:aZ)), 
                                 name.response = c("y"), 
                                 units = list(x = "days", y = "log(beta-HCG)", 
                                              covariates = c("-")), 
                                 name.X = c("time"))
 
 
plot(data.sp.8k, plot.type="data")
 
model.sp.8k<-saemixModel(model=modelsp.8k, 
                                      #type='structural',
                                      description = "Semiparametric model 8 knots", 
                                      psi0 = matrix(c(coef(nlme.fit)$fixed,rep(0,9)), ncol = 3+aZ, byrow = TRUE,
                                                    dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9"))),
                                      transform.par = c(rep(0,3),rep(0,aZ)),
                                      fixed.estim=c(rep(1,3+aZ)),
                                      covariance.model = cbind(rbind(diag(12))))
 
 
options.sp<-list(seed=39546,map=FALSE,fim=FALSE,ll.is=TRUE,nb.chains = 5,
                    nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0.9,
                    displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
 
 
sp.fit.8k<- saemix(model.sp.8k, data.sp.8k, options.sp)
 
plot(sp.fit.8k,plot.type="individual.fit")
 
 
#Plot semiparametric model 8 knots 
library(lattice)
 
xyplot(y + fitted(nlme.fit) + fitted(sp.fit.8k) ~ time | id, data = datasp.full.8k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#1AD98D"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Semiparametric (SVD)", sub="id - abnormal and normal group",
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME", "Semiparametric 8 knots"),
                       col=c("blue","#DB3B1B","#1AD98D")))


# Semiparametric model with 10 knots
modelsp.10k<-function(psi, id, x){
   RES <- x
   
   tim<-RES[,1]
   a<-psi[id,1]
   b<-psi[id,2]
   c<-psi[id,3]
   
   ypred<-a/(1+exp(-(tim-b)/c))
   return(ypred)
   
   ##
   
   b1<- psi[id, 4]
   b2<- psi[id, 5]
   b3<- psi[id, 6]
   b4<- psi[id, 7]
   b5<- psi[id, 8]
   b6<- psi[id, 9]
   b7<- psi[id, 10]
   b8<- psi[id, 11]
   b9<- psi[id, 12]
   b10<-psi[id,13]
   b11<-psi[id,14]
   
   yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + RES[,10]*b9 + RES[,11]*b10 + RES[,12]*b11
   resp <- ypred + yspline
   
   return(resp)
 }
 
ndx = 10
bdeg  = 3
pord = 2
type = "SVD"
 
x <- data$time
mm <-mixed.model(x,min(x)-0.01,max(x)+0.01,ndx,bdeg,pord,type=type)
 
X <- mm$X
Z <- mm$Z
 
aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(data)[2]
 
RES <- cbind(X,Z)
datasp.full.10k <-cbind(data,Z)
colnames(datasp.full.10k)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))
 
 
#Data 5-knots semiparametric model to run in monolix 
write.csv(datasp.full.10k, file = "datasp10k.csv", row.names = FALSE)
 
data.sp.10k <- saemixData(name.data = data.full.10k, header = TRUE,
                                 sep = " ", na = NA,
                                 name.predictors = c("time",paste0("Z",1:aZ)),
                                 name.response = c("y"),
                                 units = list(x = "days", y = "log(beta-HCG)",
                                              covariates = c("-")),
                                 name.X = c("time"))
 
 
plot(data.sp.10k, plot.type="data")
 
model.sp.10k<-saemixModel(model=modelsp.10k,
                                      #type='structural',
                                      description = "Semiparametric model 10 knots",
                                      psi0 = matrix(c(coef(nlme.fit)$fixed,rep(0,11)), ncol = 3+aZ, byrow = TRUE,
                                                    dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11"))),
                                      transform.par = c(rep(0,3),rep(0,aZ)),
                                      fixed.estim=c(rep(1,3+aZ)),
                                      covariance.model = cbind(rbind(diag(14))))
 
 
options.sp<-list(seed=39546,map=FALSE,fim=FALSE,ll.is=TRUE,nb.chains = 5,
                    nbiter.mcmc = c(2,2,2,0), nbiter.saemix = c(K1,K2),nbiter.sa=0.9,
                    displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)
 
 
sp.fit.10k <- saemix(model.sp.10k, data.sp.10k, options.sp)
 
plot(sp.fit.10k,plot.type="individual.fit")
 
#Plot semiparametric model 10 knots
library(lattice)
 
 xyplot(y + fitted(nlme.fit) + fitted(sp.fit.10k) ~ time | id, data = datasp.full.10k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#3E7008"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Semiparametric (SVD)", sub="id - abnormal and normal group",
        auto.key = list(points=FALSE, lines = FALSE, space = "right", text= c("Observations","NLME", "Semiparametric 10 knots"),
                       col=c("blue","#DB3B1B","#3E7008")))

 
#### Semiparametrics Models Plots 
 
xyplot(y + fitted(nlme.fit) + fitted(sp.fit.5k) +  fitted(sp.fit.8k) + fitted(sp.fit.10k) ~ time | id, data = data.full.10k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#E6CE17","#1AD98D","#3E7008"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME - Semiparametric (SVD)", sub="id - abnormal and normal group",
        auto.key =list(points=FALSE,lines = FALSE, space = "right",text=c("Observations","NLME","Semiparametric 5 knots","Semiparametric 8 knots","Semiparametric 10 knots"),
                       col=c("blue","#DB3B1B","#E6CE17","#1AD98D","#3E7008")))
 

#### Graficos Comparativos
 
xyplot(y + fitted(nlme.fit) + fitted(bs.fit.5k) +  fitted(sp.fit.5k) ~ time | id, data = data.full.5k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#3E7008","#1AD98D"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME-Bspline-Semiparametric", sub="id - abnormal and normal group",
        auto.key =list(points=FALSE,lines = FALSE, space = "right",text=c("Observations","NLME","Bsplines 5 knots","Semipararmetric 5 knots"),
                       col=c("blue","#DB3B1B","#3E7008","#1AD98D")))

 
xyplot(y + fitted(nlme.fit) + fitted(bs.fit.8k) +  fitted(sp.fit.8k) ~ time | id, data = data.full.8k,
        subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                   which(id==130),which(id==136), which(id==158), which(id==169)),
        strip = FALSE, aspect = "xy", pch = 16, grid = TRUE,
        type = c("p","l","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#3E7008","#1AD98D"), lwd=2,
        xlab = "days", ylab = "log(beta-HCG)", main="NLME-Bspline-Semiparametric", sub="id - abnormal and normal group",
        auto.key =list(points=FALSE,lines = FALSE, space = "right",text=c("Observations","NLME","Bsplines 8 knots","Semiparametric 8 knots"),
                       col=c("blue","#DB3B1B","#3E7008","#1AD98D")))
 
 
xyplot(y + fitted(nlme.fit) + fitted(bs.fit.10k) +  fitted(sp.fit.10k) ~ time | id, data = data.full.10k,
       subset = c(which(id==14), which(id==62), which(id==65), which(id==120),
                  which(id==130),which(id==136), which(id==158), which(id==169)),
       strip = FALSE, aspect = "xy", pch = 16, grid = TRUE, #layout=c(4,4),
       type = c("p","l","l","l"), distribute.type = TRUE,col.line=c("blue","#DB3B1B","#3E7008","#1AD98D"), lwd=2,
       xlab = "days", ylab = "log(beta-HCG)", main="NLME-Bspline-Semiparametric", sub="id - abnormal and normal group",
       auto.key =list(points=FALSE,lines = FALSE, space = "right",text=c("Observations","NLME","Bsplines 10 knots","Semiparametric 10 knots"),
                      col=c("blue","#DB3B1B","#3E7008","#1AD98D")))
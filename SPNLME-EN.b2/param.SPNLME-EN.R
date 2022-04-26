#author Maritza Márquez
#septiembre 2021

#######################################################################
     ### Codigo SAEMIX - Normal Pregnancy Semiparametric Model ###
         ######### Using different parameterization ########
#######################################################################

library(saemix)
library(ggplot2)
library(sitar)
par(mar=rep(2,4))

dataE <- read.table("C:/Users/asus/Dropbox/Mi PC (DESKTOP-INEI7US)/Documents/DOCTORADO MATEMATICAS/PROYECTO TESIS/R - Tesis/SAEMIX-EN Semiparametric Model/data/bHCG_D.txt", header=T)
colnames(dataE)[1:4]<-c("id","time","y","grupo")
dataE$Time7 <-NULL
dataE$D <- NULL
dataEN <- dataE[dataE$grupo == 0,]  #data normal pregnancy


### Longitudinal model normal pregnancies ###

nlmeEN.data <-saemixData(name.data=dataEN,header=TRUE,sep=" ",na=NA,
                         name.group=c("id"),name.predictors=c("time"),
                         name.response=c("y"),
                         units=list(x="days",y="log(beta-HCG)-EN"),name.X="time")


# Function for longitudinal data (NLME)

logEN<-function(psi,id,xidep) { 
  tim<-xidep[,1] 
  
  a<-psi[id,1]
  b<-psi[id,2]
  c<-psi[id,3]
  
  ypred<-a/(1+exp(-(tim-b)/c))
  
  return(ypred)
}


### Saemix model for longitudinal data for normal pregnancies

nlmeEN.model<-saemixModel(model=logEN,description="Modelo longitudinal EN",
                          psi0=matrix(c(4,15,7,0,0,0),ncol=3,byrow=TRUE, dimnames=list(NULL, c("a","b","c"))),
                          transform.par=c(0,0,0),omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
                          covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,
                                                  byrow=TRUE))

##RUNS
K1 = 500
K2 = 100
c1 = 5
c2 = 5
c3 = 5
iterations = 1:(K1+K2+1)
end = K1+K2


options_log <-list(seed=139546,map=F,fim=F,ll.is=T,
                   nbiter.mcmc = c(c1,c2,c3), nbiter.saemix = c(K1,K2),
                   displayProgress = TRUE, save.graphs = FALSE, nbiter.burn = 0)

nlmeEN.fit <- saemix(nlmeEN.model, nlmeEN.data, options_log) 


plot(nlmeEN.fit, plot.type = "individual.fit")


# B-splines basis for semiparametric model

library(splines)

# Obtain the B-spline Basis
bbase <- function(x, xl = min(x), xr = max(x), ndx = ndx, bdeg = bdeg) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- spline.des(knots, x, ord = bdeg + 1, outer.ok = TRUE)$design
  res <- list(B = B, knots = knots, bdeg = bdeg)
  class(res) <- "bbase"
  res
}

# We represent the mixed model for the B-splines

MMSP.basis <- function (x, xl, xr, ndx, bdeg, pord, param = param) {
  Bb <- bbase(x,xl,xr,ndx,bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B)
  n <- nrow(B)
  D <- diff(diag(m), differences=pord)
  P.svd <- svd(crossprod(D))
  U.Z <- (P.svd$u)[,1:(m-pord)]
  d <- (P.svd$d)[1:(m-pord)]
  Z <- B%*%U.Z
  U.X <- NULL
  if(param == 1) {
    U.X = ((P.svd$u)[,-(1:(m-pord))])
    X = B%*%U.X
  } else if (param == 2){
    X <- NULL
    for(i in 0:(pord-1)){
      X <- cbind(X,x^i)
    }
    U.X <- NULL
    for(i in 0:(pord-1)){
      U.X <- cbind(U.X, knots[-c((1:(bdeg - 1)),(length(knots)- (bdeg - 1) + 1):length(knots))]^i)
    }
    # Note that B%*%U.X == X
  } else if(param == 3) {
    U.X <- NULL
    for(i in 0:(pord-1)){
      U.X <- cbind(U.X, knots[-c((1:(bdeg - 1)),(length(knots)- (bdeg - 1) + 1):length(knots))]^i)
    }
    X <- B%*%U.X
  } else if(param == 4) {
    X <- B%*%((P.svd$u)[,-(1:(m-pord))])
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X <- ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
  } else if(param == 5) {
    U.Z <- (t(D)%*%solve(D%*%t(D)))
    Z <- B%*%U.Z
    X <- B%*%((P.svd$u)[,-(1:(m-pord))])
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X <- ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf					
  } else if(param == 6) {
    U.Z <- t(D)
    Z <- B%*%U.Z		
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X <- ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf			
  } else if(param == 7) {
    U.Z <- contr.sum(ncol(B))%*%contr.sum(ncol(B)-1)
    Z <- B%*%U.Z
    X <- NULL
    for(i in 0:(pord-1)){
      X <- cbind(X,x^i)
    }		
  }	
  list(Bb = Bb, X = X, Z = Z, d = d, m = m, D = D, U.X = U.X, U.Z = U.Z)
}

semipar.model <-function(psi, id, x){ 
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

# Construct B-spline basis for semiparametric model 
# (with different parametrization; argument param)


############## Models with different parametrization #################

ndx = 5
bdeg  = 3
pord = 2
x <- dataEN$time


######### Model 1 ###########

MSP1 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 1)
names(MSP1)

X <- MSP1$X
Z <- MSP1$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN1.full <-cbind(dataEN,Z)
colnames(dataEN1.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN1.full.data <- saemixData(name.data = dataEN1.full, header = TRUE, 
                               sep = " ", na = NA, 
                               name.predictors = c("time",paste0("Z",1:aZ)), 
                               name.response = c("y"), 
                               units = list(x = "days", y = "log(beta-HCG)", 
                                            covariates = c("-")), 
                               name.X = c("time"))

plot(dataEN1.full.data, plot.type="data")


SPEN1.saemix <-saemixModel(model=semipar.model,
                               description = "Semiparametric Model param = 1", 
                               psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                             dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                               transform.par = c(rep(0,3),rep(0,aZ)),
                               fixed.estim=c(rep(1,3+aZ)),
                               covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN1.fit <- saemix(SPEN1.saemix, dataEN1.full.data, options_SP)

plot(SPEN1.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=1.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN1.fit) ~ time | id, data = dataEN1.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()


######### Model 2 ###########

MSP2 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 2)
names(MSP2)

X <- MSP2$X
Z <- MSP2$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN2.full <-cbind(dataEN,Z)
colnames(dataEN2.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN2.full.data <- saemixData(name.data = dataEN2.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEN2.full.data, plot.type="data")


SPEN2.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 2", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN2.fit <- saemix(SPEN2.saemix, dataEN2.full.data, options_SP)

plot(SPEN2.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=2.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN2.fit) ~ time | id, data = dataEN2.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()

######### Model 3 ###########

MSP3 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 3)
names(MSP3)

X <- MSP3$X
Z <- MSP3$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN3.full <-cbind(dataEN,Z)
colnames(dataEN3.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN3.full.data <- saemixData(name.data = dataEN3.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEN3.full.data, plot.type="data")


SPEN3.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 3", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN3.fit <- saemix(SPEN3.saemix, dataEN3.full.data, options_SP)

plot(SPEN3.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=3.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN3.fit) ~ time | id, data = dataEN3.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()


######### Model 4 ###########

MSP4 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 4)
names(MSP4)

X <- MSP4$X
Z <- MSP4$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN4.full <-cbind(dataEN,Z)
colnames(dataEN4.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN4.full.data <- saemixData(name.data = dataEN4.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEN4.full.data, plot.type="data")


SPEN4.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 4", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN4.fit <- saemix(SPEN4.saemix, dataEN4.full.data, options_SP)

plot(SPEN4.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=4.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN4.fit) ~ time | id, data = dataEN4.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()

######### Model 5 ###########

MSP5 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 5)
names(MSP5)

X <- MSP5$X
Z <- MSP5$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN5.full <-cbind(dataEN,Z)
colnames(dataEN5.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN5.full.data <- saemixData(name.data = dataEN5.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEN5.full.data, plot.type="data")


SPEN5.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 5", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN5.fit <- saemix(SPEN5.saemix, dataEN5.full.data, options_SP)

plot(SPEN5.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=5.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN5.fit) ~ time | id, data = dataEN5.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()

######### Model 6 ###########

MSP6 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 6)
names(MSP6)

X <- MSP6$X
Z <- MSP6$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN6.full <-cbind(dataEN,Z)
colnames(dataEN6.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN6.full.data <- saemixData(name.data = dataEN6.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEN6.full.data, plot.type="data")


SPEN6.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 6", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN6.fit <- saemix(SPEN6.saemix, dataEN6.full.data, options_SP)

plot(SPEN6.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=6.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN6.fit) ~ time | id, data = dataEN6.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()

######### Model 7 ###########

MSP7 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 7)
names(MSP6)

X <- MSP7$X
Z <- MSP7$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN7.full <-cbind(dataEN,Z)
colnames(dataEN7.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN7.full.data <- saemixData(name.data = dataEN7.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEN7.full.data, plot.type="data")


SPEN7.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 7", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(9))))


##RUNS
K3 = 500
K4 = 100
c4 = 5
c5 = 5
c6 = 5
iterations = 1:(K3+K4+1)
end = K3+K4


options_SP<-list(seed=139546,map=F,fim=F,ll.is=T,nb.chains = 5,
                 nbiter.mcmc = c(c4,c5,c6), nbiter.saemix = c(K3,K4),nbiter.sa=0.9,
                 displayProgress=TRUE,save.graphs=FALSE,nbiter.burn =0)


SPEN7.fit <- saemix(SPEN7.saemix, dataEN7.full.data, options_SP)

plot(SPEN7.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-param=7.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN7.fit) ~ time | id, data = dataEN7.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()

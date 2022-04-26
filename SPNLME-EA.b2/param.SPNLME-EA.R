#author Maritza Márquez
#septiembre 2021

#######################################################################
     ### Codigo SAEMIX - Abnormal Pregnancy Semiparametric Model ###
           ######### Using different parameterization  ########
#######################################################################

library(saemix)
library(ggplot2)
library(sitar)
par(mar=rep(2,4))

dataE <- read.table("C:/Users/asus/Dropbox/Mi PC (DESKTOP-INEI7US)/Documents/DOCTORADO MATEMATICAS/PROYECTO TESIS/R - Tesis/SAEMIX-EN Semiparametric Model/data/bHCG_D.txt", header=T)
colnames(dataE)[1:4]<-c("id","time","y","grupo")
dataE$Time7 <-NULL
dataE$D <- NULL
dataEA <- dataE[dataE$grupo == 1,]  #data abnormal pregnancy


### Longitudinal model abnormal pregnancies ###

nlmeEA.data <-saemixData(name.data=dataEA,header=TRUE,sep=" ",na=NA,
                         name.group=c("id"),name.predictors=c("time"),
                         name.response=c("y"),
                         units=list(x="days",y="log(beta-HCG)-EA"),name.X="time")


# Function for longitudinal data (NLME)

logEA<-function(psi,id,xidep) { 
  tim<-xidep[,1] 
  
  a<-psi[id,1]
  b<-psi[id,2]
  c<-psi[id,3]
  
  ypred<-a/(1+exp(-(tim-b)/c))
  
  return(ypred)
}


### Saemix model for longitudinal data for normal pregnancies

nlmeEA.model<-saemixModel(model=logEA,description="Modelo longitudinal EA",
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

nlmeEA.fit <- saemix(nlmeEA.model, nlmeEA.data, options_log) 


plot(nlmeEA.fit, plot.type = "individual.fit")


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
x <- dataEA$time


######### Model 1 ###########

MSP1 <- MMSP.basis(x, min(x), max(x), ndx = ndx, bdeg = bdeg, pord = pord, param = 1)
names(MSP1)

X <- MSP1$X
Z <- MSP1$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA1.full <-cbind(dataEA,Z)
colnames(dataEA1.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA1.full.data <- saemixData(name.data = dataEA1.full, header = TRUE, 
                               sep = " ", na = NA, 
                               name.predictors = c("time",paste0("Z",1:aZ)), 
                               name.response = c("y"), 
                               units = list(x = "days", y = "log(beta-HCG)", 
                                            covariates = c("-")), 
                               name.X = c("time"))

plot(dataEA1.full.data, plot.type="data")


SPEA1.saemix <-saemixModel(model=semipar.model,
                               description = "Semiparametric Model param = 1", 
                               psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA1.fit <- saemix(SPEA1.saemix, dataEA1.full.data, options_SP)

plot(SPEA1.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=1.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA1.fit) ~ time | id, data = dataEA1.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                  which(id==156),which(id==164), which(id==168), which(id==169)),
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
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA2.full <-cbind(dataEA,Z)
colnames(dataEA2.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA2.full.data <- saemixData(name.data = dataEA2.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEA2.full.data, plot.type="data")


SPEA2.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 2", 
                           psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA2.fit <- saemix(SPEA2.saemix, dataEA2.full.data, options_SP)

plot(SPEA2.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=2.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA2.fit) ~ time | id, data = dataEA2.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                  which(id==156),which(id==164), which(id==168), which(id==169)),
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
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA3.full <-cbind(dataEA,Z)
colnames(dataEA3.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA3.full.data <- saemixData(name.data = dataEA3.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEA3.full.data, plot.type="data")


SPEA3.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 3", 
                           psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA3.fit <- saemix(SPEA3.saemix, dataEA3.full.data, options_SP)

plot(SPEA3.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=3.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA3.fit) ~ time | id, data = dataEA3.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                  which(id==156),which(id==164), which(id==168), which(id==169)),
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
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA4.full <-cbind(dataEA,Z)
colnames(dataEA4.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA4.full.data <- saemixData(name.data = dataEA4.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEA4.full.data, plot.type="data")


SPEA4.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 4", 
                           psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA4.fit <- saemix(SPEA4.saemix, dataEA4.full.data, options_SP)

plot(SPEA4.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=4.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA4.fit) ~ time | id, data = dataEA4.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                 which(id==156),which(id==164), which(id==168), which(id==169)),
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
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA5.full <-cbind(dataEA,Z)
colnames(dataEA5.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA5.full.data <- saemixData(name.data = dataEA5.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEA5.full.data, plot.type="data")


SPEA5.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 5", 
                           psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA5.fit <- saemix(SPEA5.saemix, dataEA5.full.data, options_SP)

plot(SPEA5.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=5.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA5.fit) ~ time | id, data = dataEA5.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                  which(id==156),which(id==164), which(id==168), which(id==169)),
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
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA6.full <-cbind(dataEA,Z)
colnames(dataEA6.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA6.full.data <- saemixData(name.data = dataEA6.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEA6.full.data, plot.type="data")


SPEA6.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 6", 
                           psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA6.fit <- saemix(SPEA6.saemix, dataEA6.full.data, options_SP)

plot(SPEA6.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=6.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA6.fit) ~ time | id, data = dataEA6.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                  which(id==156),which(id==164), which(id==168), which(id==169)),
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
aD <- dim(dataEA)[2]

RES <- cbind(X,Z)
dataEA7.full <-cbind(dataEA,Z)
colnames(dataEA7.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEA7.full.data <- saemixData(name.data = dataEA7.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

plot(dataEA7.full.data, plot.type="data")


SPEA7.saemix <-saemixModel(model=semipar.model, 
                           description = "Semiparametric Model param = 7", 
                           psi0 = matrix(c(coef(nlmeEA.fit)$fixed,rep(0,6)), ncol = 3+aZ, byrow = TRUE,
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


SPEA7.fit <- saemix(SPEA7.saemix, dataEA7.full.data, options_SP)

plot(SPEA7.fit, plot.type="individual.fit")


## Gráfico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EA-param=7.pdf")

xyplot(y + fitted(nlmeEA.fit) + fitted(SPEA7.fit) ~ time | id, data = dataEA7.full,
       subset = c(which(id==125), which(id==127), which(id==143), which(id==144),
                  which(id==156),which(id==164), which(id==168), which(id==169)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME"),
                      col=c("blue","black","magenta")))

dev.off()

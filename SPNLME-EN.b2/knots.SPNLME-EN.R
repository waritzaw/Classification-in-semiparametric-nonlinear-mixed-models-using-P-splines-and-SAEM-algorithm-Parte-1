#author Maritza Marquez
#march 2022

#######################################################################
### Codigo SAEMIX - Normal Pregnancy (simulated data) ###
          ######### Using different nodes  ########
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
predict.bbase <- function(object, newx, deriv = 0) {
  B <- spline.des(knots = object$knots, x = newx, derivs = deriv, ord = object$bdeg + 1, outer.ok = TRUE)$design
  B
}
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


# Construct B-spline basis for semiparametric model 
# (with different nodes; argument = ndx)


############## Models with different nodes #################

bdeg  = 3
pord = 2
param = 3
x <- dataEN$time

######### Model 1 ###########

model.5n <-function(psi, id, x){ 
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

MSP5 <- MMSP.basis(x, min(x), max(x), ndx = 5, bdeg = bdeg, pord = pord, param = param)
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

# Data modelo semiparametrico EN 5 nodos a usar en monolix  
write.csv(dataEN5.full, file = "dataEN-MSP5.csv", row.names = FALSE)

plot(dataEN5.full.data, plot.type="data")


SPEN5.saemix <-saemixModel(model=model.5n,
                           description = "Semiparametric Model nodos = 5", 
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


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=5.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN5.fit) ~ time | id, data = dataEN5.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 5n"),
                      col=c("blue","black","magenta")))

dev.off()


######### Model 2 ###########

model.6n <-function(psi, id, x){ 
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
  
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + 
    RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7
  resp <- ypred + yspline
  
  return(resp)
}

MSP6 <- MMSP.basis(x, min(x), max(x), ndx = 6, bdeg = bdeg, pord = pord, param = param)
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

# Data modelo semiparametrico EN 6 nodos a usar en monolix  
write.csv(dataEN6.full, file = "dataEN-MSP6.csv", row.names = FALSE)

plot(dataEN6.full.data, plot.type="data")


SPEN6.saemix <-saemixModel(model=model.6n, 
                           description = "Semiparametric Model nodos = 6", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,7)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(10))))


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

## Gr?fico NLME & SMNLME (6 nodos)

library(lattice)

pdf("NLME_SMNLME-EN-nodo=6.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN6.fit) ~ time | id, data = dataEN6.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 6n"),
                      col=c("blue","black","magenta")))

dev.off()



######### Model 3 ###########

model.7n <-function(psi, id, x){ 
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
  
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + 
             RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8
  resp <- ypred + yspline
  
  return(resp)
}

MSP7 <- MMSP.basis(x, min(x), max(x), ndx = 7, bdeg = bdeg, pord = pord, param = param)
names(MSP7)

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

# Data modelo semiparametrico EN 7 nodos a usar en monolix  
write.csv(dataEN7.full, file = "dataEN-MSP7.csv", row.names = FALSE)

plot(dataEN7.full.data, plot.type="data")


SPEN7.saemix <-saemixModel(model=model.7n, 
                           description = "Semiparametric Model nodos = 7", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,8)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(11))))


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


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=7.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN7.fit) ~ time | id, data = dataEN7.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 7n"),
                      col=c("blue","black","magenta")))

dev.off()



######### Model 4 ###########

model.8n <-function(psi, id, x){ 
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
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + 
             RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + 
             RES[,10]*b9
  resp <- ypred + yspline
  
  return(resp)
}

MSP8 <- MMSP.basis(x, min(x), max(x), ndx = 8, bdeg = bdeg, pord = pord, param = param)
names(MSP8)

X <- MSP8$X
Z <- MSP8$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN8.full <-cbind(dataEN,Z)
colnames(dataEN8.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN8.full.data <- saemixData(name.data = dataEN8.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

# Data modelo semiparametrico EN 8 nodos a usar en monolix  
write.csv(dataEN8.full, file = "dataEN-MSP8.csv", row.names = FALSE)

plot(dataEN8.full.data, plot.type="data")


SPEN8.saemix <-saemixModel(model=model.8n, 
                           description = "Semiparametric Model nodos = 8", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,9)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(12))))


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


SPEN8.fit <- saemix(SPEN8.saemix, dataEN8.full.data, options_SP)

plot(SPEN8.fit, plot.type="individual.fit")


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=8.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN8.fit) ~ time | id, data = dataEN8.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 8n"),
                      col=c("blue","black","magenta")))

dev.off()


######### Model 5 ###########

model.9n <-function(psi, id, x){ 
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
  b10<-psi[id, 13]
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + 
             RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + RES[,10]*b9 + 
             RES[,11]*b10
  resp <- ypred + yspline
  
  return(resp)
}

MSP9 <- MMSP.basis(x, min(x), max(x), ndx = 9, bdeg = bdeg, pord = pord, param = param)
names(MSP9)

X <- MSP9$X
Z <- MSP9$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN9.full <-cbind(dataEN,Z)
colnames(dataEN9.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN9.full.data <- saemixData(name.data = dataEN9.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

# Data modelo semiparametrico EN 9 nodos a usar en monolix  
write.csv(dataEN9.full, file = "dataEN-MSP9.csv", row.names = FALSE)

plot(dataEN9.full.data, plot.type="data")


SPEN9.saemix <-saemixModel(model=model.9n, 
                           description = "Semiparametric Model nodos = 9", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,10)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(13))))


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


SPEN9.fit <- saemix(SPEN9.saemix, dataEN9.full.data, options_SP)

plot(SPEN9.fit, plot.type="individual.fit")


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=9.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN9.fit) ~ time | id, data = dataEN9.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 9n"),
                      col=c("blue","black","magenta")))

dev.off()


######### Model 6 ###########

model.10n <-function(psi, id, x){ 
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
  b10<-psi[id, 13]
  b11<-psi[id, 14]
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + 
    RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + 
    RES[,10]*b9 + RES[,11]*b10 + RES[,12]*b11
  resp <- ypred + yspline
  
  return(resp)
}

MSP10 <- MMSP.basis(x, min(x), max(x), ndx = 10, bdeg = bdeg, pord = pord, param = param)
names(MSP10)

X <- MSP10$X
Z <- MSP10$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN10.full <-cbind(dataEN,Z)
colnames(dataEN10.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN10.full.data <- saemixData(name.data = dataEN10.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

# Data modelo semiparametrico EN 10 nodos a usar en monolix  
write.csv(dataEN10.full, file = "dataEN-MSP10.csv", row.names = FALSE)

plot(dataEN10.full.data, plot.type="data")


SPEN10.saemix <-saemixModel(model=model.10n, 
                           description = "Semiparametric Model nodos = 10", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,11)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(14))))


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


SPEN10.fit <- saemix(SPEN10.saemix, dataEN10.full.data, options_SP)

plot(SPEN10.fit, plot.type="individual.fit")


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=10.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN10.fit) ~ time | id, data = dataEN10.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 10n"),
                      col=c("blue","black","magenta")))

dev.off()


######### Model 7 ###########

model.11n <-function(psi, id, x){ 
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
  b10<-psi[id, 13]
  b11<-psi[id, 14]
  b12<-psi[id, 15]
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + 
    RES[,6]*b5 + RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + 
    RES[,10]*b9 + RES[,11]*b10 +  RES[,12]*b11 + RES[,13]*b12
  resp <- ypred + yspline
  
  return(resp)
}

MSP11 <- MMSP.basis(x, min(x), max(x), ndx = 11, bdeg = bdeg, pord = pord, param = param)
names(MSP11)

X <- MSP11$X
Z <- MSP11$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN11.full <-cbind(dataEN,Z)
colnames(dataEN11.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN11.full.data <- saemixData(name.data = dataEN11.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

# Data modelo semiparametrico EN 11 nodos a usar en monolix  
write.csv(dataEN11.full, file = "dataEN-MSP11.csv", row.names = FALSE)

plot(dataEN11.full.data, plot.type="data")


SPEN11.saemix <-saemixModel(model=model.11n, 
                           description = "Semiparametric Model nodos = 11", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,12)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(15))))


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


SPEN11.fit <- saemix(SPEN11.saemix, dataEN11.full.data, options_SP)

plot(SPEN11.fit, plot.type="individual.fit")


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=11.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN11.fit) ~ time | id, data = dataEN11.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 11n"),
                      col=c("blue","black","magenta")))

dev.off()



######### Model 8 ###########

model.12n <-function(psi, id, x){ 
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
  b10<-psi[id, 13]
  b11<-psi[id, 14]
  b12<-psi[id,15]
  b13<-psi[id,16]
  
  yspline <- RES[,2]*b1 + RES[,3]*b2 + RES[,4]*b3 + RES[,5]*b4 + RES[,6]*b5 + 
    RES[,7]*b6 + RES[,8]*b7 + RES[,9]*b8 + RES[,10]*b9 + RES[,11]*b10 + 
    RES[,12]*b11 + RES[,13]*b12 + RES[,14]*b13
  resp <- ypred + yspline
  
  return(resp)
}

MSP12 <- MMSP.basis(x, min(x), max(x), ndx = 12, bdeg = bdeg, pord = pord, param = param)
names(MSP12)

X <- MSP12$X
Z <- MSP12$Z

aX <- dim(X)[2]
aZ <- dim(Z)[2]
aD <- dim(dataEN)[2]

RES <- cbind(X,Z)
dataEN12.full <-cbind(dataEN,Z)
colnames(dataEN12.full)[(aD+1):(aD+aZ)] <- paste0("Z",1:ncol(Z))

dataEN12.full.data <- saemixData(name.data = dataEN12.full, header = TRUE, 
                                sep = " ", na = NA, 
                                name.predictors = c("time",paste0("Z",1:aZ)), 
                                name.response = c("y"), 
                                units = list(x = "days", y = "log(beta-HCG)", 
                                             covariates = c("-")), 
                                name.X = c("time"))

# Data modelo semiparametrico EN 12 nodos a usar en monolix  
write.csv(dataEN12.full, file = "dataEN-MSP12.csv", row.names = FALSE)

plot(dataEN12.full.data, plot.type="data")


SPEN12.saemix <-saemixModel(model=model.12n, 
                           description = "Semiparametric Model nodos = 12", 
                           psi0 = matrix(c(coef(nlmeEN.fit)$fixed,rep(0,13)), ncol = 3+aZ, byrow = TRUE,
                                         dimnames = list(NULL, c("a","b","c","b1","b2","b3","b4","b5","b6","b7","b8","b9","b10","b11","b12","b13"))),
                           transform.par = c(rep(0,3),rep(0,aZ)),
                           fixed.estim=c(rep(1,3+aZ)),
                           covariance.model = cbind(rbind(diag(16))))


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


SPEN12.fit <- saemix(SPEN12.saemix, dataEN12.full.data, options_SP)

plot(SPEN12.fit, plot.type="individual.fit")


## Gr?fico NLME & SMNLME

library(lattice)

pdf("NLME_SMNLME-EN-nodo=12.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN12.fit) ~ time | id, data = dataEN12.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 12n"),
                      col=c("blue","black","magenta")))

dev.off()


## Gr?fico NLME & SMNLME (NODOS PARES)

library(lattice)

pdf("NLME_SMNLME-EN-nodo=par.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN6.fit) + fitted(SPEN8.fit) + fitted(SPEN10.fit) + fitted(SPEN12.fit) ~ time | id, data = dataEN12.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l","l","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta","green","red","darkblue"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 6n", "SMNLME 8n", "SMNLME 10n", "SMNLME 12n"),
                      col=c("blue","black","magenta","green","red","darkblue")))

dev.off()


## Gr?fico NLME & SMNLME (NODOS IMPARES)

library(lattice)

pdf("NLME_SMNLME-EN-nodo=imp.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN5.fit) + fitted(SPEN7.fit) + fitted(SPEN9.fit) + fitted(SPEN11.fit) ~ time | id, data = dataEN11.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l","l","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta","green","red", "darkblue"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 5n", "SMNLME 7n", "SMNLME 9n", "SMNLME 11n"),
                      col=c("blue","black","magenta","green","red", "darkblue")))

dev.off()


### Comparamos los modelos de 11 y 12 nodos ###
### (Parecen ser los mejores ajustes hasta ahora)###

## Gr?fico NLME & SMNLME (NODOS 11 & 12)

library(lattice)

pdf("NLME_SMNLME-EN-nodo=11-12.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN11.fit) + fitted(SPEN12.fit) ~ time | id, data = dataEN12.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l","l"), distribute.type = TRUE, col.line=c("blue","black","magenta","darkblue"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SMNLME 11n", "SMNLME 12n"),
                      col=c("blue","black","magenta","darkblue")))

dev.off()


## Gr?fico NLME & SMNLME (NODOS 5 al 12)

library(lattice)

pdf("NLME_SPNLME-EN-nodos=5-12.pdf")

xyplot(y + fitted(nlmeEN.fit) + fitted(SPEN5.fit) + fitted(SPEN8.fit) + fitted(SPEN11.fit) + fitted(SPEN12.fit) ~ time | id, data = dataEN12.full,
       subset = c(which(id==14), which(id==30), which(id==41), which(id==68),
                  which(id==80),which(id==91), which(id==120), which(id==122)),
       strip = TRUE, pch = 16, grid = TRUE,
       type = c("p","l","l","l","l","l"), distribute.type = TRUE, col.line=c("blue","red","black","darkblue", "green", "magenta"),
       xlab = "Time (days)", ylab = "log(beta-HCG)",layout=c(4,2),
       auto.key =list(points = FALSE,lines = FALSE, space = "right",
                      text = c("Observaciones","NLME","SPNLME 5n", "SPNLME 8n", "SPNLME 11n", "SPNLME 12n"),
                      col=c("blue","red","black","darkblue", "green", "magenta")))

dev.off()

# Estamos separando los datos de concentración de la hormona BHCG en embarazos normales 
# (EN) y embarazos anormales (EA)


# We separate the data set
dataEN<-data[data$grupo <= 0,] # normal pregnancies 
dataEA<-data[data$grupo >= 1,] # abnormal pregnancies

plotEN=ggplot(dataEN, aes(x = time, y = y, color = factor(grupo), group = id)) + 
  geom_point()+ 
  geom_line()+ 
  theme_linedraw(base_size = 11)  + 
  labs(colour = "Group") +
  theme(legend.position = "none") + 
  xlab("days") + ylab("log(beta-HCG)") +
  ggtitle("Normal Pregnancies") +
  theme(plot.title = element_text(hjust = 0.5)) 

plotEA=ggplot(dataEA, aes(x = time, y = y, color = factor(grupo), group = id)) + 
  geom_point(color ="#4FB6BD") + 
  geom_line(color ="#4FB6BD") + 
  theme_linedraw(base_size = 11)  + 
  labs(colour = "Group") +
  theme(legend.position = "none") + 
  xlab("days") + ylab("log(beta-HCG)") +
  ggtitle("Abnormal Pregnancies") +
  theme(plot.title = element_text(hjust = 0.5)) 

library(ggpubr)
pdf("plot_EN-EA.pdf")
ggarrange(plotEN,plotEA,ncol=2)
dev.off()



# Codigo para contruir bases spline 

library(splines)

#B-spline bases
bspline <- function(x,xl,xr,ndx,bdeg){
  dx <- (xr-xl)/ndx
  knots <- seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
  B <- spline.des(knots,x,bdeg+1,0*x,
                  outer.ok=TRUE)$design
  output <- list(knots=knots,B=B)
  return(output)
}

#Mixed model for B-splines
mixed.model<-function(x,xl,xr,ndx,bdeg,pord,
                      type="Eilers"){
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


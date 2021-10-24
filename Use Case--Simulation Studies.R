### This file consists of simulation codes used in the main paper ###
library(doParallel)


### Setting 1 ###
##Training
N=5000
J=200
set.seed(1)
theta.vec=runif(N, -2, 2)
theta.vec=theta.vec-mean(theta.vec)
beta.vec=runif(J, -2, 2)
dat=generate_missingness(N, J)

temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec);
prob = 1/(1+exp(-temp));
registerDoParallel(cores = 25)  #Use 25 CPU for parallel training
cl <- makeCluster(25, type="FORK")  
result <- foreach(i=1:2000) %dopar% {#Training for each of the 2000 Monte Carlo samples
  data = matrix(0, N, J);
  data[] = rbinom(N*J, 1, prob);
  data[is.na(dat)]=NA; #Keeps the missingness pattern to be the same
  re=Rasch_missing(data, beta.vec, theta.vec);
  return(re)
}
stopCluster(cl)  

##Inference
par_est1=unlist(result)
theta1=par_est1$the
beta1=par_est1$beta
par_sample=list(theta_sample=theta1, beta_sample=beta1)
da=generate_missingness(N=5000, J=200)


#Estimation accuracy in terms of MSE
mean((as.matrix(theta1)-rep(1, 2000)%*%t(theta.vec))^2)#MSE for theta hat is 0.06372664;
mean((as.matrix(beta1)-rep(1, 2000)%*%t(beta.vec))^2)#MSE for beta hat is 0.002799754;
M1=theta.vec%*%t(rep(1, 200))-rep(1, 5000)%*%t(beta.vec)
sumSquaredErr=0
th1=as.matrix(theta1)
bt1=as.matrix(beta1)
for (i in 1:2000){
  M_c=th1[i,]%*%t(rep(1, 200))-rep(1, 5000)%*%t(bt1[i, ])
  sumSquaredErr=sumSquaredErr+sum((M1-M_c)^2)
}
sumSquaredErr/(5000*200*2000)#MSE for M hat is 0.06652639


#Evaluate oracle variance
temp=matrix(0, nrow = 5000, ncol = 200)
for (i in 1:N){
  for (j in 1:J){
    temp[i,j] = theta[i]-beta[j]
  }
}
p_mat = 1/(1+exp(-temp))

vaTrue=p_mat*(1-p_mat)
vaTrue[is.na(da)==T]=0
vaTrue_i=apply(vaTrue, 1, sum)
vaTrue_j=apply(vaTrue, 2, sum)
vaTrue_=sum(vaTrue)

#Compute estimated asymptotic variance using sample estimates
a_est= apply(par_sample$theta_sample, 2, mean)
b_est= apply(par_sample$beta_sample, 2, mean)

#Derived variance approximate
MEst= a_est %*% t(rep(1, J)) - rep(1, N) %*% t(b_est)
va=matrix(0, nrow=N, ncol=J)
for (i in 1:N){
  for (j in 1:J){
    va[i,j]=exp(-MEst[i,j])/((1+exp(-MEst[i,j]))^2)
  }
}
va[is.na(da)==T]=0
va_i=apply(va, 1, sum)
va_j=apply(va, 2, sum)
va_=sum(va)
varM11=1/va_i[1]+1/va_j[1]
varTheta1=1/va_i[1]
varBeta1=1/va_j[1]


#Plot density curves
M11Est=par_sample$theta_sample[,1]-par_sample$beta_sample[,1]#Sample estimates of M11
pdf("density_plots_corrected.pdf", width=12, height=4)
par(mfrow=c(1,3),cex=0.7, mai=c(0.6,0.6,0.6,0.6))
hist(M11Est, breaks=seq(0, 2, 0.1),xlim=c(0,2),ylim=c(0,2.5),freq=F, col='cadetblue',axes=F, xlab='', ylab='Density', main='(a)') 
lines(seq(0,2,0.01), generate_norm_density(sigma=varM11^0.5,mu=mean(M11Est),0,2,0.01))
axis(2, at=seq(0,2.5,0.5), labels =T)
axis(1, at=seq(0,2,0.5), labels = T, cex.axis=1, las=0)

hist(par_sample$theta_sample[,1], breaks=seq(-2, 0.2, 0.1),xlim=c(-2,0.2),col='cadetblue',ylim=c(0,2.5),freq=F,axes=F, xlab='', ylab='Density', main='(b)') 
lines(seq(-2,0.2,0.01), generate_norm_density(sigma=varTheta1^0.5,mu=mean(par_sample$theta_sample[,1]),-2,0.2,0.01))
axis(2, at=seq(0,2.5,0.5), labels =T)
axis(1, at=seq(-2,0.2,0.5), labels = T, cex.axis=1, las=0);

hist(par_sample$beta_sample[,1], breaks=seq(-2.2, -1.6, 0.015),xlim=c(-2.2, -1.6),col='cadetblue',ylim=c(0,10),freq=F,axes = F, xlab='', ylab='Density', main='(c)') 
lines(seq(-2.2, -1.6,0.01), generate_norm_density(sigma=varBeta1^0.5,mu=mean(par_sample$beta_sample[,1]),-2.2, -1.6,0.01))
axis(2, at=seq(0, 10,2), labels =T)
axis(1, at=seq(-2.2, -1.6,0.2), labels = T, cex.axis=1, las=0);
dev.off()






























### Setting 2 ###

##Training
N2=10000
J2=400
set.seed(2)
theta.vec2=runif(N2, -2, 2)
theta.vec2=theta.vec2-mean(theta.vec2)
beta.vec2=runif(J2, -2, 2)
dat=generate_missingness(N2, J2)
temp = theta.vec2 %*% t(rep(1, J2)) - rep(1, N2) %*% t(beta.vec2);
prob = 1/(1+exp(-temp));
registerDoParallel(cores = 25)  
cl <- makeCluster(25, type="FORK")  
result <- foreach(i=1:2000) %dopar% {
  data = matrix(0, N2, J2);
  data[] = rbinom(N2*J2, 1, prob);
  data[is.na(dat)]=NA;
  re=Rasch_missing(data, beta.vec2, theta.vec2);
  return(re)
}
stopCluster(cl)  

par_est2=unlist(result)


##Inference
theta2=par_est2$the
beta2=par_est2$beta
par_sample2=list(theta_sample=theta2, beta_sample=beta2)
beta2=beta.vec2[,1]
theta2=theta.vec2[,1]
da2=generate_missingness(N=N2, J=J2)

#Estimation accuracy in terms of MSE
mean((as.matrix(par_sample2$theta_sample)-rep(1, 2000)%*%t(theta2))^2)#MSE for theta hat is 0.03132349
mean((as.matrix(par_sample2$beta_sample)-rep(1, 2000)%*%t(beta2))^2)#MSE for beta hat is 0.001348611
M1=theta2%*%t(rep(1, 400))-rep(1, 10000)%*%t(beta2)
sumSquaredErr2=0
th2=as.matrix(par_sample2$theta_sample)
bt2=as.matrix(par_sample2$beta_sample)
for (i in 1:2000){
  M_c=th2[i,]%*%t(rep(1, 400))-rep(1, 10000)%*%t(bt2[i,])
  sumSquaredErr2=sumSquaredErr2+sum((M1-M_c)^2)
}
sumSquaredErr2/(10000*400*2000)#MSE for M hat is 0.0326721



#Evaluate the oracle variance
temp2=matrix(0, nrow = N2, ncol = J2)
for (i in 1:N2){
  for (j in 1:J2){
    temp2[i,j] = theta2[i]-beta2[j]
  }
}
p_mat2 = 1/(1+exp(-temp2))

vaTrue2=p_mat2*(1-p_mat2)
vaTrue2[is.na(da2)==T]=0
vaTrue2_i=apply(vaTrue2, 1, sum)
vaTrue2_j=apply(vaTrue2, 2, sum)
vaTrue2_=sum(vaTrue2)

#Compute estimated asymptotic variance using sample estimates
a_est2= apply(par_sample2$theta_sample, 2, mean)
b_est2=apply(par_sample2$beta_sample, 2, mean)
MEst2= a_est2 %*% t(rep(1, J2)) - rep(1, N2) %*% t(b_est2)
va2=matrix(0, nrow=N2, ncol=J2)
for (i in 1:N2){
  for (j in 1:J2){
    va2[i,j]=exp(-MEst2[i,j])/((1+exp(-MEst2[i,j]))^2)
  }
}
va2[is.na(da2)==T]=0
va2_i=apply(va2, 1, sum)
va2_j=apply(va2, 2, sum)
va2_=sum(va2)

#Theoretical variance
var2M11=1/va2_i[1]+1/va2_j[1]
var2Theta1=1/va2_i[1]
var2Beta1=1/va2_j[1]

M11Est2=par_sample2$theta_sample[,1]-par_sample2$beta_sample[,1]#sample estimates of M11



#Plot density curves
pdf("density_plots_case2_corrected.pdf", width=12, height=4)
par(mfrow=c(1,3),cex=0.7, mai=c(0.6,0.6,0.6,0.6))
hist(M11Est2[[1]], breaks=seq(-0.2, 1.2, 0.12),xlim=c(-0.1, 1.2),ylim=c(0,4),freq=F, col='cadetblue',axes=F, xlab='', ylab='Density', main='(e)') 
lines(seq(-0.2, 1.2,0.01), generate_norm_density(sigma=var2M11^0.5,mu=mean(M11Est2[[1]]),-0.2, 1.2,0.01))
axis(2, at=seq(0,3,0.5), labels =T)
axis(1, at=seq(-0.1, 1.2,0.1), labels = T, cex.axis=1, las=0)


hist(par_sample2$theta_sample[[1]], breaks=seq(-2, -0.3, 0.12),xlim=c(-2, -0.3),col='cadetblue',ylim=c(0,3),freq=F,axes=F, xlab='', ylab='Density', main='(f)') 
lines(seq(-2, -0.3,0.01), generate_norm_density(sigma=var2Theta1^0.5,mu=mean(par_sample2$theta_sample[[1]]),-2, -0.3,0.01))
axis(2, at=seq(0,3,0.5), labels =T)
axis(1, at=seq(-2, -0.3,0.5), labels = T, cex.axis=1, las=0);

hist(par_sample2$beta_sample[[1]], breaks=seq(-2.1, -1.69, 0.025),xlim=c(-2, -1.69),col='cadetblue',ylim=c(0,12),freq=F,axes = F, xlab='', ylab='Density', main='(g)') 
lines(seq(-2, -1.69,0.01), generate_norm_density(sigma=var2Beta1^0.5,mu=mean(par_sample2$beta_sample[[1]]),-2, -1.69,0.01))
axis(2, at=seq(0, 12, 3), labels =T)
axis(1, at=seq(-2, -1.5,0.1), labels = T, cex.axis=1, las=0);
dev.off()








#### Compare estimated sample variance to oracle variance and derived variance approximate in the paper for both settings ###
#oracle variance using true parameters
N1=5000
J1=200
set.seed(1)
theInd1=sample(1:N1, 100)
betInd1=sample(1:J1, 100)
MI1=sample(1:N1, 10)
MJ1=sample(1:J1, 10)
theInd2=sample(1:N2, 100)
betInd2=sample(1:J2, 100)
MI2=sample(1:N2, 10)
MJ2=sample(1:J2, 10)
optThetaTrue=1/vaTrue_i
optBetaTrue=1/vaTrue_j
optWTrue=optThetaTrue%*%t(rep(1, J1))+rep(1,N1)%*%t(optBetaTrue)
optThetaTrue=optThetaTrue[theInd1]
optBetaTrue=optBetaTrue[betInd1]
optWTrue=optWTrue[MI1, MJ1]
optThetaTrue2=1/vaTrue2_i
optBetaTrue2=1/vaTrue2_j
optWTrue2=optThetaTrue2%*%t(rep(1, J2))+rep(1,N2)%*%t(optBetaTrue2)
optThetaTrue2=optThetaTrue2[theInd2]
optBetaTrue2=optBetaTrue2[betInd2]
optWTrue2=optWTrue2[MI2, MJ2]




###Oracle variance using estimated parameters###

optThetaEst=NULL
optBetaEst=NULL
optWEst=NULL
thetaSample=as.matrix(par_sample$theta_sample)
betaSample=as.matrix(par_sample$beta_sample)
for (i in 1:2000){
  MEst=thetaSample[i,]%*%t(rep(1, J1))-rep(1, N1)%*%t(betaSample[i,])
  v=exp(-MEst)/((1+exp(-MEst))^2)
  v[is.na(da)==T]=0
  v_i=apply(v, 1, sum)
  v_j=apply(v, 2, sum)
  vM=(1/v_i)%*%t(rep(1, J1))+rep(1,N1)%*%t(1/v_j)
  vt=1/v_i
  vb=1/v_j
  optThetaEst=rbind(optThetaEst, vt[theInd1])
  optBetaEst=rbind(optBetaEst, vb[betInd1])
  optWEst=rbind(optWEst, unmatrix(vM[MI1, MJ1], byrow = T))
}
optThetaEst=apply(optThetaEst, 2, mean) 
optBetaEst=apply(optBetaEst, 2, mean)
optWEst=apply(optWEst, 2, mean)


optThetaEst2=NULL
optBetaEst2=NULL
optWEst2=NULL
thetaSample2=as.matrix(par_sample2$theta_sample)
betaSample2=as.matrix(par_sample2$beta_sample)
for (i in 1:2000){
  MEst2=thetaSample2[i,]%*%t(rep(1, J2))-rep(1, N2)%*%t(betaSample2[i,])
  v2=exp(-MEst2)/((1+exp(-MEst2))^2)
  v2[is.na(da2)==T]=0
  v2_i=apply(v2, 1, sum)
  v2_j=apply(v2, 2, sum)
  vM2=(1/v2_i)%*%t(rep(1, J2))+rep(1,N2)%*%t(1/v2_j)
  vt2=1/v2_i
  vb2=1/v2_j
  optThetaEst2=rbind(optThetaEst2, vt2[theInd2])
  optBetaEst2=rbind(optBetaEst2, vb2[betInd2])
  optWEst2=rbind(optWEst2, unmatrix(vM2[MI2, MJ2], byrow = T))
}

optThetaEst2=apply(optThetaEst2, 2, mean) 
optBetaEst2=apply(optBetaEst2, 2, mean)
optWEst2=apply(optWEst2, 2, mean)



#Sample variances
samThetaVar=apply(par_sample$theta_sample, 2, var)
samThetaVar=samThetaVar[theInd1]
samBetaVar=apply(par_sample$beta_sample, 2, var)
samBetaVar=samBetaVar[betInd1]
samWVar=NULL
for (i in 1:5000){
  c=par_sample$theta_sample[,i]%*%t(rep(1, J))-par_sample$beta_sample
  v=apply(c, 2, var)
  samWVar=rbind(samWVar, v)
}
samWVar=samWVar[MI1, MJ1]

samThetaVar2=apply(par_sample2$theta_sample, 2, var)
samThetaVar2=samThetaVar2[theInd2]
samBetaVar2=apply(par_sample2$beta_sample, 2, var)
samBetaVar2=samBetaVar2[betInd2]
samWVar2=NULL
for (i in 1:N2){
  c=par_sample2$theta_sample[[i]]%*%t(rep(1, J))-par_sample2$beta_sample
  v=apply(c, 2, var)
  samWVar2=rbind(samWVar2, v)
}
samWVar2=samWVar2[MI2, MJ2]


#plots---randomly sample 100
thetaMat=cbind(optThetaTrue, optThetaEst, samThetaVar)
colnames(thetaMat)=c('true', 'est', 'sam')
betaMat=cbind(optBetaTrue, optBetaEst, samBetaVar)
colnames(betaMat)=c('true', 'est', 'sam')
MMat=cbind(unmatrix(optWTrue, byrow = T), optWEst, unmatrix(samWVar, byrow = T))
colnames(MMat)=c('true', 'est', 'sam')
thetaMat2=cbind(optThetaTrue2, optThetaEst2, samThetaVar2)
colnames(thetaMat2)=c('true', 'est', 'sam')
betaMat2=cbind(optBetaTrue2, optBetaEst2, samBetaVar2)
colnames(betaMat2)=c('true', 'est', 'sam')
MMat2=cbind(unmatrix(optWTrue2, byrow = T), optWEst2, unmatrix(samWVar2, byrow = T))
colnames(MMat2)=c('true', 'est', 'sam')


thetaMat=rbind(thetaMat, thetaMat2)
betaMat=rbind(betaMat, betaMat2)
MMat=rbind(MMat, MMat2)

thetaMatall=cbind(thetaMat, c(rep(2, 100), rep(4, 100)), c(rep(19, 100),rep(24, 100)))
betaMatall=cbind(betaMat, c(rep(2, 100), rep(4, 100)),c(rep(19, 100),rep(24, 100)))
MMatall=cbind(MMat, c(rep(2, 100), rep(4, 100)),c(rep(19, 100),rep(24, 100)))

pdf("varPair_sampleVsEst.pdf", width=12, height=4)
par(mfrow=c(1,3),cex=0.7, mai=c(0.6,0.6,0.6,0.6))
plot(MMatall[,2], MMatall[,3], col=MMatall[,4], xlab='', ylab='', 
     main='(a)\n ',xlim=c(0, 0.1), axes=F,ylim = c(0, 0.1), pch=MMatall[,5])
abline(a=0, b=1)
axis(2, at=seq(0,0.1,0.02), labels =T)
axis(1, at=seq(0,0.1,0.02), labels = T, cex.axis=1, las=0)
legend(x=0, y=0.09, legend=c("Setting 1",'Setting 2'), col=c('red',"blue"),
       pch=c(19, 24), pt.cex=1, cex=1, text.font=2, y.intersp = 2, lwd=1, bty='n')

plot(thetaMatall[,2], thetaMatall[,3], col=thetaMatall[,4], xlab='', ylab='', 
     main='(b)',xlim=c(0, 0.1), axes=F,ylim = c(0, 0.1), pch=thetaMatall[,5])
abline(a=0, b=1)
axis(2, at=seq(0,0.1,0.02), labels =T)
axis(1, at=seq(0,0.1,0.02), labels = T, cex.axis=1, las=0)
legend(x=0, y=0.09, legend=c("Setting 1",'Setting 2'), col=c('red',"blue"),
       pch=c(19, 24), pt.cex=1, cex=1, text.font=2, y.intersp = 2, lwd=1, bty='n')

plot(betaMatall[,2], betaMatall[,3], col=betaMatall[,4], xlab='', ylab='', 
     main='(c)',xlim=c(0, 0.005), axes=F,ylim = c(0, 0.005), pch=betaMatall[,5])
abline(a=0, b=1)
axis(2, at=seq(0,0.005,0.001), labels =T)
axis(1, at=seq(0,0.005,0.001), labels = T, cex.axis=1, las=0)
legend(x=0, y=0.0045, legend=c("Setting 1",'Setting 2'), col=c('red',"blue"),
       pch=c(19, 24), pt.cex=1, cex=1, text.font=2, y.intersp = 2, lwd=1, bty='n')

dev.off()

pdf("varPair_estVstrue.pdf", width=12, height=4)
par(mfrow=c(1,3),cex=0.7, mai=c(0.6,0.6,0.6,0.6))
plot(MMatall[,1], MMatall[,2], col=MMatall[,4], xlab='', ylab='', 
     main='(d)',xlim=c(0, 0.1), axes=F,ylim = c(0, 0.1), pch=MMatall[,5])
abline(a=0, b=1)
axis(2, at=seq(0,0.1,0.02), labels =T)
axis(1, at=seq(0,0.1,0.02), labels = T, cex.axis=1, las=0)
legend(x=0, y=0.09, legend=c("Setting 1",'Setting 2'), col=c('red',"blue"),
       pch=c(19, 24), pt.cex=1, cex=1, text.font=2, y.intersp = 2, lwd=1, bty='n')

plot(thetaMatall[,1], thetaMatall[,2], col=thetaMatall[,4], xlab='', ylab='', 
     main='(e)',xlim=c(0, 0.1), axes=F,ylim = c(0, 0.1), pch=thetaMatall[,5])
abline(a=0, b=1)
axis(2, at=seq(0,0.1,0.02), labels =T)
axis(1, at=seq(0,0.1,0.02), labels = T, cex.axis=1, las=0)
legend(x=0, y=0.09, legend=c("Setting 1",'Setting 2'), col=c('red',"blue"),
       pch=c(19, 24), pt.cex=1, cex=1, text.font=2, y.intersp = 2, lwd=1, bty='n')

plot(betaMatall[,1], betaMatall[,2], col=betaMatall[,4], xlab='', ylab='', 
     main='(f)',xlim=c(0, 0.005), axes=F,ylim = c(0, 0.005), pch=betaMatall[,5])
abline(a=0, b=1)
axis(2, at=seq(0,0.005,0.001), labels =T)
axis(1, at=seq(0,0.005,0.001), labels = T, cex.axis=1, las=0)
legend(x=0, y=0.0045, legend=c("Setting 1",'Setting 2'), col=c('red',"blue"),
       pch=c(19, 24), pt.cex=1, cex=1, text.font=2, y.intersp = 2, lwd=1, bty='n')
dev.off()






##### To evaluate the coverage rate of the empirical 95% confidence intervals for theta, beta and M #####

#Setting 1
thetaCILower=NULL
thetaCIUpper=NULL
betaCILower=NULL
betaCIUpper=NULL
q95=qnorm(0.975)
th1=as.matrix(par_sample$theta_sample)
bt1=as.matrix(par_sample$beta_sample)
MTrue=theta%*%t(rep(1, J1))-rep(1, N1)%*%t(beta)
MProp=rep(0, N1*J1)
for (i in 1:2000){
  M=th1[i,]%*%t(rep(1, 200))-rep(1, 5000)%*%t(bt1[i, ])
  v=(1/(1+exp(-M)))*(1-1/(1+exp(-M)))
  v[is.na(da)==T]=0
  v_i=1/apply(v, 1, sum)
  v_j=1/apply(v, 2, sum)
  thetaCILower=rbind(thetaCILower, th1[i,]-q95*v_i^0.5)
  thetaCIUpper=rbind(thetaCIUpper, th1[i,]+q95*v_i^0.5)
  betaCILower=rbind(betaCILower, bt1[i,]-q95*v_j^0.5)
  betaCIUpper=rbind(betaCIUpper, bt1[i,]+q95*v_j^0.5)
  Mv=(v_i)%*%t(rep(1, J1))+rep(1,N1)%*%t(v_j)
  ML=M-q95*(Mv^0.5)
  MU=M+q95*(Mv^0.5)
  M_c=matrix(1, N1, J1)
  M_c[MTrue<ML]=0
  M_c[MTrue>MU]=0
  MProp=MProp+unmatrix(M_c, byrow = T)
}


proportion=function(x){
  return(sum(x)/length(x))
}

thMat=rep(1, 2000)%*%t(theta)
btMat=rep(1,2000)%*%t(beta)
thetaCov=matrix(1, 2000, 5000)
thetaCov[thMat<thetaCILower]=0
thetaCov[thMat>thetaCIUpper]=0
betaCov=matrix(1, 2000, 200)
betaCov[btMat<betaCILower]=0
betaCov[btMat>betaCIUpper]=0
thetaCovprob=apply(thetaCov, 2, proportion) #mean=0.9516239, sd=0.006270236
betaCovprob=apply(betaCov, 2, proportion)#mean=0.9416075, sd=0.007601719
MCovprob=MProp/2000#0.9510564, sd=0.00497511



#Setting 2
thetaCILower2=NULL
thetaCIUpper2=NULL
betaCILower2=NULL
betaCIUpper2=NULL
q95=qnorm(0.975)
th2=as.matrix(par_sample2$theta_sample)
bt2=as.matrix(par_sample2$beta_sample)
MTrue2=theta2%*%t(rep(1, J2))-rep(1, N2)%*%t(beta2)
MProp2=rep(0, N2*J2)
for (i in 1:2000){
  M=th2[i,]%*%t(rep(1, 400))-rep(1, 10000)%*%t(bt2[i, ])
  v=(1/(1+exp(-M)))*(1-1/(1+exp(-M)))
  v[is.na(da2)==T]=0
  v_i=1/apply(v, 1, sum)
  v_j=1/apply(v, 2, sum)
  thetaCILower2=rbind(thetaCILower2, th2[i,]-q95*v_i^0.5)
  thetaCIUpper2=rbind(thetaCIUpper2, th2[i,]+q95*v_i^0.5)
  betaCILower2=rbind(betaCILower2, bt2[i,]-q95*v_j^0.5)
  betaCIUpper2=rbind(betaCIUpper2, bt2[i,]+q95*v_j^0.5)
  
  Mv=(v_i)%*%t(rep(1, J2))+rep(1,N2)%*%t(v_j)
  ML=M-q95*(Mv^0.5)
  MU=M+q95*(Mv^0.5)
  M_c=matrix(1, N2, J2)
  M_c[MTrue2<ML]=0
  M_c[MTrue2>MU]=0
  MProp2=MProp2+unmatrix(M_c, byrow = T)
}


thMat2=rep(1, 2000)%*%t(theta2)
btMat2=rep(1,2000)%*%t(beta2)
thetaCov2=matrix(1, 2000, 10000)
thetaCov2[thMat2<thetaCILower2]=0
thetaCov2[thMat2>thetaCIUpper2]=0
betaCov2=matrix(1, 2000, 400)
betaCov2[btMat2<betaCILower2]=0
betaCov2[btMat2>betaCIUpper2]=0
thetaCovprob2=apply(thetaCov2, 2, proportion)# mean=0.9507013, sd=0.005559579
betaCovprob2=apply(betaCov2, 2, proportion)# mean=0.9452112, sd=0.005991096
MCovprob2=MProp2/2000#0.9504628, #0.004869668




#Box-plots
pdf("coverage_prob.pdf", width=12, height=4)
par(mfrow=c(1,2),cex=0.7, mai=c(0.6,0.6,0.6,0.6))
boxplot(list(MCovprob, thetaCovprob, betaCovprob), axes=F, 
        xlab='', ylab='', main='(a)', ylim=c(0.9,1), col = 'cadetblue')
axis(1, at=seq(0,4,1), labels = c('','M', expression(theta),expression(beta), ''))
axis(2, at=seq(0.9,1,0.01), labels = T, cex.axis=1, las=0)

boxplot(list(MCovprob2, thetaCovprob2,betaCovprob2), axes=F, 
        xlab='', ylab='', main='(b)', ylim=c(0.9,1), col = 'cadetblue')
axis(1, at=seq(0,4,1), labels = c('','M', expression(theta),expression(beta), ''))
axis(2, at=seq(0.9,1,0.01), labels = T, cex.axis=1, las=0)
dev.off()





















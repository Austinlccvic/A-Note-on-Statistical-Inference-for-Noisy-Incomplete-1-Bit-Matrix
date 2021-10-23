library(plotrix)


#Create data set
N=nrow(ADMneatX);
q=ncol(ADMneatX);
missing_mat=matrix(NA, nrow = N, ncol = 80);
mat1=cbind(ADMneatX[,c(41:q, 1:40)], missing_mat); #Re-align the data matrix so that anchor items from 81 to 120
mat2=cbind(missing_mat, ADMneatY);
data_mat=rbind(data_mat1, data_mat2);


#Heatplot
Z=matrix(0, 4000, 200);
Z[is.na(data_mat)==F]=1;
Z[is.na(data_mat)==T]=0;
Z=Z[1901:2100,];
colnames(Z)=NULL;
row.names(Z)=NULL;
write.table(Z, sep=",",  col.names=FALSE,row.names=FALSE,file = 'Z.csv');
heatmap(Z);

###Training
J = dim(data_mat)[2];
N = dim(data_mat)[1];
set.seed(1) #Ensure reproducibility
beta.vec = runif(J, -2, 2);
theta.vec = runif(N,-2, 2);
theta.vec=theta.vec-mean(theta.vec);
result = Rasch_missing(data_mat, beta.vec, theta.vec);
write.csv(result, file = 'real_result.csv');


#####Plot of confidence intervals constructed####
thetaHat=result$theta.vec;
betaHat=result$beta.vec;
MHat= thetaHat %*% t(rep(1, J)) - rep(1, N) %*% t(betaHat);
va=matrix(0, nrow=N, ncol=J); #Empirical version of var(Y)
for (i in 1:N){
  for (j in 1:J){
    va[i,j]=exp(-MHat[i,j])/((1+exp(-MHat[i,j]))^2);
  }
}
va[is.na(data_mat)==T]=0;
va_i=apply(va, 1, sum); #Vector of sigma_{i+}^2
va_j=apply(va, 2, sum); #Vector of sigma_{+j}^2
#va_=sum(va) #sigma_{++}^2

sigmaTheta=(1/va_i)^0.5; #Empirical version of sd of hat{theta}_i
sigmaBeta=(1/va_j)^0.5; #Empirical version of sd of hat{beta}_j

#####Plot 95% CI of theta's #####
qtl=qnorm(0.975); #0.975 quantile of N(0,1)
upper= thetaHat+qtl*sigmaTheta; #upper bound of CI
lower=thetaHat-qtl*sigmaTheta; #lower bound of CI
#For illustration purpose, we plot CI for subjects 1900 to 2100
mu=thetaHat[1976:2025];
indAsc=order(mu);
mu=mu[indAsc];
u=upper[1976:2025][indAsc];
l=lower[1976:2025][indAsc];

plotCI(x=1:50, y=mu, ui=u, li=l, axes=F, xlab="People index",ylab='Estimates of people parameter');
par(las=2);
axis(1, at=1:50, labels = as.character(indAsc[1:50]), cex.axis=0.6);
axis(2, at=seq(-3,3.5,0.5), labels =T);

#####Plot 95% CI of beta's #####
qtl=qnorm(0.975); #0.975 quantile of N(0,1)
upper= betaHat+qtl*sigmaBeta; #upper bound of CI
lower=betaHat-qtl*sigmaBeta; #lower bound of CI
#e plot CI for item 1 to 200
mu=betaHat[1:50];
indAsc=order(mu);
mu=mu[indAsc];
u=upper[1:50][indAsc];
l=lower[1:50][indAsc];

plotCI(x=1:50, y=mu, ui=u, li=l, axes=F, xlab="Item index",ylab='Estimates of item parameter');
par(las=2);
axis(1, at=1:50, labels = as.character(indAsc[1:50]), cex.axis=0.6);
axis(2, at=seq(-3,3,0.5), labels =T)






































































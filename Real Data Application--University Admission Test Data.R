library(plotrix)

#We use ADM data set. It contains binary responses from two forms of a college admission test.  
#Each form has 120 items and is answered by 2000 examinees.  
#There are 40 common items shared by the two test forms.  
#There is no missing data within each test.  
#N= 4000, J= 200, and 40% of the data entries are missing.

load(file = "ADMneatX.rda")

#Preprocessing
N=nrow(ADMneatX);
q=ncol(ADMneatX);
missing_mat=matrix(NA, nrow = N, ncol = 80);
mat1=cbind(ADMneatX[,c(41:q, 1:40)], missing_mat); #Re-align the data matrix so that anchor items from 81 to 120
mat2=cbind(missing_mat, ADMneatY);
colnames(mat1)=NULL;
colnames(mat2)=NULL;
data_mat=rbind(as.matrix(mat1), as.matrix(mat2));


#Heat-plot of the missing patterns
Z=matrix(0, 4000, 200);
Z[is.na(data_mat)==F]=1;
Z[is.na(data_mat)==T]=0;
Z=Z[1901:2100,];
colnames(Z)=NULL;
row.names(Z)=NULL;
write.table(Z, sep=",",  col.names=FALSE, row.names=FALSE, file = 'Z.csv');
heatmap(Z);

###Training
J = dim(data_mat)[2];
N = dim(data_mat)[1];
set.seed(1) 
beta.vec = runif(J, -2, 2);
theta.vec = runif(N,-2, 2);
theta.vec=theta.vec-mean(theta.vec);
result = Rasch_missing(data_mat, beta.vec, theta.vec);
#write.csv(result, file = 'real_result.csv');
#write.csv(result$beta.vec, file = 'trained_beta.csv')
#write.csv(result$theta.vec, file = 'trained_theta.csv')



#####Plot of confidence intervals constructed####
thetaHat=result$theta.vec;
betaHat=result$beta.vec;
N=4000;
J=200;
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


sigmaTheta=(1/va_i)^0.5; #Empirical sd of hat{theta}_i
sigmaBeta=(1/va_j)^0.5; #Empiricalsd of hat{beta}_j

#####Plot 95% CI of theta's #####
qtl=qnorm(0.975); #0.975 quantile of N(0,1)
upper= thetaHat+qtl*sigmaTheta; #upper bound of CI
lower=thetaHat-qtl*sigmaTheta; #lower bound of CI

#For illustration purpose, we plot CI for subjects 1900 to 2100
mu=thetaHat
set.seed(3)
group1Ind=sample(1:2000, 50)
group2Ind=sample(2001:4000, 50)
peopleMat=cbind(c(mu[group1Ind], mu[group2Ind]), c(upper[group1Ind], upper[group2Ind]), 
                c(lower[group1Ind], lower[group2Ind]), c(rep(2, 50), rep(4, 50)),
                c(rep(1,50), rep(2,50)), c(rep(19, 50), rep(24, 50)), c(rep(0, 50, 4, 50)))
colnames(peopleMat)=c('mu', 'upper', 'lower', 'color','lineType', 'pty', 'bg')
orderedIdx=order(peopleMat[,1])
peopleMat=peopleMat[orderedIdx, ]



#####Plot 95% CI of beta's #####
qtl=qnorm(0.975); #0.975 quantile of N(0,1)
upper= betaHat+qtl*sigmaBeta; #upper bound of CI
lower=betaHat-qtl*sigmaBeta; #lower bound of CI
set.seed(2)
group1Ind=sample(1:80, 40)
group2Ind=sample(121:200, 40)
anchorInd=sample(81:120, 20)
#plot CI for 100 items 
mu=betaHat
itemMat=cbind(c(mu[group1Ind],mu[anchorInd], mu[group2Ind]),c(upper[group1Ind],upper[anchorInd], upper[group2Ind]), 
              c(lower[group1Ind],lower[anchorInd], lower[group2Ind]),  c(rep(2, 40), rep(3, 20), rep(4, 40)),
              c(rep(1,40), rep(2,20), rep(3,40)), c(rep(19, 40),rep(1, 20),rep(4, 40)), 
              c(rep(0, 40), rep(1,20), rep(4, 40)))
colnames(peopleMat)=c('mu', 'upper', 'lower', 'color','lineType', 'pty', 'bg')
orderedIdx=order(itemMat[,1])
itemMat=itemMat[orderedIdx, ]



pdf("real_CI_corrected.pdf", width=12, height=4)
par(mfrow=c(1,2),cex=0.85, mai=c(0.6,0.6,0.6,0.6))
plotCI(x=1:100,y=peopleMat[,1],ui=peopleMat[,2],axes=F, li=peopleMat[,3], col = peopleMat[,4],
       xlab="",ylab='People parameters', main='(a)\n Examinees')#, slty=peopleMat[,5],pt.bg=peopleMat[,7], pch=peopleMat[,6],xlab="",ylab='People parameters');
axis(1, labels =F);
axis(2, at=seq(-3,3,0.5), labels =T)
legend(x=2, y=3, legend=c("Group 1", "Group 2"), col=c('red',"blue"),
       pch=c(1, 1), lty=c(1,1), pt.cex=1, cex=1, text.font=2,y.intersp = 1.5, lwd=1, bty='n')
#pch=c(19, 24), lty=c(1,2), pt.cex=1, cex=0.8, text.font=2,y.intersp = 0.5, lwd=1, bty='n')

plotCI(x=1:100,y=itemMat[,1],ui=itemMat[,2],axes=F, li=itemMat[,3], col = itemMat[,4], main='(b)\n Items',
       xlab="",ylab='Item parameters')#,slty=itemMat[,5],pt.bg=itemMat[,7], pch=itemMat[,6]);
axis(1, labels =F);
axis(2, at=seq(-3,3,0.5), labels =T)
legend(x=2, y=2, legend=c("Set 1",'Anchor set', "Set 2"), col=c('red','green',"blue"),
       pch=c(1, 1, 1), lty=c(1,1,1), pt.cex=1, cex=1, text.font=2,y.intersp = 1.5, lwd=1, bty='n')


dev.off()







































































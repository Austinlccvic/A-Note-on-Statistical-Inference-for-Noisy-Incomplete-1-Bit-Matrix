library(plotrix)
library(stringr) 
library(hash)
#Read in data
d1=read.dta(file='sen111kh.dta')
d2=read.dta(file='sen112kh.dta')
d3=read.dta(file='sen113kh.dta')

### Data pre-processing ###
idx1=which(d1$name=="OBAMA      ")
idx2=which(d2$name=="OBAMA      ")
idx3=which(d3$name=="OBAMA      ")
d1=d1[-idx1,]
d2=d2[-idx2,]
d3=d3[-idx3,]
dat1=cbind(d1$id, d1[,-c(1:9)])
colnames(dat1)[1]=c('id')
dat2=cbind(d2$id, d2[,-c(1:9)])
colnames(dat2)[1]=c('id')
dat3=cbind(d3$id,  d3[,-c(1:9)])
colnames(dat3)[1]=c('id')
dat=merge(dat1, dat2, by = "id", all = TRUE)
dat=merge(dat, dat3, by = "id", all = TRUE)
nam_r=c(d1$name, d2$name, d3$name)
id_r=c(d1$id, d2$id, d3$id)
state_r=c(d1$lstate, d2$lstate, d3$lstate)
party_r=c(d1$party, d2$party, d3$party)
nam=c()
state=c()
party=c()
for (i in dat$id){
  idx=which(id_r==i)[1]
  nam=c(nam, nam_r[idx])
  state=c(state, state_r[idx])
  if (party_r[idx]==100){
    party=c(party, 'Dem')
  }
  if (party_r[idx]==200){
    party=c(party, 'Rep')
  }
  if(party_r[idx]!=200 & party_r[idx]!=100){party=c(party, 'Ind')}
}
idx_Dem=which(party=='Dem')
idx_Rep=which(party=='Rep')
dat=dat[, -1] #Remove id column



#Change 6 to 0 and change 9 to NA
for (i in 1:nrow(dat)){
  dat[i,][which(dat[i,]==0)]=NA
  dat[i,][which(dat[i,]==2)]=1
  dat[i,][which(dat[i,]==3)]=1
  dat[i,][which(dat[i,]==4)]=0
  dat[i,][which(dat[i,]==5)]=0
  dat[i,][which(dat[i,]==6)]=0
  dat[i,][which(dat[i,]==7)]=NA
  dat[i,][which(dat[i,]==8)]=NA
  dat[i,][which(dat[i,]==9)]=NA
}


#Flip the response of the bills that have higher percentage of support within dem than rep.
count=0
for (j in 1:ncol(dat)){
  bill=dat[, j]
  dem_bill=bill[idx_Dem]
  rep_bill=bill[idx_Rep]
  dem_bill=dem_bill[is.na(dem_bill)==F]
  rep_bill=rep_bill[is.na(rep_bill)==F]
  if (sum(dem_bill)/length(dem_bill) > sum(rep_bill)/length(rep_bill)){
    bill[bill==0]=2
    bill[bill==1]=0
    bill[bill==2]=1
    dat[, j]=bill
    count=count+1
  }
}
count#1112

summation=function(x){
  sum(x, na.rm = T)
}


count_notna=function(x){
  sum(is.na(x)==FALSE)
}

s1=apply(dat, 1, summation)
c1=apply(dat, 1, count_notna)
idx_fullsupport=which(s1==0) #3  5 82 87

#Remove people with very small observations
idx_goodwin=which(nam=="GOODWIN    ")
idx_fullsupport=c(idx_fullsupport, idx_goodwin)
nam[idx_fullsupport]#"KENNEDY  ED" "BIDEN      " "CLINTON    " "SALAZAR    ", "GOODWIN    "
nam=nam[-idx_fullsupport]
state=state[-idx_fullsupport]
party=party[-idx_fullsupport]
which(s1==c1) #0
dat=dat[-idx_fullsupport,]
dim(dat)#139 1839

s2=apply(dat, 2, summation) 
c2=apply(dat, 2, count_notna)
length(which(s2==0)) #3
length(which(s2==c2)) #188
dat=dat[,-c(which(s2==0),  which(s2==c2))] #Remove columns with all 0, NA and all 1.
dim(dat)#139 1648
# write.csv(dat, file = 'senator_voting_flipped_dem1.csv')
# write.csv(nam, file = 'senator_names_flipped_dem1.csv')
# write.csv(state, file = 'senator_states_flipped_dem1.csv', )
# write.csv(party, file = 'senator_states_flipped_dem1.csv', )



#Plot of missingness pattern
J = dim(dat)[2];
N = dim(dat)[1];
Z=matrix(0, N, J);
Z[is.na(dat)==F]=1;
Z[is.na(dat)==T]=0;
idx_group1= which(apply(Z[,1:550], 1, sum)>100)
idx_other=c(1:139)[-idx_group1]
Z=Z[c(idx_group1, idx_other),]
idx_group3= which(apply(Z[1:118,1100:1648], 1, sum)<100)
idx_group3_other=c(c(1:118)[-c(idx_group3)], 111:139)
Z=Z[c(idx_group3, idx_group3_other),]
colnames(Z)=NULL;
row.names(Z)=NULL;
write.table(Z, sep=",",  col.names=FALSE,row.names=FALSE,file = 'Z.csv');


### Training the model on the data ###
set.seed(1) 
beta.vec = runif(J, -2, 2);
theta.vec = runif(N,-2, 2);
theta.vec=theta.vec-mean(theta.vec);
result = Rasch_missing(dat, beta.vec, theta.vec);


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
va[is.na(dat)==T]=0;
va_i=apply(va, 1, sum); #Vector of sigma_{i+}^2
va_j=apply(va, 2, sum); #Vector of sigma_{+j}^2


sigmaTheta=(1/va_i)^0.5; #Empirical version of sd of hat{theta}_i
sigmaBeta=(1/va_j)^0.5; #Empirical version of sd of hat{beta}_j




#####Plot 95% CI of theta's #####
qtl=qnorm(0.975); #0.975 quantile of N(0,1)
upper= thetaHat+qtl*sigmaTheta; #upper bound of CI
lower=thetaHat-qtl*sigmaTheta; #lower bound of CI
mu=thetaHat
co=rep(4, N)
peopleMat=cbind(mu, upper, lower, co, rep(1,N), rep(19, N))
colnames(peopleMat)=c('mu', 'upper', 'lower', 'color','lineType', 'pty')
orderedIdx=order(peopleMat[,1])
peopleMat=peopleMat[orderedIdx, ]
name_order=nam[orderedIdx]
state_order=state[orderedIdx]
party_order=party[orderedIdx]
sd=sigmaTheta[orderedIdx]


plotCI(x=1:N, y=peopleMat[,1],ui=peopleMat[,2],axes=F, li=peopleMat[,3], col = peopleMat[,4],
       xlab="",ylab='', main='Senator Conservativeness Score')#, slty=peopleMat[,5],pt.bg=peopleMat[,7], pch=peopleMat[,6],xlab="",ylab='People parameters');
axis(1, labels =F);
axis(2, at=seq(-6,8,1), labels =T)




#Top 10 most liberal senators
name_order[1:10]
state_order[1:10]
party_order[1:10]
round(sd[1:10], 3)
round(peopleMat[,1][1:10], 2)
# "GOODWIN    " "SCHATZ     " "BURRIS     " "HIRONO     "
#"BOOKER     " "BALDWIN    " "BROWN      " "UDALL      "
#"DURBIN     " "CARDIN     "

#Texas, South Carolina, Utah, Nebraska, Kentucky
#Wisconsin, Florida, Kansas, Arizona, Arkansas

#"WEST VI" "HAWAII " "ILLINOI" "HAWAII " "NEW JER" "WISCONS"
#"OHIO   " "NEW MEX" "ILLINOI" "MARYLAN"

#sd 3.268 0.469 0.297 0.383 0.572 0.352 0.168 0.165 0.164 0.163

#score -5.93 -4.70 -4.39 -4.12 -4.10 -3.86 -3.85 -3.80 -3.79 -3.78 
tab1=cbind(name_order[1:10], party_order[1:10], state_order[1:10], round(peopleMat[,1][1:10], 2),round(sd[1:10], 3))
row.names(tab1)=NULL
colnames(tab1)=c('Senator', 'Party', 'State', 'Conservativeness Score', 's.e.')
xtable(tab1, type = "latex", file = "Liberty ranking.tex")




#Top 10 most conservative senators
name_order[length(name_order): (length(name_order)-9)]
state_order[length(name_order):(length(name_order)-9)]
party_order[length(name_order):(length(name_order)-9)]
#"DEMINT     " "LEE        " "CRUZ       " "COBURN     "
#"PAUL       " "SCOTT      " "BUNNING    " "JOHNSON    "
#"RISCH      " "INHOFE     "

#"SOUTH C" "UTAH   " "TEXAS  " "OKLAHOM" "KENTUCK" "SOUTH C"
#"KENTUCK" "WISCONS" "IDAHO  " "OKLAHOM"


tab2=cbind(name_order[length(name_order): (length(name_order)-9)], 
           party_order[length(name_order): (length(name_order)-9)], 
           state_order[length(name_order): (length(name_order)-9)], 
           round(peopleMat[,1][length(name_order): (length(name_order)-9)], 2),
           round(sd[length(name_order): (length(name_order)-9)], 3))
row.names(tab2)=NULL
colnames(tab2)=c('Senator', 'Party', 'State', 'Conservativeness Score', 's.e.')
xtable(tab2, type = "latex", file = "Liberty ranking.tex")


#All the rankings
## hash-2.2.6 provided by Decision Patterns
h <- hash() 
h[['ALABAMA']]='AL'
h[['ALASKA ']]='AK'
h[['ARIZONA']]='AZ'
h[['ARKANSA']]=	'AR'
h[['CALIFOR']]=	'CA'
h[['COLORAD']]=	'CO'
h[['CONNECT']]	='CT'
h[['DELAWAR']]	='DE'
h[['FLORIDA']]	='FL'
h[['GEORGIA']]=	'GA'
h[["HAWAII "]]=	'HI'
h[["IDAHO  "]]	='ID'
h[["ILLINOI"]]	='IL'
h[["INDIANA"]]	='IN'
h[["IOWA   "]]	='IA'
h[["KANSAS "]]	='KS'
h[["KENTUCK"]]	='KY'
h[["LOUISIA"]]	='LA'
h[["MAINE  "]]	='ME'
h[["MARYLAN"]]	='MD'
h[["MASSACH"]]	='MA'
h[["MICHIGA"]]	='MI'
h[["MINNESO"]]	='MN'
h[["MISSISS"]]	='MS'
h[["MISSOUR"]]	='MO'
h[["MONTANA"]]	='MT'
h[["NEBRASK"]]	='NE'
h[["NEVADA "]]	='NV'
h[["NEW HAM"]]	='NH'
h[["NEW JER"]]	='NJ'
h[["NEW MEX"]]	='NM'
h[["NEW YOR"]]	='NY'
h[["NORTH C"]]	='NC'
h[["NORTH D"]]	='ND'
h[["OHIO   "]]	='OH'
h[["OKLAHOM"]]	='OK'
h[["OREGON "]]	='OR'
h[["PENNSYL"]]='PA'
h[["RHODE I"]]='RI'
h[["SOUTH C"]]='SC'
h[["SOUTH D"]]	='SD'
h[["TENNESS"]]	='TN'
h[["TEXAS  "]]	='TX'
h[["UTAH   "]]='UT'
h[["VERMONT"]]	='VT'
h[["VIRGINI"]]='VA'
h[["WASHING"]]='WA'
h[["WEST VI"]]	='WV'
h[["WISCONS"]]	='WI'
h[["WYOMING"]]	='WY'



sta=c()
for (s in state){
  sta=c(sta, h[[s]])
}

#A table of senators' ranking
tab=cbind(str_to_title(nam), sta, party, round(mu,2), round(sigmaTheta,3))
decreasing_idx=order(mu, decreasing = T)
tab=tab[decreasing_idx, ]
tab=cbind(1:139, tab)
tab1_odd=tab[seq(1, 62, 2), ]
tab1_even=tab[seq(2, 62, 2), ]
tab1_combined=cbind(tab1_odd, tab1_even)
tab2_odd=tab[seq(63, 139, 2),]
tab2_even=tab[seq(64, 139, 2),]
tab2_even=rbind(tab2_even, c(0,0,0,0,0,0))
tab2_combined=cbind(tab2_odd, tab2_even)
row.names(tab1_combined)=NULL
row.names(tab2_combined)=NULL
#colnames(tab1_combined)=c('Senator', 'State', 'Party', 'Score', 's.e.')
xtable(tab1_combined, type = "latex", file = "ranking.tex")
xtable(tab2_combined, type = "latex", file = "ranking.tex")
row.names(tab)=NULL
colnames(tab)=c('Senator', 'State', 'Party', 'Conservativeness Score', 's.e.')
library(xtable)
xtable(tab, type = "latex", file = "ranking.tex")


































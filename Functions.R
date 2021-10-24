### Generate missing pattern used in the simulation studies ###
generate_missingness=function(N, J){ #split into 5 groups, N is divisible by 5 and J is divisible by 2
  Z1=matrix(NA, N, J)
  Nmin=N/5
  Jmin=J/2
  Z1[1:Nmin, 1:Jmin]=1
  Z1[(Nmin+1):(2*Nmin), (Jmin/2+1):(1.5*Jmin)]=1
  Z1[(2*Nmin+1):(3*Nmin), (Jmin+1):(2*Jmin)]=1
  Z1[(3*Nmin+1):(4*Nmin), c(1:(Jmin/2), (Jmin+1):(1.5*Jmin))]=1
  Z1[(4*Nmin+1):(5*Nmin), c((Jmin/2+1):Jmin, (1.5*Jmin+1):(2*Jmin))]=1
  return(Z1)
}



### Training ###
### Alternating gradient descent to train Rasch model parameters under missing data setting ###

#Arguments:
#data: data matrix with missing entries; 
#beta.vec: initial values for beta vector;
#theta.vec: initial values for theta vector; 
#tau is learning rate, default 0.3; 
#tol: tolerance level on the convergence of the log-likelihood.

#Output: a list consisting of theta estimates and beta estimates. 

Rasch_missing = function(data, beta.vec, theta.vec, tau = 0.3, tol = 0.001){
  N=length(theta.vec);
  J=length(beta.vec);
  temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec); #M matrix at initial values
  JML0 = -Inf; #Initialize step 0 log likelihood as negative infinity
  JML = sum(data * temp - log(1+exp(temp)), na.rm=T); #Initialize step 1 log likelihood at initial values of beta and theta
  
  while(JML - JML0 > tol){ #tolerance on the consecutive change in log likelihood
    JML0 = JML;
    prob = 1/(1+exp(-temp));
    theta.vec = theta.vec + tau * rowMeans(data-prob, na.rm=T); #update for theta estimates
    theta.vec=theta.vec-mean(theta.vec); #identifiability constraint
    temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec);
    prob = 1/(1+exp(-temp));
    beta.vec = beta.vec + tau * colMeans(-data + prob, na.rm=T); #update for beta estimates
    beta.vec=beta.vec-mean(theta.vec)
    temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec);
    JML = sum(data * temp - log(1+exp(temp)), na.rm = T);
    print(JML);
  }  
  list(beta.vec = beta.vec, theta.vec = theta.vec);
}



### Results post-processing ###

#Result is a list of trained results obtained from parallel computing package doParallel
#This function converts a list of results into a list of two matrices, consisting of trained theta values and beta values
unlist=function(result){
  the_v=result[[1]][[2]]
  beta_v=result[[1]][[1]]
  for (i in 2:length(result)){
    the_v=rbind(the_v, result[[i]][[2]])
    beta_v=rbind(beta_v, result[[i]][[1]])
  }
  return(list(the=the_v, beta=beta_v))
}


#Returns the oracle variance of theta_i
#Arguments:
#p_mat: a matrix with (i,j)'th entry being P(Y_{ij}=1) under true theta and beta values;
#i: index of individual i;
#da: missing pattern matrix.
oracle_var_theta=function(p_mat, i, da){
  va=p_mat*(1-p_mat)
  va[is.na(da)==T]=0
  oracle=1/sum(va[i,])
  return(oracle)
}

#Returns the oracle variance of beta_j for no additional covariates case
#Arguments:
#p_mat: a matrix with (i,j)'th entry being P(Y_{ij}=1) under true theta and beta values;
#j: index of item j;
#da: missing pattern matrix.
oracle_var_beta=function(p_mat, j, da){
  va=p_mat*(1-p_mat)
  va[is.na(da)==T]=0
  oracle=1/sum(va[,j])
  return(oracle)
}


#Compute theoretical normal density curve with mean mu and sd sigma, within range [l,u] with jumpsize of spc
generate_norm_density=function(sigma, mu, l, u, spc){
  norm_d=c()
  for (x in seq(l, u, spc)){
    d=dnorm(x, mean = mu, sd = sigma, log = FALSE)
    norm_d=c(norm_d, d)
  }
  return(norm_d)
}

###Alternating gradient descent to train Rasch model parameters under missing data setting###

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

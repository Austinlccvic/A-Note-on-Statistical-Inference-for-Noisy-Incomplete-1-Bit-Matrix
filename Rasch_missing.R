###Alternating gradient descent to train Rasch model parameters under missing data setting
#Input: data missing data matrix; beta.vec: initial values for beta vector;
#theta.vec: initial values for theta vector; tau is learning rate; tol: tolerance level on convergence of likelihood
#Output: a list of theta estimates and beta estimates 
Rasch_missing = function(data, beta.vec, theta.vec, tau = 0.3, tol = 0.001){
  N=length(theta.vec);
  J=length(beta.vec);
  temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec); #underlying M matrix
  JML0 = -1e10; #Initial value of joint likelihood
  JML = sum(data * temp - log(1+exp(temp)), na.rm=T); #joint likelihood
  
  while(JML - JML0 > tol){ #tolerance level on the consecutive change in joint likelihood
    JML0 = JML;
    prob = 1/(1+exp(-temp));
    theta.vec = theta.vec + tau * rowMeans(data-prob, na.rm=T); #update for theta estimates
    theta.vec=theta.vec-mean(theta.vec);
    temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec);
    prob = 1/(1+exp(-temp));
    beta.vec = beta.vec + tau * colMeans(-data + prob, na.rm=T);
    beta.vec=beta.vec-mean(theta.vec)
    temp = theta.vec %*% t(rep(1, J)) - rep(1, N) %*% t(beta.vec);
    JML = sum(data * temp - log(1+exp(temp)), na.rm = T);
    print(JML);
  }  
  list(beta.vec = beta.vec, theta.vec = theta.vec);
}
# component-wise EM as in M. Figueirdo, 'Unsupervised learning of finite mixture models' (2002)

require(mvtnorm)

# INPUTS: 
#   data - input data to fit GMM to 
#   kmin - minimum number of clusters to consider (integer > 0) 
#   kmax - maximum number of clusters to consider (integer > 0)
#   epsilon - stopping criterion
#   max_iter - maximum number of iterations to do
# OUTPUTS: 
#   k - "best" number of components for GMM
#   pi - mixture weights
#   Mu - mean vector/matrix for GMM
#   C - covariance matrices for GMM
#   L - likelihood of the data at the output mixture model
#   n_iter - number of iterations done before stopping criterion met

fit_cEM <- function(data, kmin, kmax, epsilon=1e-10, max_iter=2500){
  # set up 
  n=nrow(data)
  d=ncol(data)
  N=d + d*(d+1)/2
  
  t=0
  knz = kmax
  Lmin = Inf
  
  # initialise GMM parameters
  ## random points for mean vector/matrix 
  samp = sample(1:n, size=kmax, replace=FALSE)
  Mu = as.matrix(data[samp, , drop=FALSE])
  
  ## proportional to identity for covariances
  C = array(NA, dim=c(d,d,kmax))
  sigma = 1/10 * 1/d * sum(diag(cov(data)))
  for(m in 1:kmax){
    C[,,m] = sigma*diag(nrow=d)
  }
  
  ## mixing proportions equal 
  pi = rep(1/kmax, kmax)
  
  # initialise u
  u = matrix(NA, nrow=n, ncol=kmax)
  for(m in 1:kmax){
    u[,m] = dmvnorm(data, mean=Mu[m,], sigma=C[,,m])
  }
  
  # first L (L(theta(0), Y))
  term4 = 0
  for (i in 1:n){
    sum_m = 0
    for (m in 1:knz){
      sum_m = sum_m + sum(pi[m]*u[i,m])
    }
    term4 = term4 + log(sum_m)
  }
  Lprev = (N/2)*sum(log(n*pi[pi>0]/12)) + knz/2 * log(n/12) + (knz*N + knz)/2 + term4
  
  
  # component-wise EM
  while(knz >= kmin && t < max_iter) {
    #Lprev = -Inf
    
    repeat{
      t=t+1
      
      #w = matrix(0, nrow=n, ncol=knz)
      w = matrix(0, nrow=n, ncol=kmax)
      
    for(m in 1:kmax){
      # function for the internal for loop
      updates = for_m_to_k(m, pi, u, w, Mu, C, knz, data, kmax, d, N)
      pi = updates$pi
      Mu = updates$Mu
      C = updates$C
      w = updates$w
      u = updates$u
      #knz = updates$knz
      knz = sum(pi>0)
    }

    # calculate likelihood for new parameters
    term1 = (N/2)*sum(log(n*pi[pi>0]/12))
    term2 = knz/2 * log(n/12)
    term3 = (knz*N + knz)/2
    term4 = 0
    for (i in 1:n){
      sum_m = 0
      for (m in 1:knz){
        sum_m = sum_m + sum(pi[m]*u[i,m])
      }
      term4 = term4 + log(sum_m)
    }
    L = term1 + term2 + term3 + term4
    #if (abs(Lprev - L) < epsilon*abs(L)){
    if (abs(Lprev - L) < epsilon * abs(Lprev)){
      break
    }
    if (t >= max_iter){
      break
    }
    Lprev = L
    #if(t %% 10 ==0){print(t)}
    }
    
    # store best so far
    if(L <= Lmin){
      Lmin = L
      pi_best = pi
      Mu_best = Mu
      C_best = C
      k_best = sum(pi > 0)
      #k_best = knz
    }
    
    # remove component with lowest assignment prob
    #if (knz > kmin){
      # m_star = which.min(pi)
      # pi = pi[-m_star]
      # Mu = Mu[-m_star, , drop=FALSE]
      # C = C[,,-m_star, drop=FALSE]
      # pi = pi/sum(pi)
      # knz = knz-1
    m_star = which.min(pi[pi>0])
    pi[m_star] = 0
    pi = pi/sum(pi)
    knz = knz-1

    #} else{ break}
    
    #pi = pi/sum(pi)
    #knz = sum(pi > 0)
  }
  
  pi_nz = pi_best[which(pi_best>0)]
  Mu_nz = Mu_best[which(pi_best>0), ]
  C_nz = C_best[,,which(pi_best>0)]
  return(list(k=k_best, L=Lmin, pi=pi_nz, Mu=Mu_nz, C=C_nz, n_iter=t))
}



## function for for loop inside while loop 
for_m_to_k <- function(m, pi, u, w, Mu, C, knz, data, kmax, d, N){
  # update w
  denom = 0
  #for(j in 1:knz){
  for(j in 1:kmax){
    denom = denom + pi[j]*u[,j]
  }
  w[, m] = (pi[m]*u[,m])/denom
  
  # update pi
  denom_sum = 0
  for (j in 1:knz){
    denom_sum = denom_sum + max(c(0, sum(w[,j])-N/2))
  }
  pi[m] = max(c(0, sum(w[,m])-N/2)) / denom_sum
  pi = pi/sum(pi)
  
  # if pi > 0, update GMM parameters
  if(pi[m] >0){
    theta_hat = argmax_theta(data, w[,m], d)
    Mu[m,] = theta_hat$mu
    C[,,m] = theta_hat$Sigma
    
    # update u
    u[,m] = dmvnorm(data, Mu[m,], C[,,m])
  } else {
    knz = knz-1
  }
  
  return(list(pi=pi, Mu=Mu, C=C, u=u, w=w, knz=knz))
}


argmax_theta <- function(y, w, d) {
  mu_hat = apply(y, 2, function(x){sum(x*w)}) / sum(w)
  y_centred = sweep(y, 2, mu_hat, '-')
  Sigma_hat = t(y_centred) %*% (y_centred*w) / sum(w)
  Sigma_hat = Sigma_hat + 1e-6*diag(d)
  
  return(list(mu=mu_hat, Sigma=Sigma_hat))
}


## GMM_density function
gmm_density <- function(data, pi, mu, sigma, K){
  density <- 0
  if (K==1){
    density <- dmvnorm(data, mean=mu, sigma=sigma)
  } else {
    for (k in 1:K) {
      density <- density + pi[k]*dmvnorm(data, mean=mu[k,], sigma=sigma[,,k])
    }
  }
  
  return(density)
}


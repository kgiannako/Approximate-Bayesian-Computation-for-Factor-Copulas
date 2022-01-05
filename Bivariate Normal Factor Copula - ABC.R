#--------------------- Empirical density function ---------------------

EmpiricalDensity <- function(X){
  emp <- c()
  for(i in 1:length(X)){
    emp[i] <- sum(X<X[i])/length(X)
  }
  return(emp)
}


#--------------------- Quantile dependence ---------------------

QuantileDependence <- function(G1, G2, q){
  dep <- c()
  if (q<0 | q>1){
    return(cat("q must be between 0 and 1"))
  } else if (q<=0.5){
    dep <- sum(G1<=q & G2<=q)/(length(G1)*q)
  } else {
    dep <- sum(G1>q & G2>q)/(length(G1)*(1-q))
  }
  
  return(dep)
}

#--------------------- Rank correlation ---------------------

RankCorrelation <- function(G1,G2){
  
  spearrank <- 12 * sum(G1*G2)/length(G1) - 3
  return(spearrank)
  
}

#----------------------------------------------------------------


# Normal Factor Copula
#----------------------------------------------------------------

# Length of siumulated data
n_simulation <-10000 

# Real value of beta
beta <- c()
beta[1] <- 1/sqrt(2)
beta[2] <- 1/sqrt(2)

# Correlation as function of beta
(beta[1]^2) / (beta[1]^2 + 0.5)

#----------------Simulating Data from the factor copula------------------

set.seed(71)

# Factor Z~N(mu,1)
Z <- rnorm(n = n_simulation, mean = 0, sd=1)

# Epsilons epsilon~N(0,1/2)
epsilon<-matrix(nrow=n_simulation, ncol=2)
epsilon[,1] <- rnorm(n = n_simulation, mean=0, sd=1/sqrt(2))
epsilon[,2]<- rnorm(n = n_simulation, mean=0, sd=1/sqrt(2))

# Factor copula - observed data
Simulated_data <- matrix(nrow=n_simulation, ncol=2)
Simulated_data[,1] <- beta[1]*Z + epsilon[,1]
Simulated_data[,2] <- beta[2]*Z + epsilon[,2]

# factor copula scatterplot
plot(Simulated_data, main="Simulated data from Factor copula with correlation 0.5", pch=16,
     col="red", xlab = "X1", ylab="X2", xlim=c(-4,4), ylim=c(-4,4) )

#------------------------------------------------------------------------

#---------------------sample dependence measures-------------------------

m_T <- c()
# EDFs
e1 <- ecdf(Simulated_data[,1])
e2 <- ecdf(Simulated_data[,2])
G1 <- e1(Simulated_data[,1]) 
G2 <- e2(Simulated_data[,2])
#G1 <- EmpiricalDensity(Simulated_data[,1]) - ~100 slower than ecdf
#G2 <- EmpiricalDensity(Simulated_data[,2])

# compute summary statistics
m_T[1] <- RankCorrelation(G1, G2)
m_T[2] <- QuantileDependence(G1, G2, 0.05)
m_T[3] <- QuantileDependence(G1, G2, 0.10)
m_T[4] <- QuantileDependence(G1, G2, 0.90)
m_T[5] <- QuantileDependence(G1, G2, 0.95)


#-------------------------------------------------------------------------

#---------------------------- Bootstrap ----------------------------------
n_bootstrap_sim <- 10000 # number of bootstrao simulations

set.seed(842)
m_B <- matrix(nrow = 5, ncol = n_bootstrap_sim) # store summary statistcis of bootstraps simulations

pb = txtProgressBar(min = 0, max = n_bootstrap_sim, style=3, width=50, char="=") #progress bar
start <- Sys.time() #track time

for(i in 1:n_bootstrap_sim){
  #sample from the sim data
  Bootstrap_sample_n <- sample(x = 1:nrow(Simulated_data), size = n_simulation, replace = T) 
  Bootstrap_sample <- Simulated_data[Bootstrap_sample_n,] 
  
  # computing EDFs
  e1 <- ecdf(Bootstrap_sample[,1])
  e2 <- ecdf(Bootstrap_sample[,2])
  G2 <- e1(Bootstrap_sample[,1]) 
  G1 <- e2(Bootstrap_sample[,2])
  
  # Compute bootstrap summ stats
  m_B[1,i] <- RankCorrelation(G1, G2)
  m_B[2,i] <- QuantileDependence(G1, G2, 0.05)*100
  m_B[3,i] <- QuantileDependence(G1, G2, 0.10)
  m_B[4,i] <- QuantileDependence(G1, G2, 0.90)
  m_B[5,i] <- QuantileDependence(G1, G2, 0.95)
  
  setTxtProgressBar(pb,i) #progress bar
}
end <- Sys.time()
close(pb)

#time elapsed
end-start

# means of each sim summ statistic
mean_m_B <- cbind(mean(m_B[1,]),mean(m_B[2,]),mean(m_B[3,]),mean(m_B[4,]),mean(m_B[5,]))
mean_m <- as.vector(mean_m_B)

# Covariance matrix
Sigma <- matrix(0,nrow = 5, ncol = 5)
for(j in 1:n_bootstrap_sim){
  diff <- m_B[,j]-mean_m
  s <- diff %*% t(diff)
  Sigma <- Sigma  + s
}

library(matlib)
#scaled cov matrix
Sigma_T_B <- (n_simulation/n_bootstrap_sim) * Sigma

# Efficient weight matrix
W <- inv(Sigma_T_B)
Weights_100k <- W

#-------------------------------------------------------------------------


#------------------------- Simulations - ABC -----------------------------

#vector of simulated summary stats
m_S <- matrix(nrow=5, ncol=1)

# "distance" function for ABC - see Oh & Patton
Q_theta_I <- c() #for Identity matrix
Q_theta_W <- c() # for Efficient weight matrix

beta_sim <- c() # to store simulated beta values

ind <- 10000 # Number of simulations

M <- matrix(nrow=5, ncol = ind) # to store summary statistics of all simulations

pb = txtProgressBar(min = 0, max = ind, style=3, width=50, char="=") # setting up progress bar


for (k in 1:ind){
  
  #sample beta from U(0,1)
  set.seed(seed = NULL)
  beta_sim[k] <- runif(1,0,1)
  
  # Factor Z~N(mu,1)
  Z <- rnorm(n = n_simulation, mean = 0, sd=1)
  
  # Epsilons epsilon~N(0,1/2)
  epsilon<-matrix(nrow=n_simulation, ncol=2)
  epsilon[,1] <- rnorm(n = n_simulation, mean=0, sd=1/sqrt(2))
  epsilon[,2]<- rnorm(n = n_simulation, mean=0, sd=1/sqrt(2))
  
  # X = bZ+epsilon
  X <- matrix(nrow=n_simulation, ncol=2)
  X[,1] <- beta_sim[k]*Z + epsilon[,1]
  X[,2] <- beta_sim[k]*Z + epsilon[,2]
  
  # empirical CDF
  G <- matrix(nrow=n_simulation, ncol=2) 
  e1 <- ecdf(X[,1])
  e2 <- ecdf(X[,2])
  
  G[,1] <- e1(X[,1])  #G[,1] <- EmpiricalDensity(X[,1]) 
  G[,2] <- e2(X[,2])  #G[,2] <- EmpiricalDensity(X[,2])
  
  
  # Sim Summary statistics 
  m_S[1,1] <- RankCorrelation(G[,1], G[,2])
  m_S[2,1] <- QuantileDependence(G[,1], G[,2], 0.05)
  m_S[3,1] <- QuantileDependence(G[,1], G[,2], 0.10)
  m_S[4,1] <- QuantileDependence(G[,1], G[,2], 0.90)
  m_S[5,1] <- QuantileDependence(G[,1], G[,2], 0.95)
  
  
  #Average and store Summary statistics (for loop unnecessary)
  m <- c() 
  for (i in 1:5){
    m[i] <-mean(m_S[i,]) 
  }
  M[,k]<-m
  
  # Q(theta) for id matrix - Sum of squared differences
  Q_theta_I[k] <- sum((m_T-m)^2)
  # Q(theta) for weight matrix
  diff_m <- m_T-m
  Q_theta_W[k] <- t(diff_m) %*% Weights_100k %*% diff_m
  
  setTxtProgressBar(pb,k) #progress bar
  
  
}
close(pb) #progress bar


binded_I <- cbind(beta_sim,Q_theta_I) # for id matrix
binded_W <- cbind(beta_sim, Q_theta_W) #for weight matrix

#ordering according to Q(theta)
sorted_I <- binded_I[order(binded_I[,2]),] 
sorted_W <- binded_W[order(binded_W[,2]),] 

#column names
colnames(sorted_I) <- c("Beta", "Q(theta)")
colnames(sorted_W) <- c("Beta", "Q(theta)")


# Histograms - results
real_values <- c(beta[1]) #true parametr value
limits <- matrix(c(0.4,1), nrow=1, byrow = T) #histogram limits x-axis

#par(mfrow=c(2,4))
len <- n_simulation*0.005 #accepted fraction of sim values

# histogram for id
hist(sorted_I[1:len,i], main=paste(colnames(sorted_I)[i],":",100*len/ind,"% accepted") , breaks=15, xlab="w/o weights",
     col="lightgoldenrod", xlim = limits[i,])  
abline(v=real_values[i], col="red", lty=4, lwd=3) #true value
abline(v = mean(sorted_I[1:len,i]), col="black", lwd=3, lty = 4) #posterior mean

#histogram of weight
hist(sorted_W[1:len,i], main=paste(colnames(sorted_W)[i],":",100*len/ind,"% accepted") , breaks=15, xlab="w/ weights",
     col="cornflowerblue", xlim = limits[i,])  
abline(v=real_values[i], col="red", lty=4, lwd=3)
abline(v = mean(sorted_W[1:len,i]), col="black", lwd=3, lty = 4)

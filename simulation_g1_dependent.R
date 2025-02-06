Rcpp::sourceCpp('C:/Users/26085/Desktop/Annie/gaussian_dependent.cpp')
source('C:/Users/26085/Desktop/Annie/main-gaussian_dependent.R')
library(igraph)
library(lpSolve)
library(rTensor)
library(stats)
library(extraDistr)

# HOSVD

r=c(0.5,0.6,0.7)
rho=matrix(c(0.1,0.05,0.15,0.05,0.3,0.35,0.15,0.35,0.5),3,3)
mu0=matrix(c(0.03,0.02,0.04,0.02,0.07, 0.05, 0.04,0.05, 0.1),3,3)
mu1=2*mu0
sigma0=matrix(c(0.5,0.3,0.4,0.3,0.7,0.6,0.4,0.6,1),3,3)
sigma1=2*sigma0
set.seed(1)
n=120
run_time = 0
G=3
T=9
R=100
gpr=NULL
N=matrix(T,n,n)
for(i in 1:R){
  b=block_N(n,G)
  A = get_A(n,b,rho,r)
  fit0=get_X(A,n,b,T,mu0,mu1,sigma0,sigma1)
  Ecopy = fit0$E
  Y_tensor <- as.tensor(Ecopy)
  start_time <- Sys.time()
  svd <- hosvd(Y_tensor, ranks = c(G, 3, 3))  # Keeping G leading eigenvectors for the first mode
  U1 <- svd$U[[1]]  # Matrix corresponding to the first mode
  result <- kmeans(U1, centers = G)
  end_time <- Sys.time()
  run_time = run_time + end_time - start_time
  gp=result$cluster
  gpn=matrix(0,nrow=G,ncol=G)
  for(j in 1:G){
    for(k in 1:G){
      low=(k-1)*round(n/G)+1
      up=min(k*round(n/G),n)
      gpn[j,k]=sum(gp[low:up]==j)
    }
  }
  gpr[i] = lp.assign(gpn,direction = "max")$objval/n
  print(paste("Accuracy:", gpr[i]))
}
print(mean(gpr))
print(run_time/R)

# GDBM

r=c(0.5,0.6,0.7)
rho=matrix(c(0.1,0.05,0.15,0.05,0.3,0.35,0.15,0.35,0.5),3,3)
mu0=matrix(c(0.03,0.02,0.04,0.02,0.07, 0.05, 0.04,0.05, 0.1),3,3)
mu1=2*mu0
sigma0=matrix(c(0.5,0.3,0.4,0.3,0.7,0.6,0.4,0.6,1),3,3)
sigma1=2*sigma0
set.seed(2)
n=120
G=3
T=5
R=100
N=matrix(T,n,n)
b=block_N(n,G)
Res = SR_Gauss(N,T,n,G,rho,r,mu0,mu1,sigma0,sigma1,R)
print(Res$execution_time)

# Annie's method

r=c(0.5,0.6,0.7)
rho=matrix(c(0.1,0.05,0.15,0.05,0.3,0.35,0.15,0.35,0.5),3,3)
mu0=matrix(c(0.03,0.02,0.04,0.02,0.07, 0.05, 0.04,0.05, 0.1),3,3)
mu1=2*mu0
sigma0=matrix(c(0.5,0.3,0.4,0.3,0.7,0.6,0.4,0.6,1),3,3)
sigma1=2*sigma0
accuracy = 0
set.seed(2)
n=120
G=3
T=5
run_time = 0
R=100
N=matrix(T,n,n)
b=block_N(n,G)
for (cishu in 1:R) {
  A=get_A(n,b,rho,r)
  res=get_X(A,n,b,T,mu0,mu1,sigma0,sigma1)
  E = res$E
  for (i in 1:T) {
    X = E[,,i]
    diag(X) <- NA
    mean_value <- mean(X, na.rm = TRUE)
    X[X>=mean_value] = 1 
    X[X<mean_value] = 0
    diag(X) <- 0
    E[,,i] <- X
  }
  start_time <- Sys.time()
  result = SE(E,n,G)
  accuracy = accuracy+result$acc
  end_time <- Sys.time()
  run_time = run_time + end_time - start_time
  print(result$acc)
}
print(accuracy/R)
print(run_time/R)
Rcpp::sourceCpp('C:/Users/26085/Desktop/Annie/bernoulli_dependent.cpp')
source('C:/Users/26085/Desktop/Annie/main-bernoulli_dependent.R')
library(igraph)
library(lpSolve)
library(rTensor)
library(stats)

# HOSVD

set.seed(123)
alpha=matrix(c(0.2,0.3,0.5,0.3,0.5,0.6,0.5,0.6,0.7),3,3)
rho=matrix(c(0.3,0.4,0.5,0.4,0.5,0.6,0.5,0.6,0.7),3,3)
beta=matrix(c(0.01,0.01,0.02,0.01,0.02,0.04,0.02,0.04,0.05),3,3)
r=c(0.5,0.6,0.7)
options (warn = -1)
n=120
gpr=NULL
G=3
T=5
run_time = 0
R =100
N=matrix(T,n,n)
accuracy = 0
for(i in 1:R){
  b=block_N(n,G)
  fit0=get_AE(T,n,b,rho,r,alpha,beta)
  Ecopy = fit0$E_copy
  Y_tensor <- as.tensor(Ecopy)
  start_time <- Sys.time()
  svd <- hosvd(Y_tensor, ranks = c(G, 3, 3))
  U1 <- svd$U[[1]]
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
  gpr[i]=lp.assign(gpn,direction = "max")$objval/n
  print(paste("Accuracy:", gpr[i]))
}
print(mean(gpr))
print(run_time/R)

# GDBM

set.seed(123)
alpha=matrix(c(0.2,0.3,0.5,0.3,0.5,0.6,0.5,0.6,0.7),3,3)
rho=matrix(c(0.3,0.4,0.5,0.4,0.5,0.6,0.5,0.6,0.7),3,3)
beta=matrix(c(0.01,0.01,0.02,0.01,0.02,0.04,0.02,0.04,0.05),3,3)
r=c(0.5,0.6,0.7)
options (warn = -1)
n=100
G=3
T=5
R =100
N=matrix(T,n,n)
b=block_N(n,G)
Res = SR(T,N,n,G,r,rho,alpha,beta,R)
print(Res$execution_time)

# Annie's method

set.seed(123)
r=c(0.5,0.6,0.7)
alpha=matrix(c(0.2,0.3,0.5,0.3,0.5,0.6,0.5,0.6,0.7),3,3)
rho=matrix(c(0.3,0.4,0.5,0.4,0.5,0.6,0.5,0.6,0.7),3,3)
beta=matrix(c(0.01,0.01,0.02,0.01,0.02,0.04,0.02,0.04,0.05),3,3)
accuracy = 0
n=120
G=3
T=5
R=100
run_time = 0
N=matrix(T,n,n)
b=block_N(n,G)
for (cishu in 1:R) {
  fit0=get_AE(T,n,b,rho,r,alpha,beta)
  E=fit0$E_copy
  start_time <- Sys.time()
  result = SE(E,n,G)
  end_time <- Sys.time()
  run_time = run_time + end_time - start_time
  accuracy = accuracy+result$acc
  print(result$acc)
}
print(accuracy/R)
print(run_time/R)



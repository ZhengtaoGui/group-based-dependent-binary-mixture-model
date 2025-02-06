
block_N <- function(n,G){  #grouping
  b=NULL
  for(i in 1:G)
    if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
  else
  {b[((i-1)*round(n/G)+1):n]=i}
  return(b)
}

#generate A using SBM G=2
get_A <- function(n,b,rho,r){
  A=matrix(0, nrow=n, ncol=n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      index1=b[i]
      index2=b[j]
      if(index1==index2){
        p=r[index1]+(1-r[index1])*rho[index1,index2]
        x=(r[index1]*(1-rho[index1,index2]))/(r[index1]+(1-r[index1])*rho[index1,index2])
        B_ij=rbinom(1,1,p)
        C=runif(1,min=0,max=1)
        if(B_ij*C>x){
          A[i,j]=1
        }
        else{
          A[i,j]=0
        }
      }
      else{
        A[i,j]=rbinom(1,1,rho[index1,index2])
      }
      A[j,i]=A[i,j]
    }
  }
  return(A)
}
# b=block_N(n,G)
# A=get_A(n,b,rho,r)

# rnorm(T,mu1[index1,index2],sigma1[index1,index2])

get_X<-function(A,n,b,T,mu0,mu1,sigma0,sigma1){
  m=n*(n-1)/2
  X <- array(0, c(T, 4, m))
  E <- array(0,c(n,n,T))
  k=0
  df = 5  # 自由度为5
  for(i in 1:n){
    if((i+1)<=n){
      for(j in (i+1):n){
        k=k+1
        X[,1,k]=i
        X[,2,k]=j
        index1=b[i]
        index2=b[j]
        if(A[i,j]==1){
          store = rnorm(T,mu1[index1,index2],sigma1[index1,index2])
          store = rlaplace(T, mu1[index1, index2], sigma1[index1, index2] / sqrt(2))
          store = mu1[index1, index2] + rt(T, df) * sigma1[index1, index2]
          X[,4,k]=store
          E[i,j,] = E[j,i,] = store
        }
        if(A[i,j]==0) {
          store = rnorm(T,mu0[index1,index2],sigma0[index1,index2])
       #  store = rlaplace(T, mu0[index1, index2], sigma0[index1, index2] / sqrt(2))
       #  store = mu0[index1, index2] + rt(T, df) * sigma0[index1, index2]
          X[,4,k]=store
          E[i,j,] = E[j,i,] = store
        }
      }
    }
  }
  return(list(X = X,E = E))
}
# X=get_X(A,n,b,T,mu0,mu1,sigma0,sigma1)


# initial value
EM.gaussian <- function(X, N, n, G, hrho, hmu0, hmu1,hsigma0, hsigma1,maxitr=1000, delta=1e-5){
  Y=matrix(0, n, n)
  Y2=matrix(0, n, n)
  k=0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      k=k+1
      Y[i,j] = sum(X[,4,k])
      Y[j,i] = Y[i,j]
      Y2[i,j] = sum(X[,4,k]^2)
      Y2[j,i] = Y2[i,j]
    }
  }
  
  k=1
  conv=1
  X1=t(matrix(X[,4,],nrow =N[1,1])) 
  # Iterations 
  while(k<maxitr & conv>delta){
    # save current values
    hmu00 = hmu0
    hmu01 = hmu1
    hsigma00 = hsigma0
    hsigma01 = hsigma1
    hrho0=hrho
    
    # E-step
    th <- 1e-10  # thresholding value to avoid exact zero 

    weight=weight_gic(n, N[1,1], X1, hrho,  hmu0, hmu1,hsigma0, hsigma1)
    # M-step (mu, sigma and rho)
    for(i in 1:n){
        a0=(1-weight[i,-i])*(Y[i,-i]-N[i,-i]*hmu0[-i])/((hsigma0[i]*hsigma0[-i])^2)
        b0=(1-weight[i,-i])*N[i,-i]/((hsigma0[i]*hsigma0[-i])^2)
        hmu0[i]=sum(a0)/sum(b0)
        a1=weight[i,-i]*(Y[i,-i]-N[i,-i]*hmu1[-i])/((hsigma1[i]*hsigma1[-i])^2)
        b1=weight[i,-i]*N[i,-i]/((hsigma1[i]*hsigma1[-i])^2)
        hmu1[i]=sum(a1)/sum(b1)
        a0=(1-weight[i,-i])*(Y2[i,-i]-2*Y[i,-i]*(hmu0[i]+hmu0[-i])+N[i,-i]*(hmu0[i]+hmu0[-i])^2 )/((hsigma0[-i])^2)
        b0=(1-weight[i,-i])/((hsigma0[-i])^2)
        hsigma0[i]=sqrt(sum(a0)/sum(b0))
        a1=(1-weight[i,-i])*(Y2[i,-i]-2*Y[i,-i]*(hmu1[i]+hmu1[-i])+N[i,-i]*(hmu1[i]+hmu1[-i])^2 )/((hsigma1[-i])^2)
        b1=(1-weight[i,-i])/((hsigma1[-i])^2)
        hsigma1[i]=sqrt(sum(a1)/sum(b1))
        }
    hrho = sum(weight)/(n*(n-1))
    hsigma0[hsigma0<th]=th
    hsigma1[hsigma1<th]=th
    hrho[hrho<th]=th

    # convergence check
    k=k+1
    conv <- max(max(abs(hmu0-hmu00)), max(abs(hmu1-hmu01)), max(abs(hsigma0-hsigma00)) , max(abs(hsigma1-hsigma01)))
  }
  A=weight
  A[A>0.5]=1
  A[A<1]=0
  
  return(list(mu0=hmu0,mu1=hmu1,sigma0=hsigma0,sigma1=hsigma1,rho=hrho,A=A,iter=k))
}

#fit=EM.gaussian(X, N, n, G, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n) ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=1000, delta=1e-5)

initial.GTM.gau=function(X, N, n, G, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n)
                         ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=1000, delta=1e-5){
  fit=EM.gaussian(X, N, n, G, hrho, hmu0, hmu1,hsigma0, hsigma1,maxitr, delta)
  result=cbind(fit$mu0,fit$mu1,fit$sigma0,fit$sigma1)
  gfit=kmeans(result,centers=G)
  g0=gfit$cluster
  g1=gfit$centers
  # rho=unlist(fit$rho)
  # hrho0=matrix(rho,G,G)
  hrho0=matrix(fit$rho,G,G)
  hmu00=hmu01=hsigma00=hsigma01=matrix(0,G,G)
  
  for(i in 1:G){
   for (j in 1:G) {
     hmu00[i,j]=g1[i,1]+g1[j,1]
     hmu01[i,j]=g1[i,2]+g1[j,2]
     hsigma00[i,j]=abs(g1[i,3]*g1[j,3])
     hsigma01[i,j]=abs(g1[i,4]*g1[j,4])
   }
  }
  th=1e-5
  hsigma00[hsigma00<th]=th
  hsigma01[hsigma01<th]=th
  hrho0[hrho<th]=th
  A=fit$A
  
  return(list(group=g0,mu0=hmu00,mu1=hmu01,sigma0=hsigma00,sigma1=hsigma01,rho=hrho0,A=A ))
}

#fit=initial.GTM.gau(X, N, n, G, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n) ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=1000, delta=1e-5)

GTM.Gauss <- function(X, N, n, G, hA, hg, hrho, hmu0, hmu1
                      ,hsigma0, hsigma1, maxitr, delta){

  
  Y=matrix(0, n, n)
  Y2=matrix(0, n, n)
  k=0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      k=k+1
      Y[i,j] = sum(X[,4,k])
      Y[j,i] = Y[i,j]
      Y2[i,j] = sum(X[,4,k]^2)
      Y2[j,i] = Y2[i,j]
    }
  }
  
  k=1
  conv=1
  X1=t(matrix(X[,4,],nrow =N[1,1])) 
  # Iterations 
  while(k<maxitr & conv>delta){
    # save current values
    hg00 = hg
    hmu00 = hmu0
    hmu01 = hmu1
    hsigma00 = hsigma0
    hsigma01 = hsigma1
    hrho0=hrho
    hA0=hA
    
    # E-step
    th <- 1e-10    # thresholding value to avoid exact zero 
    r=numeric(G)
    for(g in 1:G){
      A=hA[which(hg==g),which(hg==g)]
      n1=nrow(A)
      if(is.null(n1)){
        r[g]=0
      }
      else{
        A_bar=sum(A)/(n1*(n1-1))
        A_hat=A-A_bar
        r[g]=r_g(n1,A_hat)
      } 
    }
    
    c=numeric(G)
    for(g in 1:G){
      A=hA[which(hg==g),which(hg==g)]
      n1=nrow(A)
      if(is.null(n1)){
        n1=1
      }
      if(n1>1){
        A_bar=sum(A)/(n1*(n1-1))
        A_hat=A-A_bar
        c[g]=r[g]*r_g(n,A_hat)
      }
      if(c[g]<0){
        c[g]=0
      }
    }
    c1=sum(c)+1

    weight=weight_gc(n, N[1,1], c1, hg-1, X1, hrho,  hmu0, hmu1,hsigma0, hsigma1)
    hA=weight
    hA[hA>0.5]=1
    hA[hA<1]=0
    # M-step (mu, sigma and rho)
    for(i in 1:G){
      for(j in 1:G){
        index1=which(hg==i)
        index2=which(hg==j)
        hmu0[i,j]=sum(Y[index1,index2] - Y[index1,index2]*weight[index1,index2])/sum(N[index1,index2]-N[index1,index2]*weight[index1,index2])
        hmu1[i,j]=sum(Y[index1,index2]*weight[index1,index2])/sum(N[index1,index2]*weight[index1,index2])
        hsigma0[i,j]=sqrt(abs(sum((Y2[index1,index2]-2*hmu00[i,j]*Y[index1,index2]+(hmu00[i,j]^2)*N[index1,index2])*(1 - weight[index1,index2]))/sum(N[index1,index2]-N[index1,index2]*weight[index1,index2])))
        hsigma1[i,j]=sqrt(abs(sum((Y2[index1,index2]-2*hmu01[i,j]*Y[index1,index2]+(hmu00[i,j]^2)*N[index1,index2])*weight[index1,index2])/sum(N[index1,index2]*weight[index1,index2])))
        if(i==j){
          n1=length(index1)  
          hrho[i,j]=sum(weight[index1,index2])/(n1*(n1-1))
        }
        else{
          hrho[i,j]=mean(weight[index1,index2]) 
        }
      }
    }
    hsigma0[hsigma0<th]=th
    hsigma1[hsigma1<th]=th
    hrho[hrho<th]=th
    a<-rep(0,n)
    hg=hg_gc(n,G,a,hg-1,r,hA,N,Y,Y2,weight,hrho,hmu0, hmu1,hsigma0,hsigma1)+1
    # convergence check
    k=k+1
    conv <- max(max(abs(hmu0-hmu00)), max(abs(hmu1-hmu01)), max(abs(hsigma0-hsigma00)) , max(abs(hsigma1-hsigma01)))
  }
  return(list(group=hg,mu0=hmu0,mu1=hmu1,sigma0=hsigma0,sigma1=hsigma1,rho=hrho, weight=weight,iter=k))  
}

#GTM <- GTM.Gauss(X, N, n, G, hA=fit$A, hg=fit$group, hrho=fit$rho, hmu0=fit$mu0/2, hmu1=fit$mu1 ,hsigma0=fit$sigma0/sqrt(2), hsigma1=fit$sigma1, maxitr=100, delta=1e-7)

#information criterion
IC <- function(X, N, n, G, hmu0, hmu1, hsigma0, hsigma1, hrho, hg){
  IC=0
  K=0
  Y=matrix(0, n, n)
  Y2=matrix(0, n, n)
  k=0
  for(i in 1:(n-1)){
    gi=hg[i]
    for(j in (i+1):n){
      k=k+1
      gj=hg[j]
      IC=IC-2*log((1-hrho[gi,gj])*prod(dnorm(X[,4,k],hmu0[gi,gj],hsigma0[gi,gj])) + hrho[gi,gj]*prod(dnorm(X[,4,k],hmu1[gi,gj],hsigma1[gi,gj])))
      K=K+N[i,j]
    }
  }
  IC=IC+(n+5*G*(G-1)/2)*log(K)
  return(IC)
}

get_g <- function(X,N,n){
  ICG=NULL
  for(g in 2:8){
    fit=initial.GTM.gau(X, N, n, G=g, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n)
                        ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=200, delta=1e-3)
    re1=tryCatch({GTM<-GTM.Gauss(X, N, n, G=g, hA=fit$A, hg=fit$group, hrho=fit$rho, hmu0=fit$mu0/2, hmu1=fit$mu1
                                 ,hsigma0=fit$sigma0/2, hsigma1=fit$sigma1, maxitr=200, delta=1e-3)},
                 error=function(e){
                   return(NULL)},
                 fillally=function(f){
                   return(GTM)
                 })
    k=1
    while(is.null(re1)){
      fit=initial.GTM.gau(X, N, n, G=g, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n)
                          ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=200, delta=1e-3)
      re1=tryCatch({GTM<-GTM.Gauss(X, N, n, G=g, hA=fit$A, hg=fit$group, hrho=fit$rho, hmu0=fit$mu0/2, hmu1=fit$mu1
                                   ,hsigma0=fit$sigma0/2, hsigma1=fit$sigma1, maxitr=200, delta=1e-3)},
                   error=function(e){
                     return(NULL)},
                   fillally=function(f){
                     return(GTM)
                   })
      k=k+1
    }
    print(k)
    ICG[g] <- IC(X, N, n,G=g,hmu0=GTM$mu0,hmu1=GTM$mu1,hsigma0=GTM$sigma0,hsigma1=GTM$sigma1,hrho=GTM$rho,hg=GTM$group)
    print("g:")
    print(g)
  }
  print(ICG)
  g=which.min(ICG)
  print("g_final:")
  return(g)
}
#g=get_g(X,N,n)

SE <- function(E,n,G){
  accuracy = 0.0
  mu = matrix(0, G, G)
  b=block_N(n,G)
  corr = pairwise_correlations(E)
  alpha_label = spectral_clustering(E[,,1],G)
  alpha = array(0, dim = c(n, G))
  for (i in 1:n){
    alpha[i,alpha_label[i]] = 0.7
    alpha[i,-alpha_label[i]] = (1-0.7)/(G-1)
  }
  res = EM_algorithm2(E,corr,alpha)
  gpn=matrix(0,nrow=G,ncol=G)
  for(j in 1:G){
    for(k in 1:G){
      low=(k-1)*round(n/G)+1
      up=min(k*round(n/G),n)
      gpn[j,k]=sum(res$z[low:up]==j)
    }
  }
  accuracy=lp.assign(gpn,direction = "max")$objval/n
  mu = res$mu
  return(list(acc=accuracy,mu = mu))
}

#simulation result
SR_Gauss <- function(N,T,n,G,rho,r,mu0,mu1,sigma0,sigma1,R){
  b=block_N(n,G)
  p1=NULL
  p0=NULL
  p11=NULL
  execution_time = 0
  p00=NULL
  gpr=NULL
  accuracy = 0
  for(i in 1:R){
    A=get_A(n,b,rho,r)
    res=get_X(A,n,b,T,mu0,mu1,sigma0,sigma1)
    X = res$X
    E = res$E
    g=3
    start_time <- Sys.time()
    fit=initial.GTM.gau(X, N, n, G=G, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n)
                        ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=100, delta=1e-5)
    re1=tryCatch({GTM<-GTM.Gauss(X, N, n, G=G, hA=fit$A, hg=fit$group, hrho=fit$rho, hmu0=fit$mu0/2, hmu1=fit$mu1
                            ,hsigma0=fit$sigma0/2, hsigma1=fit$sigma1, maxitr=100, delta=1e-3)},
                 error=function(e){
                   return(NULL)},
                 fillally=function(f){
                   return(GTM)
                 })
    k=1
    while(is.null(re1)){
      A=get_A(n,b,rho,r)
      res=get_X(A,n,b,T,mu0,mu1,sigma0,sigma1)
      X = res$X
      E = res$E
      fit=initial.GTM.gau(X, N, n, G=G, hrho=0.1, hmu0=rep(0,n), hmu1=rep(0,n)
                          ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=100, delta=1e-5)
      re1=tryCatch({GTM<-GTM.Gauss(X, N, n, G=G, hA=fit$A, hg=fit$group, hrho=fit$rho, hmu0=fit$mu0/2, hmu1=fit$mu1
                                   ,hsigma0=fit$sigma0/2, hsigma1=fit$sigma1, maxitr=100, delta=1e-3)},
                   error=function(e){
                     return(NULL)},
                   fillally=function(f){
                     return(GTM)
                   })
      k=k+1
    }
    #estimation error of A
    w=GTM$weight
    w1=w
    w2=w
    #sum(colSums((w>0.95)==TRUE))
    w1[w>0.8]=1
    w1[w<0.8]=0
    m1=w1+A
    p1[i]=length(m1[m1==2])/length(A[A==1])
    w2[w<0.2]=0
    w2[w>0.2]=1
    m2=w2+A
    p0[i]=length(m2[m2==0])/length(A[A==0])

    w1[w>0.9]=1
    w1[w<0.9]=0
    m1=w1+A
    p11[i]=length(m1[m1==2])/length(A[A==1])
    w2[w<0.1]=0
    w2[w>0.1]=1
    m2=w2+A
    p00[i]=length(m2[m2==0])/length(A[A==0])

    #estimation error of group
    gp=GTM$group

    gpn=matrix(0,nrow=G,ncol=G)
    for(j in 1:G){
      for(k in 1:G){
        low=(k-1)*round(n/G)+1
        up=min(k*round(n/G),n)
        gpn[j,k]=sum(gp[low:up]==j)
      }
    }
    gpr[i]=lp.assign(gpn,direction = "max")$objval/n
#    accuracy = accuracy+result$acc
    result_print=c(p1[i],p0[i],p11[i],p00[i],gpr[i])
    print(result_print)
    end_time <- Sys.time()
    execution_time <- execution_time + end_time - start_time
   # print(i)
  }
  p1=mean(p1)
  execution_time = execution_time/R
  p0=mean(p0)
  p11=mean(p11)
  p00=mean(p00)
  gpr=mean(gpr)
#  accuracy = accuracy/R
  results=c(p1,p0,p11,p00,gpr)
  return(list(results=results, execution_time=execution_time))
}


RA_Gauss <- function(N,T,n,E){
  G=3
  m=n*(n-1)/2
  X <- array(0, c(T, 4, m))
  k=0
  for(i in 1:n){
    if((i+1)<=n){
      for(j in (i+1):n){
        k=k+1
        X[,1,k]=i
        X[,2,k]=j
        X[,4,k]=E[i,j,]
      }
    }
  }
  fit=initial.GTM.gau(X, N, n, G=G, hrho=0.5, hmu0=rep(0.3,n), hmu1=rep(0.3,n)
                      ,hsigma0=rep(1,n), hsigma1=rep(1,n),maxitr=200, delta=1e-3)
  GTM<-GTM.Gauss(X, N, n, G=G, hA=fit$A, hg=fit$group, hrho=fit$rho, hmu0=fit$mu0/2, hmu1=fit$mu1
                               ,hsigma0=fit$sigma0/2, hsigma1=fit$sigma1, maxitr=200, delta=1e-3)
  
  gp=GTM$group
  return(GTM)
}



pairwise_correlations <- function(Y) {
  n <- dim(Y)[1]
  layer <- dim(Y)[3]
  v <- matrix(nrow = layer, ncol = n*n)
  for (H in 1:layer) {
    v[H, ] <- as.vector(Y[,,H])
  }
  corr_matrix <- cor((v))
  corr_matrix[is.na(corr_matrix)] <- 0
  diag(corr_matrix) <- 0
  return(corr_matrix)
}



spectral_clustering <- function(adj_matrix, K) {
  require(igraph)
  D <- diag(rowSums(adj_matrix))
  L <- D - adj_matrix
  eigen_values_vectors <- eigen(L)
  eigen_vectors <- eigen_values_vectors$vectors
  eigen_vectors_selected <- eigen_vectors[, order(eigen_values_vectors$values)[1:K]]
  norm <- sqrt(rowSums(eigen_vectors_selected^2))
  U <- eigen_vectors_selected / norm
  clusters <- kmeans(U, centers = K)$cluster
  return(clusters)
}




L_tilde2 <- function(Y, Y_hat, z ,rho_store, mu, alpha){
  N <- dim(alpha)[1]
  K <- dim(alpha)[2]
  M <- dim(Y)[3]
  result_L = matrix(0, nrow = N, ncol = K)
  xy_comb <- combn(1:N, 2)
  alpha_save = alpha
  A = array(0, dim = c(K, N, N))
  cx = array(0, dim = N)
  xx = array(0, dim = N)
  for (i in 1:K) {
    cx = array(0, dim = N)
    xx = array(0, dim = N)
    cx[which(z == i)] = 1
    xx[which(z == i)] = rho_store[i]
    A[i,,] = outer(cx,xx)
  }
  
  result_L <- compute_result_L(N, K, M, as.matrix(alpha_save), as.vector(rho_store), as.vector(Y_hat), as.matrix(z), as.matrix(xy_comb), as.matrix(mu), as.vector(Y))
  # for (i in 1:N) {
  #   for (q in 1:K) {
  #     alpha = alpha_save
  #     alpha[i,] = 0
  #     alpha[i,q] = 1
  #     log_likelihood <- 0
  #     for (m in 1:M) {
  #       log_likelihood <- log_likelihood + cal_A(K,Y[,,m],alpha,xy_comb,mu)
  #     }
  # for (m in 1:M) {
  # 
  #   sum_result = 0
  #   for(k in 1:K) {
  #     result = rho_store[k]*sum(outer(alpha[which(z == k), k],alpha[which(z == k), k])*outer(alpha[which(z == k), k],alpha[which(z == k), k])*Y_hat[which(z == k),which(z == k),m]*Y_hat[which(z == k),which(z == k),m])/2
  #     max_result = max(result,0)
  #     sum_result = sum_result + max_result
  #   }
  # 
  #   log_likelihood <- log_likelihood + log(1 + sum_result)
  # }
  #     result_L[i,q] = log_likelihood
  #   }
  # }
  return(result_L)
}


EM_algorithm2 <- function(Y, rho,alpha_init, epsilon = 5e-1){
  alpha <- alpha_init
  N <- dim(alpha)[1]
  K <- dim(alpha)[2]
  M <- dim(Y)[3]
  alpha_new <- matrix(1, nrow = N, ncol = K)
  alpha_new = alpha
  mat <- matrix(1, nrow = N, ncol = N)
  diag(mat) <- 0
  
  Y_hat = compute_Y_hat(Y)
  while(TRUE) {
    alpha <- alpha_new
    mu <- compute_mu(alpha, as.vector(Y), N, M, K)
    
    z <- rep(0, N)
    for (i in 1:N){
      z[i] = which(alpha[i,]==max(alpha[i,]),arr.ind=TRUE)
    }
    
    rho_store  = array(0, dim = K)
    count = array(0, dim = K)
    
    rho_store = cal_C(rho,rho_store,count,z,N,K)
    
    result_L = L_tilde2(Y,Y_hat,z,rho_store, mu,alpha)
    
    # Expectation
    alpha_new <- update_alpha(alpha, result_L, N, K)
    
    # Checking Convergence
    epsilon = N*0.1
    if (norm(alpha - alpha_new) < epsilon){
      break
    }
  }
  
  
  # Get the final cluster assignment
  z <- rep(0, N)
  for (i in 1:N){
    z[i] = which(alpha[i,]==max(alpha[i,]),arr.ind=TRUE)
  }
  
  return(list(z = z, alpha = alpha, mu = mu))
}



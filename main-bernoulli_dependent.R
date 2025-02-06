
block_N <- function(n,G){  #grouping
  b=NULL
  for(i in 1:G)
    if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
  else
  {b[((i-1)*round(n/G)+1):n]=i}
  return(b)
}

#generate A and E using SBM
get_AE <- function(T,n,b,rho,r,alpha,beta){
  A = E = matrix(0,n,n)
  E_copy = array(0, dim = c(n, n, T))
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
      if(A[i,j]==1){
        for (layer in 1:T) {
          store=rbinom(1,1,alpha[index1,index2])
          E[i,j]=E[i,j]+store
          E_copy[i,j,layer] = store
        }
      }
      else{
        for (layer in 1:T) {
          store=rbinom(1,1,beta[index1,index2])
          E[i,j]=E[i,j]+store
          E_copy[i,j,layer] = store
        }
      }
      A[j,i]=A[i,j]
      E[j,i]=E[i,j]
      for (layer in 1:T) {
        E_copy[j,i,layer] = E_copy[i,j,layer]
      }
    }
  }
  return(list(A=A,E=E,E_copy = E_copy))
}


EM.Bernoulli1 <- function(E,N, hrho,halpha,hbeta,maxitr=500,delta=1e-7){
  m=nrow(E)
  rho=hrho
  #initial
  
  k=1
  conv=1
  while(k<maxitr & conv>delta){
    halpha0 = halpha
    hbeta0 = hbeta
    
    #E step
    Q <- matrix(0, m, m)
    for(i in 2:m){
      for(j in 1:(i-1)){
        a=rho*(halpha[i]^E[i,j])*((1-halpha[i])^(N[i,j]-E[i,j]))*(halpha[j]^E[j,i])*((1-halpha[j])^(N[j,i]-E[j,i]))
        b=(1-rho)*(hbeta[i]^E[i,j])*((1-hbeta[i])^(N[i,j]-E[i,j]))*(hbeta[j]^E[j,i])*((1-hbeta[j])^(N[j,i]-E[j,i]))
        Q[i,j] =a/(a+b)
        Q[j,i]=Q[i,j]
      }
    }
    
    #M step
    th <- 1e-20  # thresholding value to avoid exact zero 
    halpha = colSums(E*Q) / colSums(N*Q)
    halpha[halpha<th] = th
    hbeta = colSums(E-E*Q) / colSums(N-N*Q)
    hbeta[hbeta<th] = th
    rho = sum(Q)/(m*(m-1))
    conv <- max(c(abs(halpha0-halpha),abs(hbeta0 - hbeta)))
    k=k+1
  }
  return(list(alpha=halpha,beta=hbeta,rho=rho,Q=Q,iter=k))
}

initial.GTM1=function(E,N,G, hrho,halpha,hbeta, maxitr=500,delta=1e-7){
  m=nrow(E)
  fit= EM.Bernoulli1(E, N, hrho, halpha, hbeta,maxitr, delta)
  gfit=kmeans(fit$alpha,centers=G)
  g0=gfit$cluster
  gfit$centers
  rho0 <- matrix(fit$rho,G,G)
  alpha0 <-gfit$centers%*%t(gfit$centers)
  b=numeric(G)
  for(i in 1:G){
    b[i]=mean(fit$beta[which(g0==i)])
  }
  b[b<1e-3]=1e-3
  beta0 <- b%*%t(b)
  return(list(group=g0,alpha=alpha0,beta=beta0,rho=rho0 ))
}

#fit=initial.GTM(E,N,G=3,hrho=0.03,halpha=rep(0.6,n),hbeta=rep(0.06,n)) 


GTM.Bernoulli1 <- function(E, N, G,hrho,halpha,hbeta,hg, maxitr=500,delta= 1e-7 ){
  m <- nrow(E)
  k=1
  conv=1
  
  # Iterations 
  while(k<maxitr & conv>delta){
    # save current values
    halpha0 = halpha
    hbeta0 = hbeta
    hrho0 = hrho
    hg0 = hg
    
    # E-step
    th1=1e-5
    th =1e-10  # thresholding value to avoid exact zero 
    # weight <- matrix(0, m, m)
    # for(i in 2:m){
    #   for(j in 1:(i-1)){
    #     g1=hg[i]
    #     g2=hg[j]
    #     a=hrho[g1,g2]*(halpha[g1,g2]^E[i,j])*((1-halpha[g1,g2])^(N[i,j]-E[i,j]))
    #     b=(1-hrho[g1,g2])*(hbeta[g1,g2]^E[i,j])*((1-hbeta[g1,g2])^(N[i,j]-E[i,j]))
    #     weight[i,j]=a/(a+b)
    #     weight[j,i]=weight[i,j]
    #   }
    # }
    weight=weight_c1(m,hg-1, E, N, hrho, halpha, hbeta)
    # M-step (alpha, beta and rho)
    for(i in 1:G){
      for(j in 1:G){
        index1=which(hg==i)
        index2=which(hg==j)
        halpha[i,j]=sum(E[index1,index2]*weight[index1,index2])/sum(N[index1,index2]*weight[index1,index2])
        hbeta[i,j]=sum(E[index1,index2]-E[index1,index2]*weight[index1,index2])/sum(N[index1,index2]-N[index1,index2]*weight[index1,index2])
        if(i==j){
          n1=length(index1)  
          hrho[i,j]=sum(weight[index1,index2])/(n1*(n1-1))
        }
        else{
          hrho[i,j]=mean(weight[index1,index2]) 
        }
      }
    }
    halpha[halpha<th1]=th1
    hbeta[hbeta<th]=th
    
    # M-step (grouping variable)
    # val <- matrix(0,m,G)
    # for(i in 1:m){
    #   for(g in 1:G){
    #     for(j in 1:m){
    #       gj=hg[j]
    #       if(j!=i){
    #         val[i,g]=val[i,g]+(1-weight[i,j])*(E[i,j] * log(hbeta[g,gj])+(N[i,j]-E[i,j]) * log(1-hbeta[g,gj]) +log(1-hrho[g,gj]))+weight[i,j]*(E[i,j] * log(halpha[g,gj]) +(N[i,j]-E[i,j]) * log(1-halpha[g,gj]) + log(hrho[g,gj]))
    #       }
    #     }
    #   }
    #   hg[i]=which.max(val[i,])
    # }
    
    hg=hg_c1( m,G,hg-1,E,N,weight,hrho,  halpha,hbeta)+1
    # convergence check
    k=k+1
    conv <- max(max(abs(halpha-halpha0)), max(abs(hbeta-hbeta0)))
  }
  
  return(list(group=hg,alpha=halpha,beta=hbeta,rho=hrho, weight=weight,iter=k))  
}


#information criterion
IC1 <- function(N,E,G,halpha,hbeta,hrho,hg){
  IC=0
  n=length(hg)
  for(j in 2:n){
    gj=hg[j]
    for(i in 1:(j-1)){
      gi=hg[i]
      IC=IC-2*log((1-hrho[gi,gj])*hbeta[gi,gj]^E[i,j]*(1-hbeta[gi,gj])^(N[i,j]-E[i,j]) + hrho[gi,gj]*halpha[gi,gj]^E[i,j]*(1-halpha[gi,gj])^(N[i,j]-E[i,j]))
    }
  }
  K=(sum(N)-sum(diag(N)))/2
  IC=IC+(n+3*(G-1)*G/2)*log(K)
  return(IC)
}

#IC(N,E,G=3,halpha=GTM$alpha,hbeta=GTM$beta,hrho=GTM$rho,hg=GTM$group)

get_g1 <- function(E,N,n){
  ICG=NULL
  for(k in 1:3){
    fit=initial.GTM1(E,N,G=k,hrho=0.1,halpha=rep(0.3,n),hbeta=rep(0.06,n))
    GTM=GTM.Bernoulli1(E, N, G=k,hrho=fit$rho,halpha=fit$alpha,hbeta=fit$beta,hg=fit$group,maxitr=5000,delta= 1e-7 )
    ICG=c(ICG,IC1(N,E,G=k,halpha=GTM$alpha,hbeta=GTM$beta,hrho=GTM$rho,hg=GTM$group))
  }
  g=which.min(ICG)
  return(g)
}


EM.Bernoulli <- function(E,N, hrho,halpha,hbeta,maxitr=500,delta=1e-7){
  m=nrow(E)
  rho=hrho
  #initial
  
  k=1
  conv=1
  while(k<maxitr & conv>delta){
    halpha0 = halpha
    hbeta0 = hbeta
    
    #E step
    Q <- matrix(0, m, m)
    for(i in 2:m){
      for(j in 1:(i-1)){
        a=rho*(halpha[i]^E[i,j])*((1-halpha[i])^(N[i,j]-E[i,j]))*(halpha[j]^E[j,i])*((1-halpha[j])^(N[j,i]-E[j,i]))
        b=(1-rho)*(hbeta[i]^E[i,j])*((1-hbeta[i])^(N[i,j]-E[i,j]))*(hbeta[j]^E[j,i])*((1-hbeta[j])^(N[j,i]-E[j,i]))
        Q[i,j] =a/(a+b)
        Q[j,i]=Q[i,j]
      }
    }
    
    #M step
    th <- 1e-5  # thresholding value to avoid exact zero 
    halpha = colSums(E*Q) / colSums(N*Q)
    halpha[halpha<th] = th
    halpha[halpha>(1-th)] = 1-th
    hbeta = colSums(E-E*Q) / colSums(N-N*Q)
    hbeta[hbeta<th] = th
    hbeta[hbeta>(1-th)] = 1-th
    rho = sum(Q)/(m*(m-1))
    conv <- max(c(abs(halpha0-halpha),abs(hbeta0 - hbeta)))
    k=k+1
  }
  #A
  A=Q
  A[A>0.5]=1
  A[A<1]=0
  return(list(alpha=halpha,beta=hbeta,rho=rho,Q=Q,hA=A,iter=k))
}

initial.GTM=function(E,N,G, hrho,halpha,hbeta, maxitr=500,delta=1e-7){
  m=nrow(E)
  fit= EM.Bernoulli(E, N, hrho, halpha, hbeta,maxitr, delta)
  gfit=kmeans(fit$alpha,centers=G)
  g0=gfit$cluster
  g1=gfit$centers
  rho0 <- matrix(fit$rho,G,G)
  alpha0 <-g1%*%t(g1)
  b=numeric(G)
  for(i in 1:G){
    b[i]=mean(fit$beta[which(g0==i)])
  }
  b[b<1e-3]=1e-3
  beta0 <- b%*%t(b)
  A=fit$hA
  return(list(group=g0,alpha=alpha0,beta=beta0,rho=rho0,hA=A ))
}


GTM.Bernoulli <- function(E, N, G,hA,hrho,halpha,hbeta,hg ){
  m <- nrow(E)
  k=1
  conv=1
  maxitr=1000
  delta=1e-7
  
  # Iterations 
  if(k<maxitr && conv>delta){
    # save current values
    halpha0 = halpha
    hbeta0 = hbeta
    hrho0 = hrho
    hg0 = hg
    hA0=hA
    
    # E-step
    th1=1e-5
    th =1e-10  # thresholding value to avoid exact zero 
    r=numeric(G)
    for(g in 1:G){
      A=hA[which(hg==g),which(hg==g)]
      n=nrow(A)
      if(is.null(n)){
        r[g]=0
      }
      else{
        A_bar=sum(A)/(n*(n-1))
        A_hat=A-A_bar
        r[g]=r_g(n,A_hat)
        # r[g]=(sum(A_hat)^2-2*sum(A_hat*A_hat))/((n*(n-1))^2-2*n)
      } 
    }

    c=numeric(G)
    for(g in 1:G){
      A=hA[which(hg==g),which(hg==g)]
      n=nrow(A)
      if(is.null(n)){
        n=1
      }
      if(n>1){
        A_bar=sum(A)/(n*(n-1))
        A_hat=A-A_bar
        c[g]=r[g]*r_g(n,A_hat)
      }
      if(c[g]<0){
        c[g]=0
      }
    }
    c1=sum(c)+1
    weight=weight_c(m,c1,hg-1,E,N,hrho,halpha,hbeta)
    # M-step (A,alpha, beta and rho)
    hA=weight
    hA[hA>0.5]=1
    hA[hA<1]=0
    for(i in 1:G){
      for(j in 1:G){
        index1=which(hg==i)
        index2=which(hg==j)
        halpha[i,j]=sum(E[index1,index2]*weight[index1,index2])/sum(N[index1,index2]*weight[index1,index2])
        hbeta[i,j]=sum(E[index1,index2]-E[index1,index2]*weight[index1,index2])/sum(N[index1,index2]-N[index1,index2]*weight[index1,index2])
        if(i==j){
          n1=length(index1)
          hrho[i,j]=sum(weight[index1,index2])/(n1*(n1-1))
        }
        else{
          hrho[i,j]=mean(weight[index1,index2]) 
        }
      }
    }
    
    halpha[halpha<th1]=th1
    halpha[halpha>(1-th1)]=1-th1
    hbeta[hbeta<th]=th
    hbeta[hbeta>(1-th)]=1-th
    hrho[hrho<th1]=th1
    hrho[hrho>(1-th1)]=1-th1
    
    a<-rep(0,m)
    hg=hg_c(m,G,a,hg-1,r,hA,E,N,weight,hrho,halpha,hbeta)+1
    #hg=hg_c( m,G,hg-1,E,N,weight,hrho,halpha,hbeta,r,hA)+1
    # convergence check
    k=k+1
    conv <- max(max(abs(halpha-halpha0)), max(abs(hbeta-hbeta0)))
  }
  
  return(list(group=hg,alpha=halpha,beta=hbeta,rho=hrho, weight=weight,r=r,A=hA,iter=k))  
}


#information criterion
IC <- function(N,E,G,halpha,hbeta,hrho,hg,hA,weight){
  IC=0
  m=length(hg)
  
  r=numeric(G)
  for(g in 1:G){
    A=hA[which(hg==g),which(hg==g)]
    n=nrow(A)
    if(is.null(n)){
      n=1
      r[g]=0
    }
    if(n>1){
       A_bar=sum(A)/(n*(n-1))
       A_hat=A-A_bar
      # a1=b1=0
      # for(i in 2:n){
      #   for(j in 1:(i-1)){
      #     for(u in 2:n){
      #       for(v in 1:(u-1)){
      #         if(i!=u||j!=v){
      #           a1=a1+A_hat[i,j]*A_hat[u,v]
      #           b1=b1+1
      #         }
      #       }
      #     }
      #   }
      # }
      # if(b1==0){
      #   r[g]=0
      # }
      # else{
      #   rg=a1/b1
      #   if(rg<0){
      #     r[g]=0
      #   }
      #   else{
      #     r[g]=rg
      #   }
      # }
      r[g]=r_g(n,A_hat)
    }
  }
  
  # c=0
  # for(g in 1:G){
  #   A=hA[which(hg==g),which(hg==g)]
  #   n=nrow(A)
  #   if(is.null(n)){
  #     n=1
  #   }
  #   if(n>1){
  #     A_bar=sum(A)/(n*(n-1))
  #     A_hat=A-A_bar
  #     for(i in 2:n){
  #       for(j in 1:(i-1)){
  #         for(u in 2:n){
  #           for(v in 1:(u-1)){
  #             if(i!=u||j!=v){
  #               c=c+r[g]*A_hat[i,j]*A_hat[u,v]
  #             }
  #           }
  #         }
  #       }
  #     }
  #   }
  # }
  # c=c+1
  c=numeric(G)
  for(g in 1:G){
    A=hA[which(hg==g),which(hg==g)]
    n=nrow(A)
    if(is.null(n)){
      n=1
    }
    if(n>1){
      A_bar=sum(A)/(n*(n-1))
      A_hat=A-A_bar
      #for(i in 2:n){
      #  for(j in 1:(i-1)){
      #    for(u in 2:n){
      #      for(v in 1:(u-1)){
      #        if(i!=u||j!=v){
      #          c=c+r[g]*A_hat[i,j]*A_hat[u,v]
      #        }
      #      }
      #    }
      #  }
      #}
      c[g]=which_c(n,r[g],A_hat)
    }
    if(c[g]<0){
      c[g]=0
    }
  }
  c1=sum(c)+1
  for(j in 2:m){
    gj=hg[j]
    for(i in 1:(j-1)){
      gi=hg[i]
      IC=IC-2*((1-hA[i,j])*log((1-hrho[gi,gj])*hbeta[gi,gj]^E[i,j]*(1-hbeta[gi,gj])^(N[i,j]-E[i,j]))+hA[i,j]*log(hrho[gi,gj]*halpha[gi,gj]^E[i,j]*(1-halpha[gi,gj])^(N[i,j]-E[i,j]))+log(c1))
    }
  }
  # IC=which_IC(m,c1,hg,E,N,hrho,halpha,hbeta,weight)
  K=(sum(N)-sum(diag(N)))/2
  IC=IC+(n+3*(G+1)*G/2)*log(K)
  #IC=IC+(3*(G+1)*G/2)*log(K)
  return(IC)
}

#IC(N,E,G,halpha=GTM$alpha,hbeta=GTM$beta,hrho=GTM$rho,hg=GTM$group,hA=GTM$A)

get_g <- function(E,N,n){
  ICG=NULL
  for(kg in 1:3){
    fit=initial.GTM(E,N,G=kg,hrho=0.8,halpha=rep(0.6,n),hbeta=rep(0.06,n))
    GTM=GTM.Bernoulli(E, N, G=kg,hA=fit$hA,hrho=fit$rho,halpha=fit$alpha,hbeta=fit$beta,hg=fit$group )
    ICG=c(ICG,IC(N,E,G=kg,halpha=GTM$alpha,hbeta=GTM$beta,hrho=GTM$rho,hg=GTM$group,hA=GTM$A,weight=GTM$weight))
  }
  g=which.min(ICG)
  return(g)
}

# t=proc.time()
#g=get_g(E,N,n)
# proc.time()-t


SR1 <- function(T,N,n,G,r,rho,alpha,beta,R){
  b=block_N(n,G)
  m=NULL
  p1=NULL
  p0=NULL
  gpr=NULL
  rIC=NULL
  for(i in 1:R){
    fit0=get_AE(T,n,b,rho,r,alpha,beta)
    A=fit0$A
    E=fit0$E
    g=get_g1(E,N,n)
    fit=initial.GTM1(E,N,G=G,hrho=0.1,halpha=rep(0.3,n),hbeta=rep(0.06,n))
    GTM=GTM.Bernoulli1(E, N, G=G,hrho=fit$rho,halpha=fit$alpha,hbeta=fit$beta,hg=fit$group,maxitr=5000,delta= 1e-7 )
    
    #estimation error of A
    w=GTM$weight
    m[i]=sum(abs(w-A))/(n*(n-1))
    w1=w
    w2=w
    #sum(colSums((w>0.95)==TRUE))
    w1[w1>0.75]=1
    w1[w1<0.75]=0
    m1=w1+A
    p1[i]=length(m1[m1==2])/length(A[A==1])
    
    w2[w2<0.25]=0
    w2[w2>0.25]=1
    m2=w2+A
    p0[i]=length(m2[m2==0])/length(A[A==0])
    
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
    
    #estimation of group number
    if(g==G) rIC[i]=1 else rIC[i]=0
  }
  m=mean(m)
  p1=mean(p1)
  p0=mean(p0)
  gpr=mean(gpr)
  rIC=mean(rIC)
  results=c(m,p1,p0,gpr,rIC)
  return(results)
}





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
SR <- function(T,N,n,G,r,rho,alpha,beta,R){
  b=block_N(n,G)
  accuracy = 0
  execution_time = 0
  gpr=NULL
  for(i in 1:R){
    fit0=get_AE(T,n,b,rho,r,alpha,beta)
    E=fit0$E
    g=3
    start_time <- Sys.time()
    fit=initial.GTM(E,N,G=G,hrho=0.8,halpha=rep(0.6,n),hbeta=rep(0.06,n))
    GTM=GTM.Bernoulli(E, N, G=G,hA=fit$hA,hrho=fit$rho,halpha=fit$alpha,hbeta=fit$beta,hg=fit$group)
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
    print(gpr[i])
    end_time <- Sys.time()
    execution_time <- execution_time + end_time - start_time
  }
  gpr=mean(gpr)
  execution_time = execution_time/R
  return(list(gpr=gpr,execution_time=execution_time))
}


block_N <- function(n,G){  #grouping
  b=NULL
  for(i in 1:G)
    if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
  else
  {b[((i-1)*round(n/G)+1):n]=i}
  return(b)
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


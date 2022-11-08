library(MASS)
library(Matrix)
library(glmnet)
library(survival)
library(matrixStats)


#### Data generating process ####
# (H^* has diagnal structure)
dgp <- function(seed_=0, n_=100,r=50,p=50,q=4,n_main=8,n_interactions=8,
                 rho1=0.3,rho2=0.3,rho3=0.3,
                 n_tes=0
)
{
  # Params:
  #  seed_: sample size 
  #  n_: sample size ; 
  #  r: num of gene features
  #  p: num of imaging features; 
  #  q: num of environment features
  #  rho1, rho2, rho3: covariance param fo Z (G factors), U(error terms when generating I factors) and E (E factors)
  #  n_tes: number of test samples
  
  # Return:
  #  X: G factors; 
  #  E: E factors;
  #  Z: I factors;
  #  Y: censored survival time
  #  N_: censored indicator (1 if dead; 0 if censored)
  #  phi1_all: coefficients for G-E interaction model
  
  n = n_+n_tes
  #### initialization steps ####
  # initailize H (r*p)
  H_ = as(diag(1,r,p),"sparseMatrix")
    
  # initialize coefficients
  set.seed(1)  # seed for random coefficients
  idxs = 1:n_main
  beta = sparseVector(runif(n_main,0.2,0.4),
                      idxs,r)
  idxs2 = sample(0:(n_main*q-1),n_interactions,replace=F)
  j = idxs[idxs2%/%q + 1]
  i = idxs2%%q+1
  gamma = sparseMatrix(i=i,j=j,
                       x=runif(n_interactions,0.2,0.4),
                       dims=c(q,r))  
  alpha = runif(q,0,0.3)
  phi1 = c(beta,as.vector(t(gamma)))
  
  #### Generate data ####
  set.seed(seed_) # seed for generate data
  # generate gene factor (Z); Environmental factor (E); Imaging features (X)
  Sigma1 = rho1^(abs(outer(1:r,1:r,"-")))
  Sigma2 = rho2^(abs(outer(1:p,1:p,"-")))
  Sigma3 = rho3^(abs(outer(1:q,1:q,"-")))
  Z = mvrnorm(n,mu=rep(0,r),Sigma = Sigma1)
  U = mvrnorm(n,mu=rep(0,p),Sigma = Sigma2)*0.2
  X = as.matrix(Z%*%H_) + U
  E = mvrnorm(n,mu=rep(0,q),Sigma = Sigma3)
  
  # generate survival time (T_), censored survival time (Y), survival indicator (N_)
  # generate survival time
  interaction_term = rep(NA,n)
  for(i in 1:n){
    interaction_term[i] = sum(gamma*(outer(E[i,],Z[i,],"*")))
  }
  lambda = exp(1+(Z%*%beta)[,1]+(E%*%alpha)[,1]+interaction_term)
  T_ = rexp(n,rate=lambda)
  
  # generate censoring (use exponential distribution to generate censoring time C)
  # the parameter of this exponential distribution is set to met the cesoring rate requirement
  f <- function(rate) (1-mean(exp(-rate*T_)))-0.4
  c_rate = uniroot(f,interval=c(1e-10,1e10))$root
  C = rexp(n,c_rate)
  
  # generate censored survival time 
  Y = pmin(T_,C)
  
  # generate censored indicator
  N_ = as.integer(T_<=C)  # = 0 if censoring
  
  phi1_all = c(phi1,alpha)
  res = list(X=X[1:n_,],E=E[1:n_,],Z=Z[1:n_,],Y=Y[1:n_],N_=N_[1:n_],
             phi1_all=phi1_all
             )
  
  if(n_tes>0){
    res$X_tes = X[(n_+1):(n_+n_tes),]
    res$E_tes = E[(n_+1):(n_+n_tes),]
    res$Z_tes = Z[(n_+1):(n_+n_tes),]
    res$Y_tes = Y[(n_+1):(n_+n_tes)]
    res$T_tes = T_[(n_+1):(n_+n_tes)]
    res$N_tes = N_[(n_+1):(n_+n_tes)]
  }
  return(res)
}


#### Functions for the proposed method ####
#### auxiliary functions(for function ARMI_GE) ####
# estimate F 
F_est <- function(X,Z)
{
  # regress X (n*p) on Z (n*r)
  # use BIC to choose tunings
  n = dim(X)[1]
  p = dim(X)[2]
  r = dim(Z)[2]
  F_ = as(matrix(0,p,r),"sparseMatrix")
  for(j in 1:r){
    model = glmnet(X,Z[,j],nlambda = 10,lambda.min.ratio = 0.02)
    beta = model$beta
    bic = model$nulldev*(1-model$dev.ratio)+apply(model$beta!=0,2,sum)*log(n)
    idx = which.min(bic)
    F_[,j] = beta[,idx]
  }
  return(F_) # p*r
}


# partial_loglikelihood function
partial_loglike <- function(phi,N_ord,A_ord,Y_ord_raw,part=NULL)
{
  n = length(Y_ord_raw)
  tmp = (A_ord%*%phi)[,1]
  tmp2 = cumsum(exp(tmp))
  for(i in (n-1):1){
    if(Y_ord_raw[i]==Y_ord_raw[i+1]){
      tmp2[i] = tmp2[i+1]
    }
  }
  if(is.null(part)){
    pl = mean(N_ord*(tmp-log(tmp2)))
  }else{
    pl = mean((N_ord*(tmp-log(tmp2)))[part]) 
  }
  return(pl)
}


# first order derivative for partial loglikelihood (gradient)
partial_loglike.first <- function(phi,N_ord,A_ord,Y_ord_raw)
{
  n = length(Y_ord_raw)
  tmp = exp((A_ord%*%phi)[,1])
  tmp2 = cumsum(tmp)
  tmp3 = colCumsums(A_ord*tmp)
  for(i in (n-1):1){
    if(Y_ord_raw[i]==Y_ord_raw[i+1]){
      tmp2[i] = tmp2[i+1]
      tmp3[i,] = tmp3[i+1,]
    }
  }
  pl.first = apply((A_ord-tmp3/tmp2)*N_ord,2,mean)
  return(pl.first)
}

# update step I
update.I_1 <- function(phi1_all,phi2,N_ord,A_ord_all,Y_ord_all,F_,
                       rho1,lambda3,Tau1,p,q,r)
{
  idx1 = 1:(r*(q+1))
  n = dim(A_ord_all)[1]
  D_t = matrix(phi2,p,q+1)
  Tau1_ = as.vector(Tau1)
  obj <- function(phi1_all){
    B_t = matrix(phi1_all[idx1],r,q+1)
    objective = -partial_loglike(phi1_all,N_ord,A_ord_all,Y_ord_all)+
      rho1/2*sum((phi1_all[idx1]-Tau1_)^2)+
      lambda3*sum((D_t-F_%*%B_t)^2)
    return(objective)
  }
  
  obj.first <- function(phi1_all){
    B_t = matrix(phi1_all[idx1],r,q+1)
    g = rho1*(phi1_all[idx1]-Tau1_) - 2*lambda3*as.vector(t(F_)%*%(D_t-F_%*%B_t))
    gr = -partial_loglike.first(phi1_all,N_ord,A_ord_all,Y_ord_all)
    gr[idx1] = gr[idx1]+g
    return(gr)
  }
  return(optim(phi1_all,fn=obj,gr=obj.first,method="BFGS")$par)
}

update.I_2 <- function(phi1,phi2_all,N_ord,A_ord_all,Y_ord_all,F_,
                       rho2,lambda3,Tau2,p,q,r)
{
  idx2 = 1:(p*(q+1))
  n = dim(A_ord_all)[1]
  B_t = matrix(phi1,r,q+1)
  Tau2_ = as.vector(Tau2)
  obj <- function(phi2_all){
    D_t = matrix(phi2_all[idx2],p,q+1)
    objective = -partial_loglike(phi2_all,N_ord,A_ord_all,Y_ord_all)+
      rho2/2*sum((phi2_all[idx2]-Tau2_)^2)+
      lambda3*sum((D_t-F_%*%B_t)^2)
    return(objective)
  }
  
  obj.first <- function(phi2_all){
    D_t = matrix(phi2_all[idx2],p,q+1)
    g = rho2*(phi2_all[idx2]-Tau2_) + 2*lambda3*as.vector(D_t-F_%*%B_t)
    gr = -partial_loglike.first(phi2_all,N_ord,A_ord_all,Y_ord_all)
    gr[idx2] = gr[idx2]+g
    return(gr)
  }
  return(optim(phi2_all,fn=obj,gr=obj.first,method="BFGS")$par)
}


# soft-thresholding function
St <- function(x,lambda)
{
  sign(x)*pmax(abs(x)-lambda,0)
}

# first order derivative for penalty function(non-zero region)
St.d <- function(x,lambda)
{
  return(sign(x)*lambda)
}

# updatate step II
update.II <- function(a,lambda1)
{ 
  q = length(a)-1
  st_res = St(a[2:(q+1)],lambda1)
  tmp = c(a[1],st_res)
  # all group is 0
  if(sqrt(sum(tmp^2))<=lambda1*sqrt(q+1)){
    return(rep(0,q+1))
  }
  # iterative to find coefficients within group
  h = a
  for(iter in 1:100){
    h0 = h
    h_norm = sqrt(sum(h0^2))
    h[1] = a[1]/(1+lambda1*sqrt(q+1)/h_norm)
    h[2:(q+1)] = st_res/(1+lambda1*sqrt(q+1)/h_norm)
    if(mean((h-h0)^2)<=1e-3){
      break
    }
  }
  return(h)
}

#### functions of the proposed method ####
IA_GE <- function(X,E,Z,Y,N_,lambda1,lambda2,lambda3,rho=1,
                    phi1_all=NULL,phi2_all=NULL,
                    maxit1=200,
                    epsilon1=5e-4,epsilon2=5e-4
                    )
{
  # Params:
  #  X: G factors; 
  #  E: E factors;
  #  Z: I factors;
  #  Y: censored survival time;
  #  N_: censored indicator (1 if dead; 0 if censored);
  #  lamdba1-lambda3: tuning parameters ;
  #  rho: paramter for ADMM algorithm (augment factor);
  #  phi1_all: initial value for `phi1_all`;
  #  phi2_all: initial value for `phi2_all`;
  #  maxit1: maximum number of iterations for ADMM algorithm;
  #  epsilon1: stopping criteria for ADMM (primal feasibility) 
  #  epsilon2: stopping criteria for ADMM (dual feasibility) 

  # Return:
  # phi1_all: estimated coefficients in G-E model 
  # phi2_all: estimated coefficients in I-E model 
  # iters: number of iteration for ADMM algorithm
  
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(E)[2]
  r = dim(Z)[2]
  idx1 = 1:(r*(q+1))
  idx2 = 1:(p*(q+1))
  
  #### initialize #### 
  # estimate F_
  F_ = F_est(X,Z)

  # initialize parameters
  if(is.null(phi1_all)){
    phi1_all = runif(r*(q+1)+q,-0.5,0.5)
  }
  if(is.null(phi2_all)){
    phi2_all = rep(0,(q+1)*p+q)
    phi2_all[idx2] = as.vector(F_%*%matrix(phi1_all[idx1],r,q+1))
    phi2_all[-idx2] = phi1_all[-idx1]
  }
  
  # vectorize variables
  A1 = Z
  A2 = X
  for(k in 1:q){
    W = E[,k]*Z
    A1 = cbind(A1,W)
    V = E[,k]*X
    A2 = cbind(A2,V)
  }
  A1_all = cbind(A1,E)
  A2_all = cbind(A2,E)
  # order the data (follow decring order of Y)
  Y_ord = order(Y,decreasing = T)
  Y_ord_raw = Y[Y_ord]
  N_ord = N_[Y_ord]
  A1_ord_all = A1_all[Y_ord,]
  A2_ord_all = A2_all[Y_ord,]
  # Additional variable
  H = matrix(phi1_all[idx1],q+1,r,byrow = T)
  Phi = matrix(phi2_all[idx2],q+1,p,byrow = T)
  
  # Lagrangian Multiplier
  Pi1 = matrix(0,q+1,r)
  Pi2 = matrix(0,q+1,p)
  
  #### model fitting(ADMM) ####
  rho1 = rho
  rho2 = rho
  
  for(i1 in 1:maxit1){
    H0 = H      # these quantities are used to calculate stoping critiria(dual feasibility)
    Phi0 = Phi
    #### Step I ####
    Tau1 = t(H+Pi1)    # r*(q+1)
    Tau2 = t(Phi+Pi2)  # p*(q+1)
    
    for(i2 in 1:5){
      phi1_all_0 = phi1_all
      phi2_all_0 = phi2_all
      phi1_all = update.I_1(phi1_all,phi2_all[idx2],N_ord,A1_ord_all,Y_ord_raw,F_,
                            rho1,lambda3,Tau1,p,q,r)
      phi2_all = update.I_2(phi1_all[idx1],phi2_all,N_ord,A2_ord_all,Y_ord_raw,F_,
                            rho2,lambda3,Tau2,p,q,r)
      if(mean((phi1_all-phi1_all_0)^2)+mean((phi2_all-phi2_all_0)^2)<1e-5){
        break
      }
    }
    B = matrix(phi1_all[idx1],q+1,r,byrow = T)
    D = matrix(phi2_all[idx2],q+1,p,byrow = T)
    
    #### Step II ####
    for(j in 1:r){
      H[,j] = update.II(B[,j]-Pi1[,j],lambda1/rho1)
    }
    for(j in 1:p){
      Phi[,j] = update.II(D[,j]-Pi2[,j],lambda2/rho2)
    }
    
    #### step III ####
    Pi1 = Pi1+H-B
    Pi2 = Pi2+Phi-D
    
    #### stoping criteria ####
    stop.primal = mean((H-B)^2)+mean((Phi-D)^2)
    stop.dual = mean((H-H0)^2)+mean((Phi-Phi0)^2)
    if((stop.primal<=epsilon1)&(stop.dual<=epsilon2)){
      break
    }
  }
  if(i1==maxit1){
    warning(sprintf("ADMM desn't converge in %d iterations.",maxit1))
  }
  
  #### output result ####
  phi1_all[idx1] = as.vector(t(H))
  phi2_all[idx2] = as.vector(t(Phi))
  
  result = list(phi1_all=phi1_all,phi2_all=phi2_all,iters=i1)
  return(result)
}




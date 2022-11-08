source("function.R")

#### generate data ####
data = dgp(seed_=0, n_=300, r=50, p=50, q=4, n_main=8, n_inter=8,
           rho1=0.3,rho2=0.3,rho3=0.3
           )

Z = data$Z # G factors
X = data$X # E factors
E = data$E # I factors
Y = data$Y # censored survival time 
N_ = data$N_ # censored indicator (1: not censored; 0: censored)

#### fit the model #### 
set.seed(42)
model = IA_GE(X, E, Z, Y, N_, 0.05, 0.05, 0.1)
phi1_all = model$phi1_all  # estimated coefficients 

#### some simple evaluation ####
phi1_true = data$phi1_all
TPs = sum((phi1_all!=0)*(phi1_true!=0))   # TPs=19
FPs = sum((phi1_all!=0)*(phi1_true==0))   # FPs=5
SSE = sum((phi1_all-phi1_true)^2)         # SSE=0.968  


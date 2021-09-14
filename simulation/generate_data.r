N=500 # number of patients
N_test = 300
N_train = N-N_test
NY=3 # dimension of y
NZ=3 # number of covariates
NXI=1 # number of latent time-invariant variables
NU=6 # number of random effects 
NX = 0 # number of covariates in trajectory 
NETA=3 # number of latent time-variant variables 
latent_label = c(NY)
J_len = 6

t_test = c(0.5,0.75,1)
delta_test = c(0.25,0.5)
result_pre = matrix(0,REP,indi1_length*indi2_length)
result_AUC = matrix(0,REP,length(t_test)*length(delta_test))
Fail_label_tdelta = matrix(0,2,length(t_test)*length(delta_test))
################################################### True value #########################################################

phi = diag(1) # covarance matrix of XI
psi_eps = diag(rep(0.3,NY+1))[1:NY,1:NY] # variance of psi_epsilon
beta = c(1,-1,1) # coefficient of covariates
gamma = c(0) # coefficient of explanatory latent variable
c = 1 # coefficient of f(t)*u_i
b = c(1,0.6,0.7)[1:NY] # factor loading
b_label =c(0,1,1)[1:NY] # label the unfixed factor loading 
G = diag(NU) # covariance matrix of random effects
alpha = 0.5*c(-1,1,1,-1,1,-1)

cut = 0.02


K = 10 # number of grids used in quadrature method

c_mcmc = rep(0,S)
b_mcmc = array(0,dim=c(S,NY))
psi_epsmcmc = array(0,dim=c(S,NY))
G_mcmc = array(0,dim=c(S,NU,NU))
beta_mcmc = array(0,dim=c(S,NZ))
gamma_mcmc = array(0,dim=c(S,NXI))
lambda_mcmc = array(0,dim=c(S,J_len-1))
phi_mcmc = array(0,dim=c(S,NXI,NXI))
alpha_mcmc = array(0,dim=c(S,NU*(NX+1)))
u_mcmc = array(0,dim=c(S,N_train,NU))

c_rep = rep(0,REP)
b_rep = array(0,dim=c(REP,NY))
psi_epsrep = array(0,dim=c(REP,NY))
G_rep = array(0,dim=c(REP,NU,NU))
beta_rep = array(0,dim=c(REP,NZ))
gamma_rep = array(0,dim=c(REP,NXI))
lambda_rep = array(0,dim=c(REP,J_len-1))
phi_rep = array(0,dim=c(REP,NXI,NXI))
alpha_rep = array(0,dim=c(REP,NU*(NX+1)))
u_rep = array(0,dim=c(REP,N_train,NU))



c_xi = 5
c_u = 0.8
c_beta = 2
c_gamma = 18
c_c = 0.18
c_alpha = 0.6
c_u_test = 0.8
c_u_last = 0.8
mh = list()
mh[[1]] = c_u
mh[[2]] = c_beta
mh[[3]] = c_c
mh[[4]] = c_alpha
mh[[5]] = c_u_test
mh[[6]] = c_u_last






#######################################################################################################################
# generate data

observed_data = list()
test_data = list()

for(iter in 1:REP){
  kns = c(0.02,0.24,0.48,0.70,0.94,1.18)
  
  
  NT = numeric(N)
  
  FT = numeric(N) # failure time
  CT = numeric(N) # censoring time 
  OT = numeric(N) # observed time
  DELTA = numeric(N) # failure indicator
  XI = matrix(NA, nrow = N, ncol =NXI)
  u = matrix(NA, nrow = N, ncol = NU) # random effects
  Z = matrix(NA, nrow = N, ncol = NZ) # covariates
  X = matrix(NA, nrow = N, ncol = NX+1) # covariates in trajectory
  X[,NX+1] = 1
  d = rmvnorm(N,mean=rep(0,2))
  
  
  ETA = c() # response latent variable
  ft = c()
  Time = c() # measurement time
  Queue = c(0)
  Y = c()
  ############  Generate Data  #############
  for(i in 1:N){
    
    XI[i,] = rmvnorm(1,rep(0,NXI),phi)
    
    u[i,] = rmvnorm(1,rep(0,NU),G)

    Z[i,1] = rexp(1,1)-1
    Z[i,2] = rnorm(1,0,1)
    Z[i,3] = rt(1,5)
    
    if(NX>0){
      X[i,1:NX] = rnorm(1,0,1)
    }
    
    temp_ratio = -log(1-runif(1))/exp(Z[i,]%*%beta+XI[i,]%*%gamma)[1,1]
    rand = seq(0,3.5,0.001)
    
    f_rand = u[i,1]
    if(NU>1){
      f_rand = f_rand + u[i,2]*rand
    }
    if(NU>2){
      if(IF_QUADRATIC==1){
        f_rand = f_rand + u[i,3]*rand^2
      }else{
        for(j in 3:NU){
          f_rand = f_rand + u[i,j] *(-1)* ((abs(rand-kns[NU])^3+(rand-kns[NU])^3-(abs(rand-kns[j-2])^3+(rand-kns[j-2])^3))/(2*(kns[NU]-kns[j-2]))-(abs(rand-kns[NU])^3+(rand-kns[NU])^3-(abs(rand-kns[NU-1])^3+(rand-kns[NU-1])^3))/(2*(kns[NU]-kns[NU-1])))
        }
      }
    }
    
    integration = cumsum(exp(c*f_rand)*0.001) # baseline harzard function lambda_0 = 1
    
    FT[i] = which(temp_ratio<integration)[1]*0.001
    
    CT[i] = runif(1,1.1,2) # to make sure the censoring rate is about 30% and each patients had at least 3 observations.
    if(is.na(FT[i])){
      DELTA[i] = 0
      OT[i] = CT[i]
    }else if(FT[i]<CT[i]){
      DELTA[i] = 1
      OT[i] = FT[i]
    }else{
      DELTA[i] = 0
      OT[i] = CT[i]
    }# failure indicator
    
    
      if(DELTA[i]==1){
          NT[i] = min(trunc(CT[i]/cut),sample(6:9,1))
          Time_temp = sort(c(OT[i],sample(setdiff(seq(cut,trunc(CT[i]/cut)*cut,cut),OT[i]),NT[i]-2)))
          Time = c(Time,0,Time_temp)
        
      }else{
          NT[i] = sample(6:9,1)
          Time_temp = c(sort(sample(seq(cut,OT[i],cut),NT[i]-2)),OT[i])
            Time = c(Time,0,Time_temp)
      }

    
    Queue=c(Queue,sum(NT))  # No. of the first measurement of each patient (used in Rcpp)
    
    Y_temp = matrix(NA,nrow = NT[i], ncol = NY)
    
    for(j in 1:NT[i]){
      
      ETA_temp = u[i,1]
      if(NU>1){
        ETA_temp = ETA_temp + u[i,2]*Time[Queue[i]+j]
      }
      if(NU>2){
        if(IF_QUADRATIC==1){
          ETA_temp = ETA_temp + u[i,3]*Time[Queue[i]+j]^2
        }else{
          for(k in 3:NU){
            ETA_temp = ETA_temp + u[i,k] *(-1)* ((abs(Time[Queue[i]+j]-kns[NU])^3+(Time[Queue[i]+j]-kns[NU])^3-(abs(Time[Queue[i]+j]-kns[k-2])^3+(Time[Queue[i]+j]-kns[k-2])^3))/(2*(kns[NU]-kns[k-2]))-(abs(Time[Queue[i]+j]-kns[NU])^3+(Time[Queue[i]+j]-kns[NU])^3-(abs(Time[Queue[i]+j]-kns[NU-1])^3+(Time[Queue[i]+j]-kns[NU-1])^3))/(2*(kns[NU]-kns[NU-1])))
          }
        }
      }
      
      
      ft= c(ft,ETA_temp)
      
      ETA = c(ETA,ETA_temp)
      eps = as.vector(rmvnorm(1,rep(0,NY),as.matrix(psi_eps)))
      for(k in 1:latent_label[1]){
        Y_temp[j,k] = b[k]*ETA[Queue[i]+j] + eps[k]
      }
    }
    
    Y = rbind(Y,Y_temp)
    
  }
  
  cr = 1-sum(DELTA)/N
  SNT = sum(NT)
  Queue=Queue[-(N+1)]  # No. of the first measurement of each patient (used in Rcpp)
  
  # Basis spline for each subject i in time j
  B = array(0,dim=c(SNT,NU))
  B[,1] = 1
  if(NU>1){
    B[,2] = Time
  }
  if(NU>2){
    for(k in 3:NU){
      B[,k] = (-1)*((abs(Time-kns[NU])^3+(Time-kns[NU])^3-(abs(Time-kns[k-2])^3+(Time-kns[k-2])^3))/(2*(kns[NU]-kns[k-2]))-(abs(Time-kns[NU])^3+(Time-kns[NU])^3-(abs(Time-kns[NU-1])^3+(Time-kns[NU-1])^3))/(2*(kns[NU]-kns[NU-1])))
    }
    if(IF_QUADRATIC==1){
      B[,3] = Time^2
    }
  }

  B_train_cpp = as.vector(B[1:Queue[N_train+1],])
  
  OT_label = rep(NA,N)
  for(k in 1:N){
    OT_label[k] = which(Time[(Queue[k]+1):(Queue[k]+NT[k])]==OT[k])-1
  }
  
  # intervals of baseline harzard function
  J_cpp = c(0,as.vector(quantile(OT,seq(0,1,0.2)))[-c(1,6)],as.vector(quantile(OT,seq(0,1,0.2)))[6]+0.001)
  Fail_label = numeric(N)
  for(i in 1:N){
    Fail_label[i] = which(!(OT[i]>J_cpp))[1]-1
  }
      
  for(i in 1:length(t_test)){
    for(j in 1:length(delta_test)){
      Fail_label_tdelta[1,(i-1)*length(delta_test)+j] = which(!(t_test[i]>J_cpp))[1]-1
      Fail_label_tdelta[2,(i-1)*length(delta_test)+j] = which(!((t_test[i]+delta_test[j])>J_cpp))[1]-1
    }
  }


  Fail_label_last = matrix(0,indi1_length,N)
  for(i in 1:N){
    for(j in (indi1_start+1):indi1_end){
      if(NT[i]>j){
        Fail_label_last[j-indi1_start,i] = which(!(Time[Queue[i]+j]>J_cpp))[1]-1
      }
    }
  }
  
  
  y_train_cpp = as.vector(Y[1:Queue[N_train+1],])
  z_train_cpp = as.vector(Z[1:N_train,])
  x_train_cpp = as.vector(X[1:N_train,])
  
  
  y_test_cpp = Y[-(1:Queue[N_train+1]),]
  z_test_cpp = Z[-(1:N_train),]
  x_test_cpp = X[-(1:N_train),]
  
  u_label = rep(0,N)
  if(NU<3){
    u_label = rep(NU,N)
  }else{
    if(IF_QUADRATIC==1){
      u_label = rep(NU,N)
    }else{
      for(i in 1:N){
        u_label[i] = sum(Time[Queue[i]+NT[i]]>kns[1:(NU-2)])+2
      }
    }
  }
  
  u_num = c(cumsum(as.numeric(table(u_label[1:N_train]))[length(table(u_label[1:N_train])):1]),rep(N_train,NU-length(table(u_label[1:N_train]))))[NU:1]
  
  
  u_label_test = rep(0,length(t_test))
  if(NU<3){
    u_label_test = rep(NU,length(t_test))
  }else{
    if(IF_QUADRATIC==1){
      u_label_test = rep(NU,length(t_test))
    }else{
      for(j in 1:length(t_test)){
        u_label_test[j] = sum(t_test[j]>kns[1:(NU-2)])+2
      }
    }
  }
  
  

  u_label_last = matrix(0,indi1_length,N)
  if(NU<3){
    u_label_last = matrix(NU,indi1_length,N)
  }else{
    if(IF_QUADRATIC==1){
      u_label_last = matrix(NU,indi1_length,N)
    }else{
      for(i in 1:N){
        for(j in (indi1_start+1):indi1_end){
          if(NT[i]>j){
            u_label_last[j-indi1_start,i] = sum(Time[Queue[i]+j]>kns[1:(NU-2)])+2
          }
          
        }
        
      }
    }
  }
  
  
  
  num_list=17
  
  observed_data[[(iter-1)*num_list+1]] = x_train_cpp
  observed_data[[(iter-1)*num_list+2]] = OT_label[1:N_train]
  observed_data[[(iter-1)*num_list+3]] = B_train_cpp
  observed_data[[(iter-1)*num_list+4]] = kns
  observed_data[[(iter-1)*num_list+5]] = y_train_cpp
  observed_data[[(iter-1)*num_list+6]] = b_label
  observed_data[[(iter-1)*num_list+7]] = J_cpp
  observed_data[[(iter-1)*num_list+8]] = Fail_label[1:N_train]
  observed_data[[(iter-1)*num_list+9]] = z_train_cpp
  observed_data[[(iter-1)*num_list+10]] = OT[1:N_train]
  observed_data[[(iter-1)*num_list+11]] = DELTA[1:N_train]
  observed_data[[(iter-1)*num_list+12]] = NT[1:N_train]
  observed_data[[(iter-1)*num_list+13]] = Queue[1:N_train]
  observed_data[[(iter-1)*num_list+14]] = latent_label
  observed_data[[(iter-1)*num_list+15]] = u_label[1:N_train]
  observed_data[[(iter-1)*num_list+16]] = u_num
  observed_data[[(iter-1)*num_list+17]] = d
  
  num_list_test = 24
  
  test_data[[(iter-1)*num_list_test+1]] = t(x_test_cpp)
  test_data[[(iter-1)*num_list_test+2]] = OT_label[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+3]] = t(B[-(1:Queue[N_train+1]),])
  test_data[[(iter-1)*num_list_test+4]] = Fail_label_tdelta
  test_data[[(iter-1)*num_list_test+5]] = t(y_test_cpp)
  test_data[[(iter-1)*num_list_test+6]] = t_test
  test_data[[(iter-1)*num_list_test+7]] = delta_test
  test_data[[(iter-1)*num_list_test+8]] = Fail_label[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+9]] = t(z_test_cpp)
  test_data[[(iter-1)*num_list_test+10]] = OT[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+11]] = DELTA[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+12]] = NT[-(1:N_train)]
  test_data[[(iter-1)*num_list_test+13]] = Queue[-(1:N_train)]-Queue[-(1:N_train)][1]
  test_data[[(iter-1)*num_list_test+14]] = Time[-(1:Queue[N_train+1])]
  test_data[[(iter-1)*num_list_test+15]] = u_label_test
  test_data[[(iter-1)*num_list_test+16]] = matrix(0,1,N_test)
  test_data[[(iter-1)*num_list_test+17]] = matrix(0,length(t_test)*length(delta_test)*NU,N_test)
  test_data[[(iter-1)*num_list_test+18]] = matrix(u_label_last[,-(1:N_train)],indi1_length,N_test)
  test_data[[(iter-1)*num_list_test+19]] = matrix(0,NU*indi1_length,N_test)
  test_data[[(iter-1)*num_list_test+20]] = matrix(Fail_label_last[,-(1:N_train)],indi1_length,N_test)
  test_data[[(iter-1)*num_list_test+21]] = indi1_start
  test_data[[(iter-1)*num_list_test+22]] = indi1_end
  test_data[[(iter-1)*num_list_test+23]] = indi2_start
  test_data[[(iter-1)*num_list_test+24]] = indi2_end
  
  
}
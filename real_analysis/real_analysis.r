set.seed(9397)
library("mvtnorm")
library("numDeriv")
library("Rcpp")
library("RcppArmadillo")
library("RcppEigen")
S=30000
T=60000
NU = 6
dint = 10
iter=1
REP=1

indi1_start = 2
indi1_end = 3
indi2_start = 1
indi2_end = 4
indi1_length = indi1_end-indi1_start
indi2_length = indi2_end-indi2_start

N_test = 200
indicator_test = matrix(0,REP,N_test)
for(i in 1:REP){
  indicator_test[i,] = sample(1:715,N_test)
}
result_pre = matrix(0,REP,indi1_length*indi2_length)
sourceCpp("mcmc.cpp")
data=read.csv("data.csv")
DELTA=read.table("DELTA.txt")[[1]]
OT=read.table("OT.txt")[[1]]
NT=read.table("NT.txt")[[1]]
t_test = c(0.12,0.18)
delta_test = c(0.6,0.12,0.24)
result_AUC = matrix(0,REP,length(t_test)*length(delta_test))
IF_PRE = 0
IF_LAST = 0
IF_QUADRATIC = 0
IF_SEM = 1


for(i in 7:9){
  data[,i] = (data[,i]-mean(data[,i]))/sd(data[,i])
}


N = length(OT)
N_train = N-N_test


for(iter in 1:REP){
  
  label_test = indicator_test[iter,]
  label_train = c()
  if((IF_PRE == 1)|(IF_LAST==1)){
    label_train = (1:N)[-label_test]
  }else{
    label_train = 1:N
  }
  
  J_cpp = c(0,0.06,0.12,0.18,0.24,0.36,0.48,0.6,0.72,0.84,1.321)
  J_len=length(J_cpp)
  Fail_label = numeric(N)
  for(i in 1:N){
    Fail_label[i] = which(!(OT[i]>J_cpp))[1]-1
  }
  
  
  NY = 3
  NXI = 1
  NZ = 5
  NX = 2
  latent_label = c(NY)
  b_label=c(0,1,1)[1:NY]
  SNT = nrow(data)
  Queue = c(0,cumsum(NT))[-(N+1)]
  
  XI = matrix(NA, nrow = N, ncol =NXI)
  Z = matrix(NA, nrow = N, ncol = NZ) # covariates
  X = matrix(NA, nrow = N, ncol = (NX+1))
  
  for(i in 1:N){
    X[i,1] = (data$data_1.APOE4[Queue[i]+1]==1)+0
    X[i,2] = (data$data_1.APOE4[Queue[i]+1]==2)+0
    Z[i,1] = (data$data_1.PTMARRY[Queue[i]+1]=="Married")+0
    Z[i,2] = 1-(data$data_1.PTGENDER[Queue[i]+1]=="Male")
    Z[i,3] = data$data_1.PTEDUCAT[Queue[i]+1]
    Z[i,4] = (data$data_1.APOE4[Queue[i]+1]==1)+0
    Z[i,5] = (data$data_1.APOE4[Queue[i]+1]==2)+0
  }
  X[,NX+1] = 1
  Z[,3] = (Z[,3]-mean(Z[,3]))/sd(Z[,3])
  
  Y = matrix(cbind(data[,7],data[,8],data[,9])[,1:NY],ncol=NY)
  
  Time = (data$data_1.Month)/100
  
  Time_test = c()
  for(i in 1:length(label_test)){
    Time_test = c(Time_test,Time[Queue[label_test[i]]+1:NT[label_test[i]]])
  }
  
  
  OT_label = rep(NA,N)
  for(k in 1:N){
    OT_label[k] = which(round(Time[(Queue[k]+1):(Queue[k]+NT[k])],4)==round(OT[k],4))-1
  }
  
  
  Time_end = rep(0,N)
  for(i in 1:N){
    Time_end[i] = Time[Queue[i]+NT[i]]
  }
  kns = c(0.06,0.12,0.18,0.24,0.36,0.48,0.6,0.72,0.84)[1:NU]
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
  
  u_num = c(cumsum(as.numeric(table(u_label))[length(table(u_label)):1]),rep(N,NU-length(table(u_label))))[NU:1]
  
  u_num_train = c(cumsum(as.numeric(table(u_label[label_train]))[length(table(u_label[label_train])):1]),rep(length(label_train),NU-length(table(u_label[label_train]))))[NU:1]

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
  
  B_train = c()
  for(i in 1:length(label_train)){
    B_train = rbind(B_train, as.matrix(B[Queue[label_train[i]]+1:NT[label_train[i]],]) )
  }
  B_test = c()
  for(i in 1:length(label_test)){
    B_test = rbind(B_test, as.matrix(B[Queue[label_test[i]]+1:NT[label_test[i]],]))
  }
  
  B_train_cpp = as.vector(B_train)
  
  
  Fail_label_tdelta = matrix(0,2,length(t_test)*length(delta_test))
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
  
  Y_train = c()
  for(i in 1:length(label_train)){
    Y_train = rbind(Y_train,matrix(Y[Queue[label_train[i]]+1:NT[label_train[i]],],ncol=NY))
  }
  Y_test = c()
  for(i in 1:length(label_test)){
    Y_test = rbind(Y_test,matrix(Y[Queue[label_test[i]]+1:NT[label_test[i]],],ncol=NY))
  }
  
  y_train_cpp = as.vector(Y_train)
  z_train_cpp = as.vector(matrix(Z[label_train,],ncol=NZ))
  x_train_cpp = as.vector(matrix(X[label_train,],ncol=NX+1))
  
  
  y_test_cpp = Y_test
  z_test_cpp = matrix(Z[label_test,],ncol=NZ)
  x_test_cpp = matrix(X[label_test,],ncol=NX+1)
  
  
  
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
  
  
  observed_data = list()
  
  observed_data[[1]] = x_train_cpp
  observed_data[[2]] = OT_label[label_train]
  observed_data[[3]] = B_train_cpp
  observed_data[[4]] = kns
  observed_data[[5]] = y_train_cpp
  observed_data[[6]] = b_label
  observed_data[[7]] = J_cpp
  observed_data[[8]] = Fail_label[label_train]
  observed_data[[9]] = z_train_cpp
  observed_data[[10]] = OT[label_train]
  observed_data[[11]] = DELTA[label_train]
  observed_data[[12]] = NT[label_train]
  observed_data[[13]] = c(0,cumsum(NT[label_train]))[-(length(label_train)+1)]
  observed_data[[14]] = latent_label
  observed_data[[15]] = u_label[label_train]
  observed_data[[16]] = u_num_train
  
  test_data = list()
  
  test_data[[1]] = t(x_test_cpp)
  test_data[[2]] = OT_label[label_test]
  test_data[[3]] = t(B_test)
  test_data[[4]] = Fail_label_tdelta
  test_data[[5]] = t(y_test_cpp)
  test_data[[6]] = t_test
  test_data[[7]] = delta_test
  test_data[[8]] = Fail_label[label_test]
  test_data[[9]] = t(z_test_cpp)
  test_data[[10]] = OT[label_test]
  test_data[[11]] = DELTA[label_test]
  test_data[[12]] = NT[label_test]
  test_data[[13]] = c(0,cumsum(NT[label_test]))[-(length(label_test)+1)]
  test_data[[14]] = Time_test
  test_data[[15]] = u_label_test
  test_data[[16]] = matrix(0,1,N_test)
  test_data[[17]] = matrix(0,length(t_test)*length(delta_test)*NU,N_test)
  test_data[[18]] = matrix(u_label_last[,label_test],indi1_length,N_test)
  test_data[[19]] = matrix(0,NU*indi1_length,N_test)
  test_data[[20]] = matrix(Fail_label_last[,label_test],indi1_length,N_test)
  test_data[[21]] = indi1_start
  test_data[[22]] = indi1_end
  test_data[[23]] = indi2_start
  test_data[[24]] = indi2_end
  
  
  
  
  c_mcmc = rep(0,S)
  b_mcmc = array(0,dim=c(S,NY))
  psi_epsmcmc = array(0,dim=c(S,NY))
  G_mcmc = array(0,dim=c(S,NU,NU))
  beta_mcmc = array(0,dim=c(S,NZ))
  lambda_mcmc = array(0,dim=c(S,J_len-1))
  alpha_mcmc = array(0,dim=c(S,NU*(NX+1)))
  u_mcmc = array(0,dim=c(S,length(label_train),NU))
  
  c_rep = rep(0,REP)
  b_rep = array(0,dim=c(REP,NY))
  psi_epsrep = array(0,dim=c(REP,NY))
  G_rep = array(0,dim=c(REP,NU,NU))
  beta_rep = array(0,dim=c(REP,NZ))
  lambda_rep = array(0,dim=c(REP,J_len-1))
  alpha_rep = array(0,dim=c(REP,NU*(NX+1)))
  u_rep = array(0,dim=c(REP,length(label_train),NU))
  
  c_sd = rep(0,REP)
  b_sd = array(0,dim=c(REP,NY))
  psi_epssd = array(0,dim=c(REP,NY))
  G_sd = array(0,dim=c(REP,NU,NU))
  beta_sd = array(0,dim=c(REP,NZ))
  lambda_sd = array(0,dim=c(REP,J_len-1))
  alpha_sd = array(0,dim=c(REP,NU*(NX+1)))
  u_sd = array(0,dim=c(REP,length(label_train),NU))
  
  ################################################   Real Data Analysis   ##############################################
  
  #for(iter in 1:REP){
  num_list = 15
  K = 10
  ####### initial value #######
  G_cpp = diag(NU)
  u_cpp = rep(0,length(label_train)*NU)
  b_cpp = rep(1,NY)
  psi_epscpp = rep(1,NY)
  xi_cpp = rep(0,length(label_train)*NXI)
  phi_cpp = diag(1)
  lambda_cpp = rep(1,length(J_cpp)-1)
  alpha_cpp = rep(0,NU*(NX+1)) # coefficient in second level model
  c_cpp = 0
  gamma_cpp = c(0)
  beta_cpp = rep(0,NZ)
  
  ####### prior distribution ########
  
  
  
  # G
  rho0 = 7
  G0 = 4*diag(NU)
  # b
  b0 = 0
  h0 = 1
  # psi_eps
  alpha_eps0 = 9
  beta_eps0 = 4
  # lambda
  alpha1 = 2
  alpha2 = 0.01
  # beta
  beta0 = rep(0,NZ)
  Sigma_beta0 = diag(NZ)
  # c
  c0 = 0
  sigma_c0 = 1
  # alpha
  alpha0 = rep(0,NU*(NX+1))
  Sigma_alpha0 = diag(NU*(NX+1))
  
  # # G
  # rho0 = 4
  # G0 = 2*diag(NU)
  # # b
  # b0 = 0
  # h0 = 10^4
  # # psi_eps
  # alpha_eps0 = 3
  # beta_eps0 = 2
  # # lambda
  # alpha1 = 2
  # alpha2 = 0.01
  # # beta
  # beta0 = rep(0,NZ)
  # Sigma_beta0 = 10^4*diag(NZ)
  # # c
  # c0 = 0
  # sigma_c0 = 10^4
  # # alpha
  # alpha0 = rep(0,NU*(NX+1))
  # Sigma_alpha0 = 10^4*diag(NU*(NX+1))
  
  prior = list()
  prior[[1]] = rho0
  prior[[2]] = G0
  prior[[3]] = b0
  prior[[4]] = h0
  prior[[5]] = alpha_eps0
  prior[[6]] = beta_eps0
  prior[[7]] = alpha1
  prior[[8]] = alpha2
  prior[[9]] = beta0
  prior[[10]] = Sigma_beta0
  prior[[11]] = c0
  prior[[12]] = sigma_c0
  prior[[13]] = alpha0
  prior[[14]] = Sigma_alpha0
  
  
  c_u = 1
  c_beta = 2
  c_c = 0.2
  c_alpha = 0.5
  c_u_test = 0.4
  c_u_last = 0.4
  mh = list()
  mh[[1]] = c_u
  mh[[2]] = c_beta
  mh[[3]] = c_c
  mh[[4]] = c_alpha
  mh[[5]] = c_u_test
  mh[[6]] = c_u_last
  
  
  
  para = list()
  para[[1]] = u_cpp
  para[[2]] = G_cpp
  para[[3]] = b_cpp
  para[[4]] = psi_epscpp
  para[[5]] = lambda_cpp
  para[[6]] = beta_cpp
  para[[7]] = alpha_cpp
  para[[8]] = c_cpp
  para[[9]] = xi_cpp
  para[[10]] = gamma_cpp
  para[[11]] = phi_cpp
  
  result = mcmc(observed_data, test_data, para, mh, prior, K, iter, REP, T, S, IF_PRE, IF_LAST, IF_QUADRATIC, IF_SEM)
  
  c_mcmc = result$c_str
  b_mcmc = matrix(result$b_str,nrow=S,byrow=T) 
  psi_epsmcmc = matrix(result$psieps_str,nrow=S,byrow=T) 
  for(i in 1:S){
    G_mcmc[i,,] = matrix(result$G_str[((i-1)*NU*NU+1):(i*NU*NU)],nrow=NU,byrow=T)
    u_mcmc[i,,] = matrix(result$u_str[((i-1)*length(label_train)*NU+1):(i*length(label_train)*NU)],nrow=length(label_train),byrow=F)
  }
  beta_mcmc = matrix(result$beta_str,nrow=S,byrow=T) 
  lambda_mcmc = matrix(result$lambda_str,nrow=S,byrow=T) 
  alpha_mcmc = matrix(result$alpha_str,nrow=S,byrow=T)
  
  
  
    c_rep[iter] = mean(c_mcmc)
    b_rep[iter,] = apply(b_mcmc,2,mean)
    psi_epsrep[iter,] = apply(psi_epsmcmc,2,mean)
    G_rep[iter,,] = apply(G_mcmc,2:3,mean)
    beta_rep[iter,] = apply(beta_mcmc,2,mean)
    lambda_rep[iter,] = apply(lambda_mcmc,2,mean)
    alpha_rep[iter,] = apply(alpha_mcmc,2,mean)
    u_rep[iter,,] = apply(u_mcmc,2:3,mean)
    
    
    c_sd[iter] = sd(c_mcmc)
    b_sd[iter,] = apply(b_mcmc,2,sd)
    psi_epssd[iter,] = apply(psi_epsmcmc,2,sd)
    G_sd[iter,,] = apply(G_mcmc,2:3,sd)
    beta_sd[iter,] = apply(beta_mcmc,2,sd)
    lambda_sd[iter,] = apply(lambda_mcmc,2,sd)
    alpha_sd[iter,] = apply(alpha_mcmc,2,sd)
    u_sd[iter,,] = apply(u_mcmc,2:3,sd)
    
    
    c_rep
    c_sd
    b_rep
    b_sd
    psi_epsrep
    psi_epssd
    G_rep
    G_sd
    beta_rep
    beta_sd
    lambda_rep
    lambda_sd
  
  write.table(beta_mcmc,"beta_mcmc.txt")
  write.table(alpha_mcmc,"alpha_mcmc.txt")
  write.table(c_mcmc,"c_mcmc.txt")
  write.table(b_mcmc,"b_mcmc.txt")
  write.table(psi_epsmcmc,"psi_epsmcmc.txt")
  write.table(lambda_mcmc,"lambda_mcmc.txt")
  write.table(G_mcmc,"G_mcmc.txt")
  write.table(u_rep[1,,],"u_rep1.txt")
  apply(result$AUC,2,mean)
  result_AUC[iter,] = apply(result$AUC,2,mean)
  
  if(IF_LAST==1){
      cog = read.table("cognitive_impairment.txt")[,1] # cognitive impairment estimated by a seperate CFA model
      for(indi1 in (indi1_start+1):indi1_end){
        for(indi2 in indi2_start:(indi2_end-1)){
          result_pre[iter,(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1] = cor(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],cog[Queue+indi1+indi2][label_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0])
          
        }
      }
      if(REP==1){
        for(indi1 in (indi1_start+1):indi1_end){
          for(indi2 in indi2_start:(indi2_end-1)){
            plot(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],cog[Queue+indi1+indi2][label_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0])
            abline(0,1)
            print(cor(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],cog[Queue+indi1+indi2][label_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0]))
            
          }
        }
      }
    
  }
  
  
  
  
}

if(IF_LAST==1){
  print(result_pre)
  write.table(result_pre,"result_pre.txt")
}
if(IF_PRE==1){
  print(result_AUC)
  write.table(result_AUC,"result_AUC.txt")
}
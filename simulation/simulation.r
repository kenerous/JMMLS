#set.seed(9397)
library("mvtnorm")
library("numDeriv")
library("Rcpp")
library("RcppArmadillo")
library("RcppEigen")

sourceCpp("mcmc.cpp")

####################################### Simulation setting and control ###################################################
REP=1      ## number of replication
T=40000    ## number of iteration of mcmc
S=20000    ## T-S is the number of burn-in iteration
IF_PRE = 1  ## Conduct out-of-sample prediction for survival probabilities
IF_LAST = 1  ## Conduct out-of-sample prediction for latent cognitive impairment
IF_QUADRATIC = 0
IF_SEM = 1

indi1_start = 2
indi1_end = 5
indi2_start = 1
indi2_end = 4
indi1_length = indi1_end-indi1_start
indi2_length = indi2_end-indi2_start

source("generate_data.r")

#######################################################################################################################
# MCMC
for(iter in 1:REP){
  num_list = 17
  num_list_test = 24
  ####### initial value #######
    G_cpp = diag(NU)
    u_cpp = rep(0,N_train*NU)
    b_cpp = rep(1,NY)
    psi_epscpp = rep(1,NY)
    xi_cpp = as.vector(rmvnorm(N_train,rep(0,NXI),diag(NXI)))
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
  
  
  result = mcmc(observed_data[((iter-1)*num_list+1):((iter-1)*num_list+16)], test_data[((iter-1)*num_list_test+1):((iter-1)*num_list_test+num_list_test)], para, mh, prior, K, iter, REP, T, S, IF_PRE, IF_LAST, IF_QUADRATIC, IF_SEM)
  result_AUC[iter,] = apply(result$AUC,2,mean)
  if(IF_LAST==1){
    for(indi1 in (indi1_start+1):indi1_end){
      for(indi2 in indi2_start:(indi2_end-1)){
        result_pre[iter,(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1] = cor(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],ETA[Queue+indi1+indi2][N_train+1:N_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0])
        #print(cor(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],cfa$scores[][Queue+indi1+indi2][N_train+1:N_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0]))
      }
    }
    if(REP==1){
      for(indi1 in (indi1_start+1):indi1_end){
        for(indi2 in indi2_start:(indi2_end-1)){
          plot(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],ETA[Queue+indi1+indi2][N_train+1:N_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0])
          abline(0,1)
          print(cor(result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0],ETA[Queue+indi1+indi2][N_train+1:N_test][result$eta_last[(indi1-indi1_start-1)*indi2_length+indi2-indi2_start+1,]!=0]))
          
        }
      }
    }
  }
  
  c_mcmc = result$c_str
  b_mcmc = matrix(result$b_str,nrow=S,byrow=T) 
  psi_epsmcmc = matrix(result$psieps_str,nrow=S,byrow=T) 
  for(i in 1:S){
    G_mcmc[i,,] = matrix(result$G_str[((i-1)*NU*NU+1):(i*NU*NU)],nrow=NU,byrow=T)
    u_mcmc[i,,] = matrix(result$u_str[((i-1)*N_train*NU+1):(i*N_train*NU)],nrow=N_train,byrow=F)
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
  
  print(mean(c_rep[1:iter]))
  if(iter>1){
    print(apply(beta_rep[1:iter,],2,mean))
  }else{
    print(beta_rep[1,])
  }
  
}
write.table(c_rep,"c_rep.txt")
write.table(b_rep,"b_rep.txt")
write.table(psi_epsrep,"psi_epsrep.txt")
write.table(G_rep,"G_rep.txt")
write.table(beta_rep,"beta_rep.txt")
write.table(alpha_rep,"alpha_rep.txt")
write.table(lambda_rep,"lambda_rep.txt")


mean(c_rep)
apply(b_rep,2,mean)
apply(psi_epsrep,2,mean)
apply(beta_rep,2,mean)
apply(alpha_rep,2,mean)
apply(lambda_rep,2,mean)

result$at_beta/T
result$at_c/T



c_rep
b_rep
psi_epsrep
beta_rep
lambda_rep
if(IF_LAST==1){
  print(result_pre)
  write.table(result_pre,"result_pre.txt")
}
if(IF_PRE==1){
  print(result_AUC)
  write.table(result_AUC,"result_AUC.txt")
}


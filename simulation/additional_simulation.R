pi_est <- result$pi_test


pi_true <- pi_est
u_true <- u[201:500,]

betaz <- (Z%*%beta)[,1]

# loop1 <- 0
# loop2 <- 2
delta_length <- length(delta_test)

for(loop1 in 0:2){
  for(loop2 in 1:2){
    t_temp1 <- t_test[loop1+1]
    delta_temp1 <- delta_test[loop2]
    
    for(i in 201:500){
      if(pi_est[loop1*delta_length+loop2,i-200]>0){
        pi_true[loop1*delta_length+loop2,i-200] <- exp( exp(betaz[i])*( (-1) *f_integral(IF_QUADRATIC, NX-1, c(1), alpha, kns, u[i,], J_cpp, rep(1,5), Fail_label_tdelta[2,loop1*delta_length+loop2], K, NU, t_temp1+delta_temp1, c, 0, Fail_label_tdelta[2,loop1*delta_length+loop2]) + f_integral(IF_QUADRATIC, NX-1, c(1), alpha, kns, u[i,], J_cpp, rep(1,5), Fail_label_tdelta[1,loop1*delta_length+loop2], K, NU, t_temp1, c, 0, Fail_label_tdelta[1,loop1*delta_length+loop2])   ) )  
      }
      
      
      
    }
    
    
    
    
  }
}


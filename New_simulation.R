PS_simulation_new <- function(n, p1, K, beta, alpha, delta, Homo, delta_function,seeds) {
  require(MASS)
  set.seed(seeds)
  X0 = rep(1, n)
  X4 = rbinom(n, size = 1, prob = p1)
  X4_p = 0.6*X4 + 0.4*(1-X4)
  X3 = c()
  for (i in 1:length(X4_p)) {
    X3[i] = rbinom(1, size = 1, prob = X4_p[i])
  }
  mean_X1 = -X3 + X4 + 0.5*X3*X4
  mean_X2 = X3 - X4 + X3*X4
  M1 = matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
  M2 = matrix(c(2,0.25,0.25,2), nrow = 2, ncol = 2)
  X1_X2 = matrix(nrow = n, ncol = 2, byrow = T)
  for (i in 1:n) {
    X1_X2[i,] = mvrnorm(1, mu = c(mean_X1[i], mean_X2[i]), Sigma = X3[i]*M1 + (1-X3[i])*M2)
  }
  if(K == 0) {
    X_1to4 = cbind(X1_X2, X3, X4)
    Data_X = cbind(X0, X_1to4, X_1to4[,1]*X_1to4[,2], 
                   X_1to4[,1]*X_1to4[,3], X_1to4[,1]*X_1to4[,4],
                   X_1to4[,2]*X_1to4[,3], X_1to4[,2]*X_1to4[,4],
                   X_1to4[,3]*X_1to4[,4], X_1to4^2, log(abs(X_1to4[,1:2])))
  } else if (K != 0) {
    X_normal = matrix(nrow = n, ncol = round(K/2), byrow = T)
    for(i in 1:n) {
      X_normal[i,] = mvrnorm(1, mu = rep(0, round(K/2)), 
                             Sigma = diag(x = 1.6, round(K/2))+matrix(rep(0.4,round(K/2)*round(K/2)), round(K/2)))
    }
    X_binary = matrix(nrow = n, ncol = round(K/2), byrow = T)
    for (i in 1:round(K/2)) {
      X_binary[,i] = rbinom(n, size = 1, prob = p1)
    }
    X_1to4 = cbind(X1_X2, X3, X4)
    Data_X = cbind(X0, X_1to4, X_1to4[,1]*X_1to4[,2], 
                   X_1to4[,1]*X_1to4[,3], X_1to4[,1]*X_1to4[,4],
                   X_1to4[,2]*X_1to4[,3], X_1to4[,2]*X_1to4[,4],
                   X_1to4[,3]*X_1to4[,4], X_1to4^2, log(abs(X_1to4[,1:2])), X_normal, X_binary)
  }
  e_0 = exp(Data_X%*%beta)
  PS = c()
  for (i in 1:n) {
    PS[i] = (e_0[i])/(1+e_0[i])
  }
  Z = c()
  for (i in 1:n) {
    Z[i] = rbinom(1, size = 1, prob = PS[i])
  }
  epsilon = rnorm(n)
  if ( Homo == 0) {
    delta = rep(1*delta, n)
  } else if (Homo == 1) {
    delta = delta_function(PS)
  } 
  Y = delta*Z + Data_X%*%alpha + epsilon
  Data_final = cbind(Y, Z, Data_X)
  if (K==0) {
    colnames(Data_final) = c("Y", "Z", "X0", "X1", "X2","X3","X4", 
                             "X1*X2", "X1*X3", "X1*X4", "X2*X3", 
                             "X2*X4", "X3*X4", "X1^2", "X2^2", "X3^2",
                             "X4^2","logX1","logX2")
  } else if (K!=0) {
    colnames(Data_final) = c("Y", "Z", "X0", "X1", "X2","X3","X4", 
                             "X1*X2", "X1*X3", "X1*X4", "X2*X3", 
                             "X2*X4", "X3*X4", "X1^2", "X2^2", "X3^2",
                             "X4^2", "logX1","logX2", paste0("Norm", seq(1, round(K/2))),
                             paste0("Bin", seq(1, round(K/2))))
  }
  return(list(Data=Data_final, propensity_score = PS))
}

delta_function_power <- function(e) {
  D = -4*(e^2)+3.94*e+0.6864
  return(D)
}

delta_function_1 <- function(e) {
  D = 4*(e^2)+3.94*e+0.6864
  return(D)
}

delta_function_2 <- function(e) {
  D = 4*(e^2)+0.6864
  return(D)
}

delta_function_linear <- function(e) {
  D = 3.94*e+0.6864
  return(D)
}

Hetero_f1 = function(e) {
  D = e + 1
  return(D)
}

Hetero_f2 = function(e) {
  D = e^2 + 2*e + 1
  return(D)
}

PS_simulation_OneCov <- function(n, p1, beta, alpha, delta, Homo, delta_function,seeds) {
  require(MASS)
  set.seed(seeds)
  X0 = rep(1, n)
  X1 = rnorm(n, mean = 0, sd = 2)
  X1_2 = X1^2
  logX1 = log(abs(X1))
  Data_X = cbind(X0,X1,X1_2,logX1)
  e_0 = exp(Data_X%*%beta)
  PS = c()
  for (i in 1:n) {
    PS[i] = (e_0[i])/(1+e_0[i])
  }
  Z = c()
  for (i in 1:n) {
    Z[i] = rbinom(1, size = 1, prob = PS[i])
  }
  epsilon = rnorm(n)
  if ( Homo == 0) {
    delta = rep(1*delta, n)
  } else if (Homo == 1) {
    delta = delta_function(PS)
  } 
  Y = delta*Z + Data_X%*%alpha + epsilon
  Data_final = cbind(Y, Z, Data_X)
  colnames(Data_final) = c("Y", "Z", "X0", "X1", "X1^2","logX1")
  return(list(Data=Data_final, propensity_score = PS))
}

#DM_12 = PS_simulation_new(1000, 0.5, 8, 
 #    beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,1,1,rep(0,8),beta_12),nrow = 23),
  #   alpha = matrix(c(0.5,1,0.6,2.2,-1.2,rep(0,10),alpha_12),nrow = 23),
   #  Homo = 0, delta_function = delta_function_1, 1234)

Find_true_delta_new = function(n, p1, K, beta,delta_function) {
  require(MASS)
  X0 = rep(1, n)
  X4 = rbinom(n, size = 1, prob = p1)
  X4_p = 0.6*X4 + 0.4*(1-X4)
  X3 = c()
  for (i in 1:length(X4_p)) {
    X3[i] = rbinom(1, size = 1, prob = X4_p[i])
  }
  mean_X1 = -X3 + X4 + 0.5*X3*X4
  mean_X2 = X3 - X4 + X3*X4
  M1 = matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
  M2 = matrix(c(2,0.25,0.25,2), nrow = 2, ncol = 2)
  X1_X2 = matrix(nrow = n, ncol = 2, byrow = T)
  for (i in 1:n) {
    X1_X2[i,] = mvrnorm(1, mu = c(mean_X1[i], mean_X2[i]), Sigma = X3[i]*M1 + (1-X3[i])*M2)
  }
  if(K == 0) {
    X_1to4 = cbind(X1_X2, X3, X4)
    Data_X = cbind(X0, X_1to4, X_1to4[,1]*X_1to4[,2], 
                   X_1to4[,1]*X_1to4[,3], X_1to4[,1]*X_1to4[,4],
                   X_1to4[,2]*X_1to4[,3], X_1to4[,2]*X_1to4[,4],
                   X_1to4[,3]*X_1to4[,4], X_1to4^2, log(abs(X_1to4[,1:2])))
  } else if (K != 0) {
    X_normal = matrix(nrow = n, ncol = round(K/2), byrow = T)
    for(i in 1:n) {
      X_normal[i,] = mvrnorm(1, mu = rep(0, round(K/2)), 
                             Sigma = diag(x = 1.6, round(K/2))+matrix(rep(0.4,round(K/2)*round(K/2)), round(K/2)))
    }
    X_binary = matrix(nrow = n, ncol = round(K/2), byrow = T)
    for (i in 1:round(K/2)) {
      X_binary[,i] = rbinom(n, size = 1, prob = p1)
    }
    X_1to4 = cbind(X1_X2, X3, X4)
    Data_X = cbind(X0, X_1to4, X_1to4[,1]*X_1to4[,2], 
                   X_1to4[,1]*X_1to4[,3], X_1to4[,1]*X_1to4[,4],
                   X_1to4[,2]*X_1to4[,3], X_1to4[,2]*X_1to4[,4],
                   X_1to4[,3]*X_1to4[,4], X_1to4^2, log(abs(X_1to4[,1:2])), X_normal, X_binary)
  }
  e_0 = exp(Data_X%*%beta)
  PS = c()
  for (i in 1:n) {
    PS[i] = (e_0[i])/(1+e_0[i])
  }
  Z = c()
  for (i in 1:n) {
    Z[i] = rbinom(1, size = 1, prob = PS[i])
  }
  delta = delta_function(PS)
  DM_good = cbind(Z, PS, delta, Data_X)
  W_good = Weight_generate(DM_good, PS, 1)
  TrueDelta_ipw = mean(W_good[[1]]*delta)/mean(W_good[[1]])
  return(True_delta_ipw = TrueDelta_ipw)
}

Weight_generate = function(DM, fitted, treat_I) {
  ipw = 1/(DM[,treat_I]*fitted + (1-DM[,treat_I])*(1-fitted))
  MW = c()
  for (i in 1:length(fitted)) {
    MW[i] = min(fitted[i], (1-fitted[i]))/(DM[i,treat_I]*fitted[i] + (1-DM[i,treat_I])*(1-fitted[i]))
  }
  OW = (fitted*(1-fitted))/(DM[,treat_I]*fitted + (1-DM[,treat_I])*(1-fitted))
  return(list(ipw, MW, OW))
}








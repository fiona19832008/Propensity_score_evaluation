#R_list = list.files("/Users/yan/Desktop/Propensity score/Report/Evaluation/R_code/Qingyuan/covalign/R"
 #                  , full.names = T)
#for(i in 1:length(R_list)) {
 # source(R_list[i])
#}

# Kernel transformation function
Kernel_transform = function(X, r) {
  X_mean = colMeans(X)
  X_sd = apply(X, 2, sd)
  X_std = matrix(nrow = nrow(X), ncol = ncol(X), byrow = T)
  for (i in 1:nrow(X)) {
    X_std[i, ] = (X[i, ] - X_mean)/X_sd
  }
  sigma = sigest(X_std)
  rbf = rbfdot(sigma = as.numeric(sigma[2]))
  K = kernelMatrix(rbf, X_std)
  K.eigen = eigen(K)
  features = K.eigen$vectors
  d = K.eigen$values
  r = min(r, ncol(features))
  features = features[, 1:r]
  d = 1/d[1:r]
  return(list(features, d))
}

Kernel_transform_choose_sigma = function(X, n_sigma, r) {
  X_mean = colMeans(X)
  X_sd = apply(X, 2, sd)
  X_std = matrix(nrow = nrow(X), ncol = ncol(X), byrow = T)
  for (i in 1:nrow(X)) {
    X_std[i, ] = (X[i, ] - X_mean)/X_sd
  }
  sigma = sigest(X_std)
  sigma_range = c(as.numeric(sigma[1]),  as.numeric(sigma[3]))
  log_sigma_range = log(sigma_range)
  sigma_seq = seq(from = log_sigma_range[1], to = log_sigma_range[2], length.out =  n_sigma)
  sigma_seq = exp(sigma_seq)
  features_out = list()
  for(i in 1:length(sigma_seq)) {
    rbf = rbfdot(sigma = sigma_seq[i])
    K = kernelMatrix(rbf, X_std)
    K.eigen = eigen(K)
    features = K.eigen$vectors
    d = K.eigen$values
    r = min(r, ncol(features))
    features = features[, 1:r]
    d = 1/d[1:r]
    features_out[[i]] = list(features, d, sigma_seq[i])
  }
  return(features_out)
}

Kernel_choose_ridge = function(Z, features, d, params, lambda.seq, normalize = normalize,
                               alpha, beta, std) {
  # std = c("std.diff", "std.norm")
  result <- list()
  for (setting in 1:nrow(params)) {
    result[[setting]] <- tryCatch(kernel.balance.path(Z, features, K = NULL, d = d,
                                             alpha = params[setting, 1],
                                             beta = params[setting, 2],
                                             lambda = lambda.seq,
                                             tol = 1e-10),
                                  error=function(err) vector("list",length(lambda.seq)))
  }
  Bern_weight = matrix(nrow = nrow(features), ncol = length(lambda.seq), byrow = F)
  Treat_effect_weight = Bern_weight
  for (i in 1:length(lambda.seq)) {
    if(sum(is.null(result[[1]][[i]]))==1) {
      Bern_weight[,i] = rep(NA, nrow(features))
    } else {
      Bern_weight[,i] = compute.weights(Z, result[[1]][[i]]$p,
                                        alpha = 0, beta = 0,
                                        normalize = normalize)
    }
    if(sum(is.null(result[[2]][[i]]))==1) {
      Treat_effect_weight[,i] = rep(NA, nrow(features))
    } else {
      Treat_effect_weight[,i]= compute.weights(Z, result[[2]][[i]]$p,
                                               alpha = alpha, beta = beta,
                                               normalize = normalize)
    }
  }
  imbalance <- matrix(0, length(result[[setting]]), ncol(features))
  #p.value <- imbalance
  imbalance_bern <- matrix(0, length(result[[setting]]), ncol(features))
  #p.value_bern <- imbalance
  for (i in 1:length(result[[setting]])) {
    # Calculate the weighted difference of X between case and control
    if(sum(is.na(Treat_effect_weight[,i])==nrow(Treat_effect_weight))) {
      imbalance[i,] = rep(NA, ncol(features))
    } else {
      imbalance[i,] <- ATE_Weight_diff(cbind(Z, features), Treat_effect_weight[,i], X_index = c(2:(ncol(features)+1)), Z_index = 1, std = std)[,3]
    }
    if(sum(is.na(Treat_effect_weight[,i])==nrow(Treat_effect_weight))) {
      imbalance_bern[i,] = rep(NA, ncol(features))
    } else {
      imbalance_bern[i,] <- ATE_Weight_diff(cbind(Z, features), Bern_weight[,i], X_index = c(2:(ncol(features)+1)), Z_index = 1, std = std)[,3]
    }
  }
  norm_imbalance = rowSums(imbalance^2)
  norm_imbalance_bern = rowSums(imbalance_bern^2)
  if(sum(is.na(norm_imbalance))==length(norm_imbalance)&
     sum(is.na(norm_imbalance_bern)) != length(norm_imbalance_bern)) {
    best = c(NA, which.min(norm_imbalance_bern))
    best_out = cbind(best, c(NA, lambda.seq[best[2]]))
    best_diff = c(min(norm_imbalance), min(norm_imbalance_bern))
  } else if(sum(is.na(norm_imbalance_bern))==length(norm_imbalance_bern)&
            sum(is.na(norm_imbalance)) != length(norm_imbalance)) {
    best = c(which.min(norm_imbalance), NA)
    best_out = cbind(best, c(lambda.seq[best[1]],NA))
    best_diff = c(min(norm_imbalance), min(norm_imbalance_bern))
  } else if(sum(is.na(norm_imbalance_bern)) == length(norm_imbalance_bern)&
            sum(is.na(norm_imbalance)) == length(norm_imbalance)) {
    best = c(NA,NA)
    best_out = cbind(best, c(NA,NA))
    best_diff = c(NA,NA)
  } else {
    best = c(which.min(norm_imbalance), 
             which.min(norm_imbalance_bern))
    best_out = cbind(best, lambda.seq[best])
    best_diff = c(min(norm_imbalance), min(norm_imbalance_bern))
  }
  if (sum(is.na(best[1]))==1&sum(is.na(best[2]))!=1) {
    weight = cbind(rep(0,length(Z)), 
                   compute.weights(Z, result[[1]][[best[2]]]$p,
                                   alpha = 0, beta = 0,
                                   normalize = normalize),
                   rep(0,length(Z)),  result[[1]][[best[2]]]$p)
    Coef = cbind(rep(0, length(result[[1]][[best[2]]]$eta)), result[[1]][[best[2]]]$eta)
    colnames(Coef) = c("treat_coef","bernoulli_coef")
  } else if (sum(is.na(best[2]))==1&sum(is.na(best[1])!=1)) {
    weight = cbind(compute.weights(Z, result[[2]][[best[1]]]$p,
                                   alpha = alpha, beta = beta,
                                   normalize = normalize),
                   rep(0,length(Z)), 
                   result[[2]][[best[1]]]$p,rep(0,length(Z)))
    Coef = cbind(result[[2]][[best[1]]]$eta,rep(0,length(result[[2]][[best[1]]]$eta)))
    colnames(Coef) = c("treat_coef","bernoulli_coef")
  } else if (sum(is.na(best))==2) {
    weight = cbind(rep(0,length(Z)), rep(0,length(Z)), rep(0,length(Z)), rep(0,length(Z)))
    Coef = NULL
  } else {
    weight = cbind(compute.weights(Z, result[[2]][[best[1]]]$p,
                                   alpha = alpha, beta = beta,
                                   normalize = normalize),
                   compute.weights(Z, result[[1]][[best[2]]]$p,
                                   alpha = 0, beta = 0,
                                   normalize = normalize), 
                   result[[2]][[best[1]]]$p,  result[[1]][[best[2]]]$p)
    Coef = cbind(result[[2]][[best[1]]]$eta, result[[1]][[best[1]]]$eta)
    colnames(Coef) = c("treat_coef", "bernoulli_coef")
  }
  colnames(weight) = c("treat_weight", "bernoulli_weight","treat_PS", "bernoulli_PS")
  return(list(best_diff = best_diff, weight = weight, Coef =  Coef, best_out = best_out))
}

Tailored_loss_kernel_weight = function(DM, X_index, Z_index, estimand, kernel_TF, r=100, 
                                lambda_test = c(0, -10, 22), normalize = FALSE, n_sigma, std = std) {
  # estimand: type of estimated outcome Y ("ATE", "ATT", "ATC")
  # kernel_TF: TRUE for using Gaussian kernel; FALSE for using the original X matrix
  # r: the first number of features
  # lambda_test: lambda.seq start value, end value, length [0, -10, 22] as default
  # kernel_TF: "kernel", "finite", "choose_sigma"
  #R_list = list.files(Rfile_path, full.names = T)
  #for(i in 1:length(R_list)) {
   # source(R_list[i])
  #}
  library(kernlab)
  library(Matching)
  
  Z = DM[,Z_index]
  X = as.matrix(DM[,X_index])
  alpha = switch (estimand, ATE = -1, ATT = 0, ATC = -1)
  beta = switch(estimand, ATE = -1, ATT = -1, ATC = 0)
  params = matrix(c(0, 0, alpha, beta), byrow = TRUE, ncol = 2)
  lambda.seq <- c(0, 10^seq(lambda_test[1], lambda_test[2], length = lambda_test[3]))
  
  if (kernel_TF == "kernel") {
    K = Kernel_transform(X, r)
    features = K[[1]]
    d = K[[2]]
    kernel_result =  Kernel_choose_ridge(Z, features, d, params, lambda.seq, 
                                         normalize = normalize, alpha, beta, std = std)
    weight = kernel_result$weight
    best_out = kernel_result$best_out
    Coef = kernel_result$Coef
  } else if (kernel_TF == "finite") {
    features = X
    d = NULL
    kernel_result =  Kernel_choose_ridge(Z, features, d, params, lambda.seq, 
                                         normalize = normalize, alpha, beta, std = std)
    weight = kernel_result$weight
    best_out = kernel_result$best_out
    Coef = kernel_result$Coef
  } else if (kernel_TF == "choose_sigma") {
    K = Kernel_transform_choose_sigma(X, n_sigma = n_sigma, r)
    features = list()
    d = list()
    sigma_seq = c()
    for(i in 1:n_sigma) {
      features[[i]] = K[[i]][[1]]
      d[[i]] = K[[i]][[2]]
      sigma_seq[i] = K[[i]][[3]]
    }
    kernel_list = list()
    best_diff = matrix(nrow = length(features), ncol = 2, byrow = T)
    for(i in 1:length(features)) {
      kernel_choose = Kernel_choose_ridge(Z, features[[i]], d[[i]], params, lambda.seq, 
                                          normalize = normalize, alpha, beta, std = std)
      kernel_list[[i]] = kernel_choose
      best_diff[i,] = kernel_choose[[1]]
    }
    best_sigma = c(which.min(best_diff[,1]), which.min(best_diff[,2]))
    treat_result = kernel_list[[best_sigma[1]]]
    bern_result = kernel_list[[best_sigma[2]]]
    weight = cbind(treat_result$weight[,1],bern_result$weight[,2],
                   treat_result$weight[,3],bern_result$weight[,4])
    colnames(weight) = c("treat_weight", "bernoulli_weight","treat_PS", "bernoulli_PS")
    Coef = cbind(treat_result$Coef[,1], bern_result$Coef[,2])
    colnames(Coef) = c("treat_coef", "bernoulli_coef")
    best_out = cbind(c(best_sigma[1], sigma_seq[best_sigma[1]], treat_result$best_out[1,]),
                     c(best_sigma[2], sigma_seq[best_sigma[2]], bern_result$best_out[2,]))
  }
  return(list(weight = weight, best = best_out, coef = Coef, data = DM))
}

Tailor_loss_finite_beta = function(DM, X_index, Z_index, estimand) {
  # estimand: type of estimated outcome Y ("ATE", "ATT", "ATC")
  
  Z = DM[,Z_index]
  X = DM[,X_index]
  alpha = switch (estimand, ATE = -1, ATT = 0, ATC = -1)
  beta = switch(estimand, ATE = -1, ATT = -1, ATC = 0)
  result = kernel.balance.path(Z, X, K = NULL, d = NULL,
                               alpha = alpha, beta = beta,
                               lambda = 0, tol = 1e-10)[[1]]
  coef = result$eta
  return(list(beta = coef, result = result))
}


TL_kernel_weight_all = function(DM, X_index, Z_index, r, lambda_test, Rfile_path) {
  ATE_weight = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand = "ATE",
                                           r = r, lambda_test, Rfile_path)
  ATT_weight = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand = "ATT",
                                           r = r, lambda_test, Rfile_path)
  ATC_weight = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand = "ATC",
                                           r = r, lambda_test, Rfile_path)
  weight = cbind(ATE_weight[[1]][,1], ATT_weight[[1]][,1], ATC_weight[[1]][,1:2])
  PS = cbind(ATE_weight[[1]][,3], ATT_weight[[1]][,3], ATC_weight[[1]][,3:4])
  colnames(weight) = c("ATE_weight", "ATT_weight", "ATC_weight", "bernoulli_weight")
  colnames(PS) = c("ATE_PS", "ATT_PS", "ATC_PS", "bernoulli_PS")
  return(list(weight = weight, propensity_score = PS, data = ATE_weight[[2]]))
}

#TailoredLoss_kernel_test = TL_kernel_weight_all(DM[[1]], c(4:7), 2,
 #                                                      r = 200, lambda_test = c(0, -10, 22),
  #                                                     Rfile_path = "/Users/yan/Desktop/Propensity score/Report/Evaluation/R_code/Qingyuan/covalign/R")

# Qingyuan's likelihood function

CovBalancing_ScoringRule = function(theta, X, Z, alpha, beta) {
  X = cbind(rep(1, nrow(X)), X)
  p_x = exp(X%*%theta)/(1+exp(X%*%theta))
  n = nrow(X)
  if (alpha == 0 & beta == 0) {
    theta_hat = (1/n)*sum(Z*log(p_x) + (1-Z)*(log(1-p_x)))
  } else if (alpha == -1 & beta == -1) {
    theta_hat = (1/n)*sum(Z*(log(p_x/(1-p_x))-(1/p_x)) + (1-Z)*(log((1-p_x)/p_x)-(1/(1-p_x))))
  } else if (alpha == -1 & beta == 0) {
    theta_hat = (1/n)*sum(Z*(-1/p_x) + (1-Z)*log((1-p_x)/p_x))
  } else if (alpha == 0 & beta == -1) {
    theta_hat = (1/n)*sum(Z*log(p_x/(1-p_x)) + (1-Z)*(-1/(1-p_x)))
  }
  return(theta_hat)
}

# Qingyuan's tailored loss with CV to choose lambda

TL_kernel_weight_all = function(DM, X_index, Z_index, r, lambda_test, Rfile_path) {
  ATE_weight = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand = "ATE",
                                           r = r, lambda_test, Rfile_path)
  ATT_weight = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand = "ATT",
                                           r = r, lambda_test, Rfile_path)
  ATC_weight = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand = "ATC",
                                           r = r, lambda_test, Rfile_path)
  weight = cbind(ATE_weight[[1]][,1], ATT_weight[[1]][,1], ATC_weight[[1]][,1:2])
  PS = cbind(ATE_weight[[1]][,3], ATT_weight[[1]][,3], ATC_weight[[1]][,3:4])
  colnames(weight) = c("ATE_weight", "ATT_weight", "ATC_weight", "bernoulli_weight")
  colnames(PS) = c("ATE_PS", "ATT_PS", "ATC_PS", "bernoulli_PS")
  return(list(weight = weight, propensity_score = PS, data = ATE_weight[[2]]))
}

# Estimate likelihood

Est_ScoringRule = function(X, Z, alpha, beta, init_theta, est_method) {
  #Est = optimize(CovBalancing_ScoringRule, interval = c(-1000, 1000), X = X, Z = Z, 
  #              alpha = alpha, beta = beta, maximum = T)
  Est = optim(init_theta, CovBalancing_ScoringRule, X = X, Z = Z,method = est_method, 
              alpha = alpha, beta = beta,lower = -Inf, upper = Inf,
              control=list(fnscale=-1, maxit = 1000))
  return(Est)
}

# Calculate propensity score

PS_use_ScoringRule = function(beta, X) {
  X0 = rep(1, nrow(X))
  X = cbind(X0, X)
  PS = c()
  for (i in 1:nrow(X)) {
    PS[i] = 1/(1+exp(-sum(X[i,]*beta)))
  }
  return(PS)
}
#"BFGS"
# Compute weight (from Qingyuan's R file named as "Kernel.R")
compute.weights <- function(Z, prob, alpha, beta, normalize = TRUE) {
  weights <- rep(0, length(Z))
  weights[Z == 1] <- prob[Z == 1]^alpha * (1 - prob[Z == 1])^(beta + 1)
  weights[Z == 0] <- prob[Z == 0]^(alpha + 1) * (1 - prob[Z == 0])^beta
  if (normalize) {
    n <- length(Z)
    weights[Z == 1] <- n * weights[Z == 1] / sum(weights[Z == 1])
    weights[Z == 0] <- n * weights[Z == 0] / sum(weights[Z == 0])
  }
  weights
}

Tailored_loss_finite_space = function(DM, X_index, Z_index, estimand, est_method) {
  if (estimand == "ATE") {
    alpha = -1
    beta = -1
  } else if (estimand == "ATT") {
    alpha = 0
    beta = -1
  } else if (estimand == "ATC") {
    alpha = -1
    beta = 0
  } else if (estimand == "OWATE") {
    alpha = 0
    beta = 0
  }
  X = DM[,X_index]
  Z = DM[,Z_index]
  parameters = Est_ScoringRule(X, Z, alpha, beta, rep(0,(ncol(X)+1)), est_method)
  PS = PS_use_ScoringRule(parameters$par, X)
  weight = compute.weights(Z, PS, alpha, beta, normalize = FALSE)
  return(list(weight = weight, propensity_score = PS, data = DM, parameters = parameters))
}

# Evaluate with logitic regression and CBPS and kbal
Qingyuan_finite_est = function(DM, Formula_fit, X_index, Z_index, K, est_method){
  logistic_reg = logistic_weight(DM, Formula_fit, X_index, Z_index)
  cbps_fit= CBPS_weight(DM, Formula_fit, X_index, Z_index)
  Qingyuan_fit = Tailored_loss_finite_space(DM, X_index, Z_index, "ATE", est_method)
  if (K > 1) {
    ATE_est = c(ATE_infer(logistic_reg[[2]], logistic_reg[[1]][,1], X_index, Z_index),
                ATE_infer(cbps_fit[[2]], cbps_fit[[1]][,1], X_index, Z_index),
                ATE_infer(Qingyuan_fit[[3]], Qingyuan_fit[[1]], X_index, Z_index),
                0)
    names(ATE_est) = c("logistic", "CBPS", "Qingyuan_finite", "null")
  } else {
    kbl_fit = KBAL_weight(DM, X_index, Z_index, "ebal")
    ATE_est = c(ATE_infer(logistic_reg[[2]], logistic_reg[[1]][,1], X_index, Z_index),
                ATE_infer(cbps_fit[[2]], cbps_fit[[1]][,1], X_index, Z_index),
                ATE_infer(Qingyuan_fit[[3]], Qingyuan_fit[[1]], X_index, Z_index))
    Z = DM[, Z_index]
    kbl_att = ATE_infer(kbl_fit[[2]], kbl_fit[[1]][,1], X_index, Z_index)
    kbl_atc = ATE_infer(kbl_fit[[2]], kbl_fit[[1]][,2], X_index, Z_index)
    kbl_ate = mean(Z)*kbl_att + (1-mean(Z))*kbl_atc
    ATE_est = c(ATE_est, kbl_ate)
    names(ATE_est) = c("logistic", "CBPS", "Qingyuan_finite","kbal")
  }
  logistic_sd = Standardized_diff(logistic_reg[[2]], logistic_reg[[1]][,1], X_index, Z_index, "ATE")
  cbps_sd = Standardized_diff(cbps_fit[[2]], cbps_fit[[1]][,1], X_index, Z_index, "ATE")
  Qingyuan_sd = Standardized_diff(Qingyuan_fit[[3]], Qingyuan_fit[[1]], X_index, Z_index, "ATE")
  Coef_all = cbind(logistic_reg[[3]], cbps_fit[[4]], Qingyuan_fit[[4]]$par)
  Coef_qingyuan = Qingyuan_fit[[4]]
  colnames(Coef_all) = c("logistic_reg", "exact CBPS ATE", "exact CBPS ATT", 
                         "exact CBPS ATC", "over CBPS ATE", 
                       "over CBPS ATT", "over CBPS ATC", "Qingyuan_finite")
  return(list(ATE_est = ATE_est, logistic_sd = logistic_sd[[1]],
              CBPS_sd = cbps_sd[[1]], Qingyuan_sd = Qingyuan_sd[[1]], 
              Coef_all = Coef_all, Coef_qingyuan = Coef_qingyuan))
}


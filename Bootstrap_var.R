CBPS_ATE_est = function(DM, Formula_fit, X_index, Z_index) {
  require(CBPS)
  DM_XZ = data.frame(DM[,c(Z_index, X_index)])
  #form = formula(Formula_fit)
  cbps_ATE = CBPS(Formula_fit, data = DM_XZ, ATT = 0, method = "exact", standardize = FALSE)
  weight = cbps_ATE$weights
  ATE = ATE_infer(DM, weight, X_index, Z_index, normalize = T) 
  return(list(ATE_est = ATE, weight = weight, fit = cbps_ATE))
}

CBPS_boot_ATE_variance = function(DM, Formula_fit, X_index, Z_index, B, seeds) {
  ATE_point = CBPS_ATE_est(DM, Formula_fit, X_index, Z_index)$ATE_est
  n = nrow(DM)
  set.seed(seeds)
  ATE_vec = c()
  for(i in 1:B) {
    I = sample(1:n, replace = T)
    DM_B = DM[I,] 
    ATE_vec[i] = CBPS_ATE_est(DM_B, Formula_fit, X_index, Z_index)$ATE_est
  }
  ATE_sd = sd(ATE_vec)
  ATE_CI = ATE_point + c(-1.96,1.96)*ATE_sd
  ATE = c(ATE_CI, ATE_point, ATE_sd)
  names(ATE) = c("CI lower", "CI upper", "point", "sd")
  return(list(ATE = ATE, ATE_vec = ATE_vec))
}

CBPS_boot_wrapper = function(iter, n, K, delta, beta, h, delta_function,
                             Formula_fit, X_index, Z_index, B, seeds) {
  ATE_CI = matrix(nrow = iter, ncol = 4, byrow = T)
  for(i in 1:iter) {
    DM = PS_simulation_new(n, 0.5, K, beta, alpha = matrix(c(0.5,1,0.6,2.2,-1.2,rep(0,10),1,0,-1,0,1,0,-1,0),nrow = 23), 
                           delta = delta, Homo = h, delta_function_1, seeds[i])
    fit = CBPS_boot_ATE_variance(DM[[1]], Formula_fit, X_index, Z_index, B, seeds[i])
    ATE_CI[i, ] = fit$ATE
  }
  return(ATE_CI)
}



#True_homo = CBPS_boot_wrapper(100, 1000, 8, h = 0, beta = matrix(c(-1.5,0.5,-0.75,2,-0.5, rep(0,10), beta_12),nrow = 23), 
 #                     delta_function_1, Cov_12, c(4:7,18:25), 2, 200, seeds)











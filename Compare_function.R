CBPS_Qingyuan_compare = function(DM, Formula_fit, X_index, Z_index, treatI,
                                 estimand, r=100, lambda_test = c(0, -10, 22), 
                                 normalize = FALSE, n_sigma) {
  logistic_reg = logistic_weight(DM, Formula_fit, X_index, Z_index)
  cbps_fit= CBPS_weight(DM, Formula_fit, X_index, Z_index)
  Qingyuan_normdiff = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand, kernel_TF= "finite", r, 
                                         lambda_test, normalize,n_sigma, std = "std.norm") # kernel_TF: "kernel", "finite", "choose_sigma"
  Qingyuan_stddiff = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand, kernel_TF= "finite", r, 
                                                  lambda_test, normalize,n_sigma, std = "std.diff") # kernel_TF: "kernel", "finite", "choose_sigma"
  Qingyuan_kern_sigma_normdiff = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand, kernel_TF="choose_sigma", r, 
                                              lambda_test, normalize, n_sigma, std = "std.norm") # std = c("std.diff", "std.norm")
  Qingyuan_kern_sigma_stddiff = Tailored_loss_kernel_weight(DM, X_index, Z_index, estimand, kernel_TF="choose_sigma", r, 
                                                    lambda_test, normalize, n_sigma, std = "std.diff")
  kbal = KBAL_weight(DM, X_index, Z_index, "ebal")
  entropy = entropy_balancing_weight(DM, X_index, Z_index)
  Coef_all = cbind(logistic_reg[[3]], cbps_fit[[4]][,1], 
                   Qingyuan_normdiff[[3]][,1], Qingyuan_stddiff[[3]][,1])
  colnames(Coef_all) = c("logistic_reg", "exact CBPS ATE", "QY_norm", "QY_stddiff")
  ATE_est = c(ATE_infer(logistic_reg[[2]], logistic_reg[[1]][,1], X_index, Z_index),
              ATE_infer(cbps_fit[[2]], cbps_fit[[1]][,1], X_index, Z_index),
              ATE_infer(Qingyuan_normdiff[[4]], Qingyuan_normdiff[[1]][,1], X_index, Z_index),
              ATE_infer(Qingyuan_stddiff[[4]], Qingyuan_stddiff[[1]][,1], X_index, Z_index),
              ATE_infer(Qingyuan_kern_sigma_normdiff[[4]], Qingyuan_kern_sigma_normdiff[[1]][,1], X_index, Z_index),
              ATE_infer(Qingyuan_kern_sigma_stddiff[[4]], Qingyuan_kern_sigma_stddiff[[1]][,1], X_index, Z_index),
              ATE_est_averageATT.ATC(kbal, X_index, Z_index)[1], ATE_est_averageATT.ATC(entropy, X_index, Z_index)[1])
  names(ATE_est) = c("logistic", "CBPS","QY_norm","QY_stddiff","QY_kernel_norm","QY_kernel_stddiff","kbal","entropy")
  W = cbind(logistic_reg[[1]][,treatI], cbps_fit[[1]][,treatI], Qingyuan_normdiff[[1]][,1],
            Qingyuan_stddiff[[1]][,1], Qingyuan_kern_sigma_normdiff[[1]][,1], 
            Qingyuan_kern_sigma_stddiff[[1]][,1], kbal[[1]], entropy[[1]])
  colnames(W) = c("log_IPW","cbps_IPW", "QY_ipw_norm","QY_ipw_stddiff","QY_kern_ipw_norm",
                  "QY_kern_ipw_stddiff","kbal_att","kbal_atc","entropy_att","entropy_atc")
  PS = cbind(logistic_reg[[1]][,4], cbps_fit[[3]][,treatI], Qingyuan_normdiff[[1]][,3], 
             Qingyuan_stddiff[[1]][,3],Qingyuan_kern_sigma_normdiff[[1]][,3], 
             Qingyuan_kern_sigma_stddiff[[1]][,3])
  colnames(PS) = c("logistic_PS","cbps_PS","QYnorm_PS","QYstddiff_PS","QY_kern_norm_PS","QY_kern_stddiff_PS")
  logistic_sd = Standardized_diff(logistic_reg[[2]], logistic_reg[[1]][,1], X_index, Z_index, "ATE")
  cbps_sd = Standardized_diff(cbps_fit[[2]], cbps_fit[[1]][,1], X_index, Z_index, "ATE")
  QY_sd_1 = Standardized_diff(Qingyuan_normdiff[[4]], Qingyuan_normdiff[[1]][,1], X_index, Z_index, "ATE")
  QY_sd_2 = Standardized_diff(Qingyuan_stddiff[[4]], Qingyuan_stddiff[[1]][,1], X_index, Z_index, "ATE")
  QY_sd_3 = Standardized_diff(Qingyuan_kern_sigma_normdiff[[4]], Qingyuan_kern_sigma_normdiff[[1]][,1], X_index, Z_index, "ATE")
  QY_sd_4 = Standardized_diff(Qingyuan_kern_sigma_stddiff[[4]], Qingyuan_kern_sigma_stddiff[[1]][,1], X_index, Z_index, "ATE")
  QY_sd = cbind(QY_sd_1, QY_sd_2[,6], QY_sd_3[,6], QY_sd_4[,6])
  colnames(QY_sd) = c("Control_mean","Control_sd","Case_mean","Case_sd",
                      "Before_S/D","QY_Weighted_norm_S/D", "QY_Weighted_stddiff_S/D",
                      "QY_kernel_norm_Weighted_S/D","QY_kernel_stddiff_Weighted_S/D")
  return(list(ATE = ATE_est, Coef_all = Coef_all, propensity_score = PS,
              logistic_sd = logistic_sd, CBPS_sd = cbps_sd, QY_sd = QY_sd, weight = W))
}

# treatI = 1 for ATE

CBPS_Qingyuan_wrapper = function(iter, n, K, delta, delta_function, seeds, beta, alpha,
                                 fit_form, h, X_index, Z_index, 
                                 estimand, r=100, lambda_test = c(0, -10, 22), 
                                 normalize = FALSE,
                                 specific_kernel=Normal_kernel, c=c(0.05,0.1), n_sigma) {
  require(generalhoslem)
  # h = 0 for homogeneous effect
  # h = 1 for heterogeneous effect
  log_sd = list()
  CBPS_sd = list()
  QY_sd = list()
  ATE_est = matrix(nrow = iter, ncol = 8, byrow = T)
  Hosmer_Lemeshow = matrix(nrow = iter, ncol = 6, byrow = T)
  Specific_test = matrix(nrow = iter, ncol = 12, byrow = T)
  data_out = list()
  Coef_all = list()
  PS_weight = list()
  for(i in 1:iter) {
    DM = PS_simulation_new(n, 0.5, K, beta, alpha, delta, Homo = h, delta_function, seeds[i])
    data_out[[i]] = cbind(DM[[1]], DM[[2]])
    Fit_comp = CBPS_Qingyuan_compare(DM[[1]], fit_form, X_index, Z_index, 1,
                                     estimand, r, lambda_test, normalize, n_sigma)
    ATE_est[i,] = Fit_comp[[1]]
    Coef_all[[i]] = Fit_comp[[2]]
    PS_weight[[i]] = cbind(Fit_comp[[3]], Fit_comp[[7]])
    colnames(PS_weight[[i]]) = c("log_PS","cbps_PS","QYnorm_PS","QYstddiff_PS","QY_kern_norm_PS","QY_kern_stddiff_PS",
                                 "log_ipw","cbps_ipw","QYnorm_ipw","QYstddiff_ipw","QY_kern_norm_ipw", "QY_kern_stddiff_ipw",
                                 "kbal_att","kbal_atc","entropy_att","entropy_atc")
    log_sd[[i]] = Fit_comp[[4]]
    CBPS_sd[[i]] = Fit_comp[[5]]
    QY_sd[[i]] = Fit_comp[[6]]
    Hosmer_Lemeshow[i,] = c(as.numeric(logitgof(DM[[1]][,Z_index], Fit_comp[[3]][,1], g = 10)$stat),
                            as.numeric(logitgof(DM[[1]][,Z_index], Fit_comp[[3]][,2], g = 10)$stat),
                            as.numeric(logitgof(DM[[1]][,Z_index], Fit_comp[[3]][,3], g = 10)$stat),
                            as.numeric(logitgof(DM[[1]][,Z_index], Fit_comp[[3]][,4], g = 10)$stat),
                            as.numeric(logitgof(DM[[1]][,Z_index], Fit_comp[[3]][,5], g = 10)$stat),
                            as.numeric(logitgof(DM[[1]][,Z_index], Fit_comp[[3]][,6], g = 10)$stat))
    Specific_test[i,] = c(TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[1], Fit_comp[[3]][,1], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[1], Fit_comp[[3]][,2], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[1], Fit_comp[[3]][,3], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[1], Fit_comp[[3]][,4], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[1], Fit_comp[[3]][,5], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[1], Fit_comp[[3]][,6], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[2], Fit_comp[[3]][,1], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[2], Fit_comp[[3]][,2], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[2], Fit_comp[[3]][,3], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[2], Fit_comp[[3]][,4], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[2], Fit_comp[[3]][,5], Z_index),
                          TestforPS_CDFonParticipation(DM[[1]], specific_kernel,c = c[2], Fit_comp[[3]][,6], Z_index))
    colnames(ATE_est) = c("logistic", "CBPS", "QY_norm","QY_stddiff",
                          "QY_kern_norm","QY_kern_stddiff", "kbal", "entropy")
    colnames(Hosmer_Lemeshow) =  c("logistic", "CBPS", "QY_norm","QY_stddiff","QY_kern_norm","QY_kern_stddiff")
    colnames(Specific_test) = c(paste0(0.05, colnames(Hosmer_Lemeshow)), paste0(0.1, colnames(Hosmer_Lemeshow)))
  }
  return(list(ATE_est = ATE_est,Coef = Coef_all,PS_weight = PS_weight, 
              logistic_sd = log_sd, cbps_sd = CBPS_sd, QY_sd = QY_sd, HS = Hosmer_Lemeshow, 
              Specific_test = Specific_test, data_out = data_out))
}

Normal_kernel = function(u) {(1/sqrt(2*pi)*exp(-(u^2)/2))}

Cov_generate = function(k, fitted_model = c("linear","power","2nd"),
                        Cov_vec = c("X1","X2","X3","X4","Norm1","Norm2","Norm3","Norm4","Bin1","Bin2","Bin3","Bin4")) {
  if (fitted_model == "linear") {
    Cov = paste0("Z~X1+X2+X3+X4")
    for(i in 1:k) {
      Cov = paste0(Cov, "+", "Norm", i)
    }
    for(i in 1:k) {
      Cov = paste0(Cov, "+","Bin", i)
    }
  } else if(fitted_model == "power") {
    Cov = paste0("Z~X1+X2+X3+X4")
    Cov = paste0("Z~X1+X2+X3+X4")
    for(i in 1:k) {
      Cov = paste0(Cov, "+", "Norm", i)
    }
    for(i in 1:k) {
      Cov = paste0(Cov, "+","Bin", i)
    }
    for(i in 1:2) {
      Cov = paste0(Cov, "+", "I(X", i, "^2)")
    }
    for(i in 1:k) {
      Cov = paste0(Cov, "+", "I(Norm", i, "^2)")
    }
  } else if(fitted_model == "2nd") {
    Cov = paste0("Z~X1+X2+X3+X4")
    for(i in 1:k) {
      Cov = paste0(Cov, "+", "Norm", i)
    }
    for(i in 1:k) {
      Cov = paste0(Cov, "+","Bin", i)
    }
    for(i in 1:2) {
      Cov = paste0(Cov, "+", "I(X", i, "^2)")
    }
    for(i in 1:k) {
      Cov = paste0(Cov, "+", "I(Norm", i, "^2)")
    }
    for (i in 2:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[1],":",Cov_vec[i])
    }
    for (i in 3:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[2],":",Cov_vec[i])
    }
    for (i in 4:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[3],":",Cov_vec[i])
    }
    for (i in 5:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[4],":",Cov_vec[i])
    }
    for (i in 6:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[5],":",Cov_vec[i])
    }
    for (i in 7:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[6],":",Cov_vec[i])
    }
    for (i in 8:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[7],":",Cov_vec[i])
    }
    for (i in 9:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[8],":",Cov_vec[i])
    }
    for (i in 10:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[9],":",Cov_vec[i])
    }
    for (i in 11:length(Cov_vec)) {
      Cov = paste0(Cov, "+", Cov_vec[10],":",Cov_vec[i])
    }
    Cov = paste0(Cov, "+", Cov_vec[11],":",Cov_vec[12])
  }
  return(Cov)
}

linear_cov_generate = function(K) {
  Cov_comp = c(paste0("X",seq(1,4)), paste0("Norm", seq(1,K/2)), paste0("Bin", seq(1,K/2)))
  if(K == 0) {
    Cov = "Z~X1+X2+X3+X4"
  } else {
    Cov = paste0("Z~", Cov_comp[1])
    for(i in 2:length(Cov_comp)) {
      Cov = paste0(Cov, "+", Cov_comp[i])
    }
  }
  return(Cov)
}





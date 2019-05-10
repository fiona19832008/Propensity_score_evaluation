# Weighted variance
weightedVariance <- function(x, w) {
  w <- w[i <- !is.na(x)]
  x <- x[i]
  sum.w <- sum(w)
  out = (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
  return(out)
}

# Standardized difference
Standardized_diff <- function(DM, weight, X_index, Z_index, estimand) {
  # weight: a vector of weight for each sample
  #library(SciencesPo)
  library(Matching)
  W_X = cbind(weight, DM[,X_index])
  Treated = W_X[which(DM[,Z_index]==1), ]
  Control = W_X[which(DM[,Z_index]==0), ]
  Out_before = matrix(nrow = (ncol(W_X) - 1), ncol = 5, byrow = T)
  #s = seq(from = I1, to = I2)
  for (i in 1:(ncol(W_X) - 1)) {
    Out_before[i,1] = mean(Control[,i+1])
    Out_before[i,2] = sd(Control[,i+1])
    Out_before[i,3] = mean(Treated[,i+1])
    Out_before[i,4] = sd(Treated[,i+1])
    Out_before[i,5] = (abs(mean(Control[,i+1]) - mean(Treated[,i+1]))/sd(W_X[,i+1]))*100
  }
  colnames(Out_before) <- c("Control_mean", "Control_sd", "Case_mean", "Case_sd", "Before_S/D")
  if (estimand == "ATE") {
    Out_after = ATE_Weight_SD(Treated, Control, 2, ncol(W_X), 1)
  } else if (estimand == "ATT") {
    Out_after = matrix(nrow = (col(W_X) - 1), ncol = 5, byrow = T)
    for (i in 1:(col(W_X) - 1)) {
      StdD = balanceUV(Treated[,i+1], control[,i+1], estimand = "ATT",
                       weights.Tr = Treated[,i], weights.Co = Control[,1])
      Out_after[i,] = c(StdD$mean.Co, sqrt(StdD$var.Co), StdD$mean.Tr, 
                        sqrt(StdD$var.Tr), StdD$sdiff.pooled)
    }
    colnames(Out_after) = c("control.mean", "control.sd", "treated.mean", "treated.sd", "weighted.S/D")
  } else if (estimand == "ATC") {
    Out_after = matrix(nrow = (col(W_X) - 1), ncol = 5, byrow = T)
    for (i in 1:(col(W_X) - 1)) {
      StdD = balanceUV(Treated[,i+1], control[,i+1], estimand = "ATC",
                       weights.Tr = Treated[,i], weights.Co = Control[,1])
      Out_after[i,] = c(StdD$mean.Co, sqrt(StdD$var.Co), StdD$mean.Tr, 
                        sqrt(StdD$var.Tr), StdD$sdiff.pooled)
    }
    colnames(Out_after) = c("control.mean", "control.sd", "treated.mean", "treated.sd", "weighted.S/D")
  }
  Out = cbind(Out_before, Out_after[,5])
  colnames(Out) = c("Control_mean", "Control_sd", "Case_mean", 
                    "Case_sd", "Before_S/D", "Weighted_S/D")
  return(Standerdized_diff = Out)
}

ATE_Weight_SD <- function(case, control, I1, I2, W_I) {
  #library(SciencesPo)
  s = seq(from = I1, to = I2)
  out = matrix(nrow = (I2-I1+1), ncol = 4, byrow = T)
  for (i in 1:(I2-I1+1)) {
    out[i,1] = sum(control[,s[i]]*control[,W_I])/sum(control[,W_I])
    out[i,2] = sqrt(weightedVariance(control[,s[i]], control[,W_I]))
    out[i,3] = sum(case[,s[i]]*case[,W_I])/sum(case[,W_I])
    out[i,4] = sqrt(weightedVariance(case[,s[i]], case[,W_I]))
  }
  Std_Diff = c()
  for (i in 1:nrow(out)) {
    Std_Diff[i] = (out[i,1] - out[i,3])/sqrt(((out[i,2]^2) + (out[i,4]^2))/2)
  }
  Std_Diff = Std_Diff*100
  out = cbind(out, Std_Diff)
  colnames(out) = c("control.mean", "control.sd", "treated.mean", "treated.sd", "weighted.S/D")
  return(out)
}

ATE_Weight_diff <- function(DM, weight, X_index, Z_index, std) {
  # std = c("std.diff", "std.norm")
  W_X = cbind(weight, DM[,X_index])
  case = W_X[which(DM[,Z_index]==1), ]
  control = W_X[which(DM[,Z_index]==0), ]
  I1 = 2 
  I2 = ncol(W_X)
  W_I = 1
  s = seq(from = I1, to = I2)
  out = matrix(nrow = (I2-I1+1), ncol = 2, byrow = T)
  sd = out
  for (i in 1:(I2-I1+1)) {
    out[i,1] = sum(control[,s[i]]*control[,W_I])/sum(control[,W_I])
    out[i,2] = sum(case[,s[i]]*case[,W_I])/sum(case[,W_I])
    sd[i,1] = weightedVariance(control[,s[i]], control[,W_I])
    sd[i,2] = weightedVariance(case[,s[i]], case[,W_I])
  }
  #Std_Diff = out[,1] - out[,2]
  if (std == "std.diff") {
    Std_Diff = c()
    for (i in 1:nrow(out)) {
      Std_Diff[i] = (out[i,1] - out[i,2])/sqrt(((sd[i,1]) + (sd[i,2]))/2)
    }
    Std_Diff = Std_Diff*100
  } else if (std == "std.norm") {
    Std_Diff = out[,1] - out[,2]
  }
  out = cbind(out, Std_Diff)
  colnames(out) = c("control.mean", "treated.mean", "weighted.diff")
  return(out)
}

# ATE estimation
ATE_infer= function(DM, weight, X_index, Z_index, normalize = T) {
  # always put y at the first column
  DM_weight = cbind(DM, weight)
  W_I = ncol(DM_weight)
  Treated = DM_weight[which(DM_weight[,Z_index]==1), ]
  control = DM_weight[which(DM_weight[,Z_index]==0), ]
  if(normalize == T) {
    ATE_case = sum(Treated[,1]*Treated[,W_I])/sum(Treated[,W_I]) 
    ATE_control = sum(control[,1]*control[,W_I])/sum(control[,W_I])
  } else if(normalize == F) {
    ATE_case = mean(Treated[,1]*Treated[,W_I]) 
    ATE_control = mean(control[,1]*control[,W_I])
  }
  ATE_after = ATE_case - ATE_control
  return(ATE_after)
}

ATE_est_averageATT.ATC = function(kbl_fit, X_index, Z_index, est_method){
  DM = kbl_fit[[2]]
  Z = DM[, Z_index]
  kbl_att = ATE_infer(DM, kbl_fit[[1]][,1], X_index, Z_index)
  kbl_atc = ATE_infer(DM, kbl_fit[[1]][,2], X_index, Z_index)
  kbl_ate = mean(Z)*kbl_att + (1-mean(Z))*kbl_atc
  out = c(kbl_ate, kbl_att, kbl_atc)
  return(out)
}

Find_bias_variance_rmse = function(DataVector, True_value) {
  bias = 100*(mean(DataVector) - True_value)/True_value
  variance = var(DataVector)
  rmse = sqrt(mean((DataVector - True_value)^2))
  stat_out = c(bias, rmse, variance)
  names(stat_out) = c("bias", "RMSE", "Var")
  return(out = stat_out)
}

Col_bias_variance_rmse = function(DM, True_value) {
  if (length(True_value) == 1) {
    True_value_vec = rep(True_value, ncol(DM))
  } else {
    True_value_vec = True_value
  }
  out = matrix(nrow = 1, ncol = 3*ncol(DM))
  I = seq(1,3*ncol(DM), 3)
  for (i in 1:ncol(DM)) {
    out[,I[i]:(I[i]+2)] = Find_bias_variance_rmse(DM[,i], True_value_vec[i]) 
  }
  colnames(out) = rep(c("bias","RMSE","Var"), ncol(DM))
  return(out)
}

Get_convergence = function(Dlist) {
  p = c()
  for(i in 1:length(Dlist[[6]])) {
    p[i] = Dlist[[6]][[i]]$convergence
  }
  return(p)
}

# Calculate Hosmer-Lemeshow statistics for a list of data
Hosmer_Lemeshow = function(DM, no_quantile, Z_index, glm_ps, cbps_ps) {
  # k: number of parameters
  require(generalhoslem)
  n = length(DM)
  HL_stat = matrix(nrow = n, ncol = 2, byrow = T)
  for (i in 1:n) {
    HL_stat[i,] = c(as.numeric(logitgof(DM[[i]][,Z_index], DM[[i]][,glm_ps], g = no_quantile)$stat),
                    as.numeric(logitgof(DM[[i]][,Z_index], DM[[i]][,cbps_ps], g = no_quantile)$stat))
  }
  colnames(HL_stat) = c("glm", "cbps")
  return(list(Hosmer_Lemeshow = HL_stat))
}

# A specification test for the propensity score using its distribution conditional on participation

Normal_kernel = function(u) {(1/sqrt(2*pi)*exp(-(u^2)/2))}

TestforPS_CDFonParticipation <- function(DM, Kernel_function, c, PS_hat, Z) {
  n = nrow(DM)
  h = c*(n^(-1/8))
  Q = PS_hat
  e = DM[,Z] - Q
  V = 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (i!=j) {
        Vij = Kernel_function((Q[i] - Q[j])/h)*e[i]*e[j]
        V_ij = Vij/h
        V = V + V_ij
      }
    }
  }
  V_n = V/(n*(n-1))
  return(test_stat = V_n)
}




# Logistic regression
logistic_weight = function(DM, Formula_fit, X_index, Z_index) {
  form = formula(Formula_fit)
  X = data.frame(DM[,X_index])
  Z = data.frame(DM[,Z_index])
  DM_XZ = data.frame(DM[,c(Z_index, X_index)])
  glm = glm(form, data = DM_XZ, family = binomial(link = "logit"))
  propensity_score = glm$fitted.values
  IPW = Z*(1/propensity_score) + (1-Z)*(1/(1-propensity_score))
  ATT_weight = Z+(1-Z)*(propensity_score/(1-propensity_score))
  ATC_weight = Z*((1-propensity_score)/propensity_score) + (1-Z)
  Out_weight = cbind(IPW, ATT_weight, ATC_weight, propensity_score)
  colnames(Out_weight) = c("IPW", "ATT_weight", "ATC_weight", "PS")
  Coef = as.numeric(glm$coefficients)
  return(list(weight = Out_weight, data = DM, coefficient = Coef))
}

# CBPS
CBPS_weight = function(DM, Formula_fit, X_index, Z_index) {
  require(CBPS)
  DM_XZ = data.frame(DM[,c(Z_index, X_index)])
  cbps1_ATE = CBPS(Formula_fit, data = DM_XZ, ATT = 0, method = "exact", standardize = FALSE)
  cbps2_ATE = CBPS(Formula_fit, data = DM_XZ, ATT = 0, method = "over", standardize = FALSE)
  cbps1_ATT = CBPS(Formula_fit, data = DM_XZ, ATT = 1, method = "exact", standardize = FALSE)
  cbps2_ATT = CBPS(Formula_fit, data = DM_XZ, ATT = 1, method = "over", standardize = FALSE)
  cbps1_ATC = CBPS(Formula_fit, data = DM_XZ, ATT = 2, method = "exact", standardize = FALSE)
  cbps2_ATC = CBPS(Formula_fit, data = DM_XZ, ATT = 2, method = "over", standardize = FALSE)
  Coef = cbind(as.numeric(cbps1_ATE$coefficients), as.numeric(cbps1_ATT$coefficients),
               as.numeric(cbps1_ATC$coefficients), as.numeric(cbps2_ATE$coefficients),
               as.numeric(cbps2_ATT$coefficients), as.numeric(cbps2_ATC$coefficients))
  propensity_score = cbind(cbps1_ATE$fitted.values, cbps1_ATT$fitted.values,
                           cbps1_ATC$fitted.values, cbps2_ATE$fitted.values,
                           cbps2_ATT$fitted.values, cbps2_ATC$fitted.values)
  colnames(Coef) = colnames(propensity_score) = c("exact CBPS ATE", "exact CBPS ATT", 
                                                  "exact CBPS ATC", "over CBPS ATE", 
                                                  "over CBPS ATT", "over CBPS ATC")
  Out_weight = cbind(cbps1_ATE$weights, cbps1_ATT$weights, cbps1_ATC$weights,
                     cbps2_ATE$weights, cbps2_ATT$weights, cbps2_ATC$weights)
  colnames(Out_weight) = c("IPW", "ATT_weight", "ATC_weight",
                           "IPW_over", "ATT_weight_over", "ATC_weight_over")
  return(list(weight = Out_weight, data_set = DM, propensity_score = propensity_score,
              coef = Coef))
}


#source("/Users/yan/Desktop/Propensity score/Report/Report_4/PS_functions_new.R")

#log_test = logistic_weight(DM[[1]], "Z~X1+X2+X3+X4", c(4:7), 2)
#cbps_test = CBPS_weight(DM[[1]], "Z~X1+X2+X3+X4", c(4:7), 2)

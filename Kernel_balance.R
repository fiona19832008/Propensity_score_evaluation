#install.packages('devtools', repos = 'http://cran.us.r-project.org')
#devtools::install_github('chadhazlett/KBAL')
#devtools::install_github('xuyiqing/tjbal')

KBAL_weight = function(DM, X_index, Z_index, est_method) {
  library(KBAL)
  Z_att = DM[,Z_index]
  X = as.matrix(DM[,X_index])
  kbal_att = kbal(D=Z_att, X=X, method = est_method)
  Z_atc = ifelse(Z_att==1, 0 , 1)
  kbal_atc = kbal(D=Z_atc, X=X, method = est_method)
  Out_weight = cbind(kbal_att$w, kbal_atc$w)
  colnames(Out_weight) = c("ATT_weight", "ATC_weight")
  return(list(weight = Out_weight, data = DM, 
              kbal_ATT = kbal_att, kbal_ATC = kbal_atc))
}

#kbal_test = KBAL_weight(DM[[1]], c(4:7), 2, "ebal")
#kbal_test_1 = KBAL_weight(DM[[1]], c(4:7), 2, "el")

# entropy balancing
entropy_balancing_weight = function(DM, X_index, Z_index) {
  require(ebal)
  Z_att = DM[,Z_index]
  X = as.matrix(DM[,X_index])
  ebal_att = ebalance(Z_att, X)
  Z_atc = ifelse(Z_att==1, 0 , 1)
  ebal_atc = ebalance(Z_atc,X)
  att_weight = rep(1, length(Z_att))
  atc_weight = att_weight
  att_weight[which(Z_att==0)] = ebal_att$w
  atc_weight[which(Z_att==1)] = ebal_atc$w
  Out_weight = cbind(att_weight, atc_weight)
  colnames(Out_weight) = c("ATT_weight", "ATC_weight")
  return(list(weight = Out_weight, data = DM, 
              ebal_ATT = ebal_att, ebal_ATC = ebal_atc))
}




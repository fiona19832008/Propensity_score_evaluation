True_delta_f1_4 = Find_true_delta_new(1000000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,10)),nrow = 15),delta_function_1)
True_delta_f2_12 = Find_true_delta_new(1000000, 0.5, 8,
                                     beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,10),beta_12),nrow = 23), delta_function_2)

True_delta_f2_100 = Find_true_delta_new(1000000, 0.5, 96,
                                      beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,10),beta_100),nrow = 111), delta_function_2)

Gamma_1_f2_4 = Find_true_delta_new(1000000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,0,1,rep(0,8)),nrow = 15),delta_function_1)
Gamma_1_f2_12 = Find_true_delta_new(1000000, 0.5, 8,
                                  beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,0,1,rep(0,8),beta_12),nrow = 23), delta_function_2)
Gamma_1_f1_100 = Find_true_delta_new(1000000, 0.5, 96,
                                    beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,0,1,rep(0,8),beta_100),nrow = 111), delta_function_1)

Gamma_5_f1_4 = Find_true_delta_new(1000000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,0,5,rep(0,8)),nrow = 15),delta_function_1)
Gamma_5_f1_12 = Find_true_delta_new(1000000, 0.5, 8,
                                   beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,0,5,rep(0,8),beta_12),nrow = 23), delta_function_1)
Gamma_5_f1_100 = Find_true_delta_new(1000000, 0.5, 96,
                                    beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,0,5,rep(0,8),beta_100),nrow = 111), delta_function_1)

X1_X2_f1_12 = Find_true_delta_new(1000000, 0.5, 8,
                                 beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,2,rep(0,5),1,1,0,0,beta_12),nrow = 23), delta_function_1)
X1_X2_f1_100 = Find_true_delta_new(1000000, 0.5, 96,
                                  beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,2,rep(0,5),1,1,0,0,beta_100),nrow = 111), delta_function_1)

#n, p1, K, beta,delta_function

E = matrix(c(True_delta_f1_4, True_delta_f1_12,
            True_delta_f1_100, Gamma_1_f1_4,
           Gamma_1_f1_12, Gamma_1_f1_100,
          Gamma_5_f1_4,
         Gamma_5_f1_12, Gamma_5_f1_100,
        X1_X2_f1_12, X1_X2_f1_100), nrow = 11)
rownames(E) = c("True_delta_f1_4","True_delta_f1_12",
               "True_delta_f1_100", "Gamma_1_f1_4",
              "Gamma_1_f1_12", "Gamma_1_f1_100",
             "Gamma_5_f1_4",
            "Gamma_5_f1_12", "Gamma_5_f1_100",
           "X1_X2_f1_12", "X1_X2_f1_100")

beta_12 = rep(c(1,-1,0,0),2)
alpha_12 = rep(c(1,0,-1,0),2)

beta_20 = rep(c(1,-1,1,-1,rep(0,4)),2)
alpha_20 = rep(c(1,-1,0,0,1,-1,0,0),2)

#beta_32 = rep(c(rep(1,2),rep(-1,2),rep(1,2),rep(-1,2),rep(0,8)),2)
#alpha_32 = rep(c(rep(1,2),rep(-1,2), rep(0,4), rep(1,2),rep(-1,2), rep(0,4)),2)

#beta_64 = rep(c(rep(1,4),rep(-1,4),rep(1,4),rep(-1,4),rep(0,16)),2)
#alpha_64 = rep(c(rep(1,4),rep(-1,4), rep(0,8), rep(1,4),rep(-1,4), rep(0,8)),2)

#beta_80 = rep(c(rep(1,5),rep(-1,5),rep(1,5),rep(-1,5),rep(0,20)),2)
#alpha_80 = rep(c(rep(1,5),rep(-1,5), rep(0,10), rep(1,5),rep(-1,5), rep(0,10)),2)

#beta_100 = rep(c(rep(1,6),rep(-1,6),rep(1,6),rep(-1,6),rep(0,24)),2)
#alpha_100 = rep(c(rep(1,6),rep(-1,6), rep(0,12), rep(1,6),rep(-1,6), rep(0,12)),2)

# Find true delta new simulation scheme

beta_model = rbind(rep(0,12), c(rep(0,6),2,rep(0,5)), c(0,5,rep(0,10)), c(4,rep(0,5),2,2,0,0,0,0), c(rep(0,10),1,0))

True_delta_h1_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,12)),nrow = 17),Hetero_f1)
True_delta_h1_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,12),beta_12),nrow = 17+8), Hetero_f1)
True_delta_h2_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,12)),nrow = 17),Hetero_f2)
True_delta_h2_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,rep(0,12),beta_12),nrow = 17+8), Hetero_f2)

NonLinear_delta_h1_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[2,]),nrow = 17),Hetero_f1)
NonLinear_delta_h1_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[2,],beta_12),nrow = 17+8), Hetero_f1)
NonLinear_delta_h2_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[2,]),nrow = 17),Hetero_f2)
NonLinear_delta_h2_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[2,],beta_12),nrow = 17+8), Hetero_f2)
NonAdd_delta_h1_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[3,]),nrow = 17),Hetero_f1)
NonAdd_delta_h1_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[3,],beta_12),nrow = 17+8), Hetero_f1)
NonAdd_delta_h2_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[3,]),nrow = 17),Hetero_f2)
NonAdd_delta_h2_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[3,],beta_12),nrow = 17+8), Hetero_f2)
Both_delta_h1_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[4,]),nrow = 17),Hetero_f1)
Both_delta_h1_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[4,],beta_12),nrow = 17+8), Hetero_f1)
Both_delta_h2_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[4,]),nrow = 17),Hetero_f2)
Both_delta_h2_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[4,],beta_12),nrow = 17+8), Hetero_f2)

log_delta_h1_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[5,]),nrow = 17),Hetero_f1)
log_delta_h1_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[5,],beta_12),nrow = 17+8), Hetero_f1)
log_delta_h2_4 = Find_true_delta_new(1000, 0.5, 0,matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[5,]),nrow = 17),Hetero_f2)
log_delta_h2_12 = Find_true_delta_new(1000, 0.5, 8,
                                       beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[5,],beta_12),nrow = 17+8), Hetero_f2)














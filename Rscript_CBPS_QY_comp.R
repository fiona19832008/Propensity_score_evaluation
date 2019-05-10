main = function() {
  
  source("New_simulation.R")
  source("Evaluation_function.R")
  source("Logistic_CBPS.R")
  source("Qingyuan_tailored_loss.R")
  source("Tailor_loss_kernel_CV_function.R")
  source("Compare_function.R")
  source("Kernel_balance.R")
  # Run the script as follows:
  # Rscript 
  # n: number of cell sampled
  # seed: seed used for random sampling
  # Fname: file name to be saved
  
  args = commandArgs(T)
  iter = args[1]
  n = args[2]
  K = args[3]
  delta = args[4]
  h = args[5]
  n_sigma = args[6]
  r = args[7]
  seed = args[8]
  QY_R_fold = args[9]
  saveName = args[10]
  m = args[11]
  a = args[12]
  Hetero = args[13]
  overlap = args[14]
  
  R_list = list.files(QY_R_fold, full.names = T)
  for(i in 1:length(R_list)) {
    source(R_list[i])
  }
  
  iter = as.numeric(iter)
  n = as.numeric(n)
  K = as.numeric(K)
  d = as.numeric(delta)
  h = as.numeric(h)
  n_sigma = as.numeric(n_sigma)
  r = as.numeric(r)
  m = as.numeric(m)
  a = as.numeric(a)
  
  #set.seed(as.numeric(seed))
  set.seed(1234)
  seeds = sample(seq(1,100000), 1000)
  Cov_linear = linear_cov_generate(K)
  print(seeds[1:10])
  
  if(overlap == "good") {
  	beta_12 = rep(c(-0.3,0.5,0,0),2)
  	b = c(-1.5,0.5,-0.75,2,-0.5)
  } else if (overlap == "poor") {
  	beta_12 = rep(c(1,-1,0,0),2)
  	b = c(-1,0.4,-1.5,2,-1.5)
  }
  
  alpha_12 = rep(c(1,0,-1,0),2)
  
  beta_model = rbind(rep(0,12), c(rep(0,6),2,rep(0,5)), c(0,5,rep(0,10)), c(4,rep(0,5),2,2,0,0,0,0), c(rep(0,10),1,0))
  alpha_model = rbind(rep(0,12), c(2,rep(0,5),1,1,rep(0,4)))
  
  if(Hetero == "h1") {
    delta_function = function(e) {
      D = e + 1
      return(D)
    }
  } else if (Hetero == "h2") {
    delta_function = function(e) {
      D = e^2 + 2*e + 1
      return(D)
    }
  }
  if(K == 0) {
    CBPS_QY_fit = CBPS_Qingyuan_wrapper(iter = iter, n, K, delta = d, delta_function, seeds = seeds,
                                        beta = matrix(c(b,beta_model[m,]),nrow = 17),
                                        alpha = matrix(c(0.5,1,0.6,2.2,-1.2,alpha_model[a,]),nrow = 17),
                                        Cov_linear, h, c(4:7), 2, "ATE", r = r, lambda_test = c(0, -10, 22), 
                                        normalize = FALSE, Normal_kernel, c=c(0.05,0.1), n_sigma = n_sigma)
  } else {
    CBPS_QY_fit = CBPS_Qingyuan_wrapper(iter = iter, n, K, delta = d, delta_function, seeds = seeds,
                                        beta = matrix(c(-1.5,0.5,-0.75,2,-0.5,beta_model[m,],beta_12),nrow = (17+K)),
                                        alpha = matrix(c(0.5,1,0.6,2.2,-1.2,alpha_model[a,], alpha_12), nrow = (17+K)),
                                        Cov_linear, h, c(4:7,20:27), 2, "ATE", r = r, lambda_test = c(0, -10, 22), 
                                        normalize = FALSE, Normal_kernel, c=c(0.05,0.1), n_sigma = n_sigma)
  }
  save(CBPS_QY_fit, file = saveName)
}

main()

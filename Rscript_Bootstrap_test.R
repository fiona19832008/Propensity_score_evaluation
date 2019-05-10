main = function() {
  source("New_simulation.R")
  source("Bootstrap_var.R")
  source("Evaluation_function.R")
  source("Compare_function.R")
  
  args = commandArgs(T)
  iter = args[1]
  n = args[2]
  B = args[3]
  h = args[4]
  m = args[5]
  delta = args[6]
  outfile = args[7]
  
  iter = as.numeric(iter)
  n = as.numeric(n)
  B = as.numeric(B)
  h = as.numeric(h)
  m = as.numeric(m)
  delta = as.numeric(delta)
  
  beta_12 = c(1,-1,0,0,1,-1,0,0)
  Cov_linear = linear_cov_generate(8)
  beta_model = rbind(rep(0,10), c(0,5,rep(0,8)), c(2,rep(0,5),1,1,0,0), c(5,rep(0,9)))
  
  set.seed(123)
  seeds = sample(seq(1,1000000), 1000)
  
  model.fit = CBPS_boot_wrapper(iter, n, 8, delta = delta, beta = matrix(c(-1.5,0.5,-0.75,2,-0.5, beta_model[m,], beta_12),nrow = 23), 
                                h = h, delta_function_1, Cov_linear, c(4:7,18:25), 2, B, seeds)
  
  write.csv(model.fit, outfile)
}
main()
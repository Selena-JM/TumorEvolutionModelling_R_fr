#date:03/10/2024, author : Séléna Jean-Mactoux
#Optimization pb 2 : Parameter identifiability 

# ---- Computing resulting functions : solve ode with the estimated parameters ----
compute_odes = function(pat_nb, Data, bounds, plot=FALSE){
  p_min = bounds[1:6]
  p_max = bounds[7:12]
  
  parameters_opt = Data$parameters_opt[[pat_nb]]
  
  # Parameters
  x1 = parameters_opt[1]
  
  para_min = list(
    sigma = p_min[1],
    rho = p_min[2],
    eta = p_min[3],
    mu = p_min[4],
    delta = p_min[5],
    alpha = p_min[6])
  
  para_max = list(
    sigma = p_max[1],
    rho = p_max[2],
    eta = p_max[3],
    mu = p_max[4],
    delta = p_max[5],
    alpha = p_max[6])
  
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  T0 = 10^9
  
  #Initial conditions
  init = c("X"=x1, "Y"=TC_pat[1]/T0) #y1 = y(tau=tau1) with tau1 = k2*K*T0*t/100, I take y1 = TC_pat/T0
  
  # Normalised time
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  time = seq(from=min(time_pat), to=max(time_pat), by=0.01) #*(max(time_pat)-min(time_pat))
  
  #Solve differential equations
  ODE_sol_min = ode(y=init, times = time, func = derive, parms = para_min, method="bdf")
  y_min = ODE_sol_min[,3]
  
  ODE_sol_max = ode(y=init, times = time, func = derive, parms = para_max, method="bdf")
  y_max = ODE_sol_max[,3]
  
  ## Find the right value for the comparison in cost function
  # Interpolation
  y_est_min = spline(time, y_min, xout = time_pat)$y
  y_est_max = spline(time, y_max, xout = time_pat)$y
  
  if (plot == TRUE){
    sol_opti1(pat_nb, Data, parameters_opt, TRUE)
    lines(time, y_max, type = "l", col='blue')
    lines(time, y_min, type = "l", col='red')
  }
  res = list(y_min = y_min, y_max = y_max, y_est_min=y_est_min, y_est_max = y_est_max)
  return(res)
}

# ---- Optimization pb 2 ----
opti2_pat = function(pat_nb, Data){
  epsilon = 0.2
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  time = Data$time[[pat_nb]]
  y_opt = Data$y_opt[[pat_nb]]
  parameters_opt = Data$parameters_opt[[pat_nb]]
  
  f_minimize_OP2 = function(bounds) {
    #y_pat = data for the patient
    #temps_pat = times of data for the patient, normalised
    #para = parameters for the patient, [sigma, rho, eta, mu, delta, alpha]
    
    # Parameters
    x1 = parameters_opt[1]
    para <- list(
      sigma = parameters_opt[2],
      rho = parameters_opt[3],
      eta = parameters_opt[4],
      mu = parameters_opt[5],
      delta = parameters_opt[6],
      alpha = parameters_opt[7])
    
    T0 = 10^9
    
    p_min = bounds[1:6]
    p_max = bounds[7:12]
    
    #Cost function
    Cost = -sum(p_max-p_min) #-sum bc minimize the function instead of maximize, no abs() and sum in what they wrote but weird
    return(Cost)
  }
  
  eval_jac_g_ineq = function(bounds){
    jacobian(eval_g_ineq, bounds)
  }
  
  eval_g_ineq = function(bounds){
    odes = compute_odes(pat_nb, Data, bounds)
    
    y_est_min = odes$y_est_min
    y_est_max = odes$y_est_max
    y_est_opt = spline(time, y_opt, xout = time_pat)$y
    
    ineq_values = c(
      y_est_max - (1 + epsilon) * y_est_opt,   # y_est_max <= (1 + epsilon) * y_est_opt
      (1 - epsilon) * y_est_opt - y_est_max,   # y_est_max >= (1 - epsilon) * y_est_opt
      y_est_min - (1 + epsilon) * y_est_opt,   # y_est_min <= (1 + epsilon) * y_est_opt
      (1 - epsilon) * y_est_opt - y_est_min    # y_est_min >= (1 - epsilon) * y_est_opt
    )
    # ineq_values = c(
    #   -(epsilon - 1)*y_est_opt - y_est_max,  
    #   -(epsilon + 1)*y_est_opt + y_est_max,  
    #   -(epsilon - 1)*y_est_opt - y_est_min,  
    #   -(epsilon + 1)*y_est_opt + y_est_min   
    # )
    
    # print(paste("Inégalités:", paste(round(ineq_values, 6), collapse = ", ")))
    
    return(ineq_values)
  }
  
  # Computing of gradient
  gradient_OP2 <- function(bounds) {
    grad <- grad(func = f_minimize_OP2, x = bounds)
    print(paste("Gradient:", paste(round(grad, 6), collapse=", ")))  # Display the gradient
    return(grad)
  }
  
  # Limits for x1 and parameters
  lower_bounds <- c(rep(10^(-2),6), parameters_opt[2:7]) #lower bound of p_min is lower bound of all para : 10^-2 and lower bound of p_max is p_opt
  upper_bounds <- c(parameters_opt[2:7], rep(c(10^2), 6)) 
  
  #Starting values for the parameter optimization
  # start_para = c("p_min" = c("sigma" = 10^(-2), "rho" = 10^(-2), "eta" = 10^(-2), "mu" = 10^(-2), "delta" = 10^(-2), "alpha" = 10^(-2)), "p_max" = c("sigma" = 10^2, "rho" = 10^2, "eta" = 10^2, "mu" = 10^2, "delta" = 10^2, "alpha" = 10^2))
  # start_para = c("p_min" = pmin(parameters_opt[2:7], rep(10^(-2), 6)), "p_max" = pmax(parameters_opt[2:7], rep(10^2, 6)))
  # start_para = c("p_min" = rep(10^(-2), 6), "p_max" = rep(10^2, 6))
  # start_para = c("p_min" = pmax(parameters_opt[2:7]-10, rep(10^(-2), 6)), "p_max" = pmin(parameters_opt[2:7]+10, rep(10^2, 6)))
  # start_para = c("p_min" = c("sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 10^(-2)), "p_max" = c("sigma" = 10^2, "rho" = 10^2, "eta" = 10^2, "mu" = 10^2, "delta" = 10^2, "alpha" = 1))
  # start_para = c("p_min" = runif(6, max=parameters_opt[2:7], min=rep(10^(-2),6)), "p_max" = runif(6, min=parameters_opt[2:7], max=rep(10^2,6)))
  start_para = c("p_min" = upper_bounds[1:6], "p_max" = lower_bounds[7:12])
  
  # Optimization
  result <- nloptr(
    x0 = start_para,
    eval_f = f_minimize_OP2,
    eval_grad_f = gradient_OP2, # Using the defined gradient
    eval_g_ineq = eval_g_ineq,
    eval_jac_g_ineq = eval_jac_g_ineq,
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-8, "maxeval"=500, "print_level"=0) #NLOPT_LD_MMA, NLOPT_LD_SLSQP
  )
  
  # print(result)
  return(result)
}


# ---- Building data base for opti pb 2 ----
add_op2_Data = function(Data_post_op1){
  Data_ter = Data_post_op1
  
  nb_pat = length(Data_ter$Patient_Anonmyized)
  
  Data_ter$y_min = rep(NA,nb_pat)
  Data_ter$y_max = rep(NA,nb_pat)
  Data_ter$parameters_min = rep(NA,nb_pat)
  Data_ter$parameters_max = rep(NA,nb_pat)
  
  for (i in 1:50){
    print("")
    print(i)
    # if (i!= 4 & i!=11 & i!=14 & i!=25){
      result = opti2_pat(i, Data_ter)
      
      bounds = result$solution
      parameters_min = bounds[1:6]
      parameters_max = bounds[7:12]
      
      odes = compute_odes(i, Data_ter, bounds)
      
      Data_ter$y_min[i] = I(list(odes$y_min))
      Data_ter$y_max[i] = I(list(odes$y_max))
      Data_ter$parameters_min[i] = I(list(parameters_min))
      Data_ter$parameters_max[i] = I(list(parameters_max))
      save(Data_ter, file="./Data_processed/Data_ter.Rda")
    # }
  }
  save(Data_ter, file="./Data_processed/Data_ter.Rda")
  return(Data_ter)
}

# ---- Plotting results ----
plot_op1_op2 = function(pat_nb, Data_ter){
  par(mar = c(5, 4, 4, 5.5))
  plot(Data_ter$time[[pat_nb]], Data_ter$y_opt[[pat_nb]], type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP2 results for patient :", pat_nb))
  lines(Data_ter$time[[pat_nb]], Data_ter$y_min[[pat_nb]], type = 'l', col = 'red')
  lines(Data_ter$time[[pat_nb]], Data_ter$y_max[[pat_nb]], type = 'l', col = 'blue')
  
  time_pat = Data_ter$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  points(time_pat, Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  legend("topright", 
         legend = c("y_opt", "y_min", "y_max", "observations"), 
         col = c("black", "red", "blue", "black"), 
         lty = c(1, 1, 1, NA),   # lty=NA pour les points
         pch = c(NA, NA, NA, 1), # pch=16 pour les points des mesures
         cex = 0.6, 
         xpd=TRUE, 
         inset = c(-0.25, 0)) 
}

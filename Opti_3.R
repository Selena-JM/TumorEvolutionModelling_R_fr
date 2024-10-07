#date:03/10/2024, author : Séléna Jean-Mactoux
#Predictions & OP3: removing some points at the end and running OP1, then finding worst case scenarios for predictions

# ---- Building data base for predictions ----
add_pred_Data = function(Data, nb_points_omitted){
  Data_qua = Data
  
  nb_pat = length(Data_qua$Patient_Anonmyized)
  
  Data_qua$y_pred = rep(NA,nb_pat)
  Data_qua$parameters_pred = rep(NA,nb_pat)
  Data_qua$nb_points_omitted = rep(NA,nb_pat)
  
  for (i in 1:20){ #pat14
    print("")
    print(i)
    result = opti1_pat(i, Data_qua, nb_points_omitted = nb_points_omitted)
    parameters = result$solution
    
    sol = sol_opti1(i, Data_qua, parameters)
    
    Data_qua$y_pred[i] = I(list(sol$y))
    Data_qua$parameters_pred[i] = I(list(parameters))
    Data_qua$nb_points_omitted[i] = nb_points_omitted
    save(Data_qua, file="./Data_processed/Data_qua.Rda")
  }
  save(Data_qua, file="./Data_processed/Data_qua.Rda")
  return(Data_qua)
}

# ---- Computing resulting functions : solve ode with the estimated parameters ----
compute_odes_op3 = function(pat_nb, Data, bounds, plot=FALSE){
  x1_low = bounds[1]
  p_min = bounds[2:7]
  x1_upp = bounds[7]
  p_max = bounds[9:14]
  
  # Parameters
  p_low = list(
    sigma = p_min[1],
    rho = p_min[2],
    eta = p_min[3],
    mu = p_min[4],
    delta = p_min[5],
    alpha = p_min[6])
  
  p_upp = list(
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
  init_low = c("X"=x1_low, "Y"=TC_pat[1]/T0) #y1 = y(tau=tau1) with tau1 = k2*K*T0*t/100, I take y1 = TC_pat/T0
  init_upp = c("X"=x1_upp, "Y"=TC_pat[1]/T0) 
  
  # Normalised time
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  time = seq(from=min(time_pat), to=max(time_pat), by=0.01) #*(max(time_pat)-min(time_pat))
  
  #Solve differential equations
  ODE_sol_low = ode(y=init_low, times = time, func = derive, parms = p_low, method="bdf")
  y_low = ODE_sol_low[,3]
  
  ODE_sol_upp = ode(y=init_upp, times = time, func = derive, parms = p_upp, method="bdf")
  y_upp = ODE_sol_upp[,3]
  
  ## Find the right value for the comparison in cost function
  # Interpolation
  y_est_low = spline(time, y_low, xout = time_pat)$y
  y_est_upp = spline(time, y_upp, xout = time_pat)$y
  
  if (plot == TRUE){
    plot(time, Data$y_pred[[pat_nb]], type = 'l', col = 'black')
    lines(time, y_upp, type = "l", col='blue')
    lines(time, y_low, type = "l", col='red')
    point(time_pat, TC_pat/T0)
  }
  
  return(list(y_low = y_low, y_upp = y_upp, y_est_low=y_est_low, y_est_upp = y_est_upp))
}

# ---- Optimization pb 3 ----
opti3_pat = function(pat_nb, Data){
  theta = 0.1
  T0 = 10^9
  
  
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  time = Data$time[[pat_nb]]
  y_pred = Data$y_pred[[pat_nb]]
  parameters_pred = Data$parameters_pred[[pat_nb]]
  nb_point_omitted = Data$nb_points_omitted[pat_nb]
  
  parameters_opt = Data$parameters_opt[[pat_nb]]
  
  f_minimize_OP3 = function(bounds) {
    # Parameters
    x1_low = bounds[1]
    p_low = list(
      sigma = bounds[2],
      rho = bounds[3],
      eta = bounds[4],
      mu = bounds[5],
      delta = bounds[6],
      alpha = bounds[7])
    
    x1_upp = bounds[8]
    p_upp = list(
      sigma = bounds[9],
      rho = bounds[10],
      eta = bounds[11],
      mu = bounds[12],
      delta = bounds[13],
      alpha = bounds[14])
    
    #Initial conditions
    init_low = c("X"=x1_low, "Y"=TC_pat[1]/T0) #y1 = y(tau=tau1) with tau1 = k2*K*T0*t/100, I take y1 = TC_pat/T0
    init_upp = c("X"=x1_upp, "Y"=TC_pat[1]/T0)
    
    #Solve differential equations
    ODE_sol_low = ode(y=init_low, times = time, func = derive, parms = p_low, method="bdf")
    y_low = ODE_sol_low[,3]
    
    ODE_sol_upp = ode(y=init_upp, times = time, func = derive, parms = p_upp, method="bdf")
    y_upp = ODE_sol_upp[,3]
    
    # Interpolation
    y_est_low <- spline(time, y_low, xout = time_pat)$y
    y_est_upp <- spline(time, y_upp, xout = time_pat)$y
    
    #Cost function
    Cost = -y_est_upp[length(y_est_upp)] + y_est_low[length(y_est_low)]
    return(Cost)
  }
  
  eval_jac_g_ineq = function(bounds){
    jacobian(eval_g_ineq, bounds)
  }
  
  eval_g_ineq = function(bounds){
    odes = compute_odes_op3(pat_nb, Data, bounds)
    
    y_est_low = spline(time, odes$y_low, xout = time_pat[1:(length(time_pat)-nb_point_omitted)])$y
    y_est_upp = spline(time, odes$y_upp, xout = time_pat[1:(length(time_pat)-nb_point_omitted)])$y
    y_est_pred = spline(time, y_pred, xout = time_pat[1:(length(time_pat)-nb_point_omitted)])$y
    
    ineq_values = c(
      y_est_upp - (1 + theta) * y_est_pred,   # y_est_max <= (1 + epsilon) * y_est_opt
      (1 - theta) * y_est_pred - y_est_upp,   # y_est_max >= (1 - epsilon) * y_est_opt
      y_est_low - (1 + theta) * y_est_pred,   # y_est_min <= (1 + epsilon) * y_est_opt
      (1 - theta) * y_est_pred - y_est_low    # y_est_min >= (1 - epsilon) * y_est_opt
    )
    
    return(ineq_values)
  }
  
  # Computing of gradient
  gradient_OP3 <- function(bounds) {
    grad <- grad(func = f_minimize_OP3, x = bounds)
    print(paste("Gradient:", paste(round(grad, 6), collapse=", ")))  # Display the gradient
    return(grad)
  }
  
  # Limits for x1 and parameters
  lower_bounds <- rep(10^(-2),14) #lower bound of p_min is lower bound of all para : 10^-2 and lower bound of p_max is p_opt
  upper_bounds <- rep(c(10^2), 14) 
  
  #Starting values for the parameter optimization
  start_para = c("x1_low" = 1, "p_low" = c("sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 1), "x1_upp" = 1, "p_upp" = c("sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 1))
  # start_para = c("x1_low" = parameters_pred[1], "p_min" = parameters_pred[2:7], "x1_upp" = parameters_pred[8], "p_max" = parameters_pred[9:14])
  # start_para = c("x1_low" = parameters_opt[1], "p_min" = parameters_opt[2:7], "x1_upp" = parameters_opt[8], "p_max" = parameters_opt[9:14])
  
  # Optimization
  result <- nloptr(
    x0 = start_para,
    eval_f = f_minimize_OP3,
    eval_grad_f = gradient_OP3, # Using the defined gradient
    eval_g_ineq = eval_g_ineq,
    eval_jac_g_ineq = eval_jac_g_ineq,
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-8, "maxeval"=500, "print_level"=0) #NLOPT_LD_MMA, NLOPT_LD_SLSQP
  )
  
  # print(result)
  return(result)
}

# ---- Building data base for opti pb 3 ----
add_op3_Data = function(Data){
  Data_5 = Data
  
  nb_pat = length(Data_5$Patient_Anonmyized)
  
  Data_5$y_low = rep(NA,nb_pat)
  Data_5$y_upp = rep(NA,nb_pat)
  Data_5$parameters_low = rep(NA,nb_pat)
  Data_5$parameters_upp = rep(NA,nb_pat)
  
  for (i in 1:13){
    print("")
    print(i)
    # if (i!= 4 & i!=11 & i!=14 & i!=25){
      result = opti3_pat(i, Data_5)

      bounds = result$solution
      parameters_low = bounds[1:7]
      parameters_upp = bounds[8:14]
      # parameters_low = Data_5$parameters_low[[i]]
      # parameters_upp = Data_5$parameters_low[[i]]
      
      odes = compute_odes_op3(i, Data, bounds, plot=FALSE)
      
      Data_5$y_low[i] = I(list(odes$y_low))
      Data_5$y_upp[i] = I(list(odes$y_upp))
      Data_5$parameters_low[i] = I(list(parameters_low))
      Data_5$parameters_upp[i] = I(list(parameters_upp))
      save(Data_5, file="./Data_processed/Data_5.Rda")
    # }
  }
  save(Data_5, file="./Data_processed/Data_5.Rda")
  return(Data_5)
}

# ---- Plotting results ----
# Plot predictions
plot_pred = function(pat_nb, Data){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-Data$nb_points_omitted[pat_nb]+1]-0.01
  
  ymax = max(max(Data$y_opt[[pat_nb]]), max(Data$y_pred[[pat_nb]]), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(Data$y_opt[[pat_nb]]), min(Data$y_pred[[pat_nb]]), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  # plot environment
  par(mar = c(5, 5, 3, 6))
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'n', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Prediction results for patient :", pat_nb))
  
  # Rectangle
  rect(xleft = x_col, ybottom = par("usr")[3], xright = par("usr")[2], 
       ytop = par("usr")[4], col = rgb(0.9, 0.8, 0.8, 0.5), border = NA) 
  
  
  lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'l', col = 'black')
  lines(Data$time[[pat_nb]], Data$y_pred[[pat_nb]], type = 'l', col = 'red')
  
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  legend("topright", 
         legend = c("y_opt", "y_pred", "observations"), 
         col = c("black", "red", "black"), 
         lty = c(1, 1, NA),   # lty=NA pour les points
         pch = c(NA, NA, 1), # pch=16 pour les points des mesures
         cex = 0.6, 
         xpd=TRUE, 
         inset = c(-0.25, 0))
}

# plot predictions and op3
plot_pred_op3 = function(pat_nb, Data){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-Data$nb_points_omitted[pat_nb]+1]-0.01
  
  ymax = max(max(Data$y_low[[pat_nb]]), max(Data$y_upp[[pat_nb]]), max(Data$y_pred[[pat_nb]]), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymax = min(min(Data$y_low[[pat_nb]]), min(Data$y_upp[[pat_nb]]), min(Data$y_pred[[pat_nb]]), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  
  # plot environment
  par(mar = c(5, 5, 3, 6))
  plot(Data$time[[pat_nb]], Data$y_pred[[pat_nb]], type = 'n', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP3 results for patient :", pat_nb))
  
  
  # Rectangle
  rect(xleft = x_col, ybottom = par("usr")[3], xright = par("usr")[2], 
       ytop = par("usr")[4], col = rgb(0.9, 0.8, 0.8, 0.5), border = NA)
  
  # Plot 
  lines(Data$time[[pat_nb]], Data$y_pred[[pat_nb]], col = 'black')
  lines(Data$time[[pat_nb]], Data$y_low[[pat_nb]], col = 'red')
  lines(Data$time[[pat_nb]], Data$y_upp[[pat_nb]], col = 'blue')
  
  # Observations
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9, col = 'black')
  
  # legend
  legend("topright", 
         legend = c("y_pred", "y_low", "y_upp", "observations"), 
         col = c("black", "red", "blue", "black"), 
         lty = c(1, 1, 1, NA),   # lty=NA for observations
         pch = c(NA, NA, NA, 1), # pch=1 for observations
         cex = 0.6, 
         xpd=TRUE, 
         inset = c(-0.25, 0))
}

# plot predictions without data base, for testing
plot_pred_test = function(pat_nb, Data_ter, nb_points_omitted=2, maxeval=500, precision=10^(-8)){
  result = opti1_pat(pat_nb, Data_ter, nb_points_omitted = nb_points_omitted, maxeval=maxeval, precision=precision)
  parameters = result$solution

  sol = sol_opti1(pat_nb, Data_ter, parameters)

  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-nb_points_omitted+1]-0.01
  
  ymax = max(max(Data_ter$y_opt[[pat_nb]]), max(sol$y), max(Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(Data_ter$y_opt[[pat_nb]]), min(sol$y), min(Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  #plot environment
  par(mar = c(5, 5, 3, 6))
  plot(sol$time, Data_ter$y_opt[[pat_nb]], type = 'n', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Prediction results for patient :", pat_nb), ylim=c(ymin,ymax))
  
  # Rectangle
  rect(xleft = x_col, ybottom = par("usr")[3], xright = par("usr")[2], 
       ytop = par("usr")[4], col = rgb(0.9, 0.8, 0.8, 0.5), border = NA)
  
  #plot results
  lines(sol$time, Data_ter$y_opt[[pat_nb]], type = 'l', col = 'black')
  lines(sol$time, sol$y, type = 'l', col = 'red')

  #observations
  points(time_pat, Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9)

  legend("topright",
         legend = c("y_opt", "y_pred", "observations"),
         col = c("black", "red", "black"),
         lty = c(1, 1, NA),   # lty=NA pour les points
         pch = c(NA, NA, 1), # pch=16 pour les points des mesures
         cex = 0.6,
         xpd=TRUE,
         inset = c(-0.25, 0))
  return(sol)
}

#plot op3 without data bases, for tests
plot_pred_op3_test = function(pat_nb, Data_qua, sol){  
  
  result = opti3_pat(pat_nb, Data_qua)
  
  bounds = result$solution
  
  odes = compute_odes_op3(pat_nb, Data_qua, bounds, plot=FALSE)
  y_pred = sol_pred$y
  time = Data_qua$time[[pat_nb]]
  y_low = odes$y_low
  y_upp = odes$y_upp
  
  plot(time, y_pred, type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP3 results for patient :", pat_nb))
  lines(time, y_low, type = 'l', col = 'red')
  lines(time, y_upp, type = 'l', col = 'blue')
  
  time_pat = Data_qua$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  points(time_pat, Data_qua$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  legend("topright",
         legend = c("y_pred", "y_low", "y_upp", "observations"),
         col = c("black", "red", "blue", "black"),
         lty = c(1, 1, 1, NA),   # lty=NA pour les points
         pch = c(NA, NA, NA, 1), # pch=16 pour les points des mesures
         cex = 0.6,
         xpd=TRUE,
         inset = c(-0.25, 0))

}
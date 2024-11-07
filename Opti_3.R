#date:03/10/2024, author : Séléna Jean-Mactoux
#OP3: finding worst case scenarios for predictions

source("derive.R")

# ---- Optimization pb 3 ----
opti3_pat = function(pat_nb, Data, TC_pat=NA, y_pred=NA){
  theta = 0.1
  T0 = 10^9
  
  
  if (length(TC_pat)==1){
    TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  }
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  time = Data$time[[pat_nb]]
  
  if (length(y_pred)==1){
    y_pred = Data$y_pred[[pat_nb]]
  }
  
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
    y_est_low = spline(time, y_low, xout = time_pat)$y
    y_est_upp = spline(time, y_upp, xout = time_pat)$y
    
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
    # print(paste("Gradient:", paste(round(grad, 6), collapse=", ")))  # Display the gradient
    return(grad)
  }
  
  # Limits for x1 and parameters
  lower_bounds <- rep(10^(-2),14) #lower bound of p_min is lower bound of all para : 10^-2 and lower bound of p_max is p_opt
  upper_bounds <- rep(c(10^2), 14) 
  
  #Starting values for the parameter optimization
  # start_para = c("x1_low" = 1, "p_low" = c("sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 1), "x1_upp" = 1, "p_upp" = c("sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 1))
  start_para = c("x1_low" = parameters_pred[1], "p_min" = parameters_pred[2:7], "x1_upp" = parameters_pred[1], "p_max" = parameters_pred[2:7])

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
  Data_OP3 = Data
  
  nb_pat = length(Data_OP3$Patient_Anonmyized)

  Data_OP3$y_low = rep(NA,nb_pat)
  Data_OP3$y_upp = rep(NA,nb_pat)
  Data_OP3$parameters_low = rep(NA,nb_pat)
  Data_OP3$parameters_upp = rep(NA,nb_pat)
  
  for (i in 1:nb_pat){
    print(paste("Patient",i))
    if (i!= 14 & i!= 60 & i!= 104 & i!= 168 & i!=189){
      result = opti3_pat(i, Data_OP3)

      bounds = result$solution
      parameters_low = bounds[1:7]
      parameters_upp = bounds[8:14]
      
      odes = compute_odes_op3(i, Data, bounds, plot=FALSE)
      
      Data_OP3$y_low[i] = I(list(odes$y_low))
      Data_OP3$y_upp[i] = I(list(odes$y_upp))
      Data_OP3$parameters_low[i] = I(list(parameters_low))
      Data_OP3$parameters_upp[i] = I(list(parameters_upp))
      save(Data_OP3, file="./Data_processed/Data_OP3.Rda")
    }
  }
  save(Data_OP3, file="./Data_processed/Data_OP3.Rda")
  return(Data_OP3)
}

# ---- Analysis ----
goodness_prediction_intervals = function(pat_nb, Data){
  if (length(which(is.na(Data$y_pred[[pat_nb]])))==0 &length(which(is.na(Data$y_low[[pat_nb]])))==0){
    T0 = 10^9
    TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
    
    y_obs = TC_pat/T0
    
    time_pat = Data$Treatment_Day[[pat_nb]]
    time_pat = time_pat/(max(time_pat) - min(time_pat)) 
    
    time = Data$time[[pat_nb]]
    
    y_low = Data$y_low[[pat_nb]]
    y_upp = Data$y_upp[[pat_nb]]
    y_pred = Data$pred[[pat_nb]]
    
    y_low_est = spline(time, y_low, xout = time_pat)$y
    y_upp_est = spline(time, y_upp, xout = time_pat)$y
    y_pred_est = spline(time, y_pred, xout = time_pat)$y
    
    # Only compare the results for the predicted points
    nb_points_omitted = Data$nb_points_omitted[pat_nb]

    y_obs = y_obs[(length(y_obs)-nb_points_omitted +1):length(y_obs)]
    y_low_est = y_low_est[(length(y_low_est)-nb_points_omitted +1):length(y_low_est)]
    y_upp_est = y_upp_est[(length(y_upp_est)-nb_points_omitted +1):length(y_upp_est)]
    y_pred_est = y_pred_est[(length(y_pred_est)-nb_points_omitted +1):length(y_pred_est)]
    
    nb_points_in_interval = length(which(y_obs<=y_upp_est & y_obs>=y_low_est))

    # Find the closest to the data points between y_low, y_pred and y_upp and compute the normal analysis, 
    diff_low = abs(y_low_est - y_obs)
    diff_upp = abs(y_upp_est - y_obs)
    diff_pred = abs(y_pred_est - y_obs)

    y_est = rep(NA, nb_points_omitted)
    for (k in 1:nb_points_omitted){
      y_est[k] = min(diff_low[k], diff_upp[k], diff_pred[k])
    }
    
    
    # Metrics
    # R²
    ss_total = sum((y_obs - mean(y_est))^2)
    ss_residual = sum((y_obs - y_est)^2)
    r_squared = 1 - (ss_residual / ss_total)
    
    # MAE
    mae = mean(abs(y_obs - y_est))
    
    # MSE
    mse = mean((y_obs - y_est)^2)
  }
  else {
    r_squared = NA
    mae = NA
    mse = NA
    nb_points_in_interval = NA
  }
  return(list(MAE = mae,MSE = mse, R2 = r_squared, nb_points_in_interval=nb_points_in_interval))
}


goodness_prediction_intervals_analysis = function(Data){
  
  Fit_OP3 = data.frame(MAE = numeric(), MSE = numeric(), R2 = numeric(), l_nb_in_intervals=numeric())
  
  for (i in 1:length(Data$Patient_Anonmyized)){
    Fit_OP3_i = goodness_prediction_intervals(i, Data)
    
    Fit_OP3 = rbind(Fit_OP3, Fit_OP3_i)
  }  
  
  save(Fit_OP3, file = "./Data_processed/Fit_OP3.Rda")
  return(Fit_OP3)
}


# ---- Plotting results ----

# plot predictions and op3
plot_pred_op3 = function(pat_nb, Data){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-Data$nb_points_omitted[pat_nb]+1]-0.01
  
  ymax = max(max(Data$y_low[[pat_nb]]), max(Data$y_pred[[pat_nb]]), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9), max(Data$y_opt[[pat_nb]]), max(Data$y_upp[[pat_nb]]))
  ymin = min(min(Data$y_low[[pat_nb]]), min(Data$y_upp[[pat_nb]]), min(Data$y_pred[[pat_nb]]), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9), min(Data$y_opt[[pat_nb]]))
  
  
  # plot environment
  par(mar = c(5, 4, 4, 5), mfrow=c(1,1))
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'n', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP3 results for patient :", pat_nb), ylim=c(ymin, ymax))
  
  
  # Rectangle
  rect(xleft = x_col, ybottom = par("usr")[3], xright = par("usr")[2], 
       ytop = par("usr")[4], col = rgb(0.9, 0.8, 0.8, 0.5), border = NA)
  
  # Plot 
  lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], col = 'black')
  lines(Data$time[[pat_nb]], Data$y_pred[[pat_nb]], col = 'black', lty=4)
  lines(Data$time[[pat_nb]], Data$y_low[[pat_nb]], col = 'red')
  lines(Data$time[[pat_nb]], Data$y_upp[[pat_nb]], col = 'blue')
  
  # Observations
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9, col = 'black')
  
  # legend
  legend("topright", 
         legend = c("y_opt", "y_pred", "y_low", "y_upp", "observations"), 
         col = c("black","black", "red", "blue", "black"), 
         lty = c(1, 4, 1, 1, NA),   # lty=NA for observations
         pch = c(NA, NA, NA, NA, 1), # pch=1 for observations
         cex = 0.5, 
         xpd=TRUE, 
         inset = c(-0.25, 0))
}


#plot op3 without data bases, for tests
plot_pred_op3_test = function(pat_nb, Data_Pred){  
  
  result = opti3_pat(pat_nb, Data_Pred)
  
  bounds = result$solution
  
  odes = compute_odes_op3(pat_nb, Data_Pred, bounds, plot=FALSE)
  y_pred = Data_Pred$y_pred[[pat_nb]]
  time = Data_Pred$time[[pat_nb]]
  y_low = odes$y_low
  y_upp = odes$y_upp
  
  ymax = max(max(y_low), max(y_upp), max(y_pred),max(Data_Pred$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(y_low), min(y_upp), min(y_pred),min(Data_Pred$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  plot(time, y_pred, type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP3 results for patient :", pat_nb), ylim=c(ymin,ymax))
  lines(time, y_low, type = 'l', col = 'red')
  lines(time, y_upp, type = 'l', col = 'blue')
  
  time_pat = Data_Pred$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  points(time_pat, Data_Pred$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  legend("topright",
         legend = c("y_pred", "y_low", "y_upp", "observations"),
         col = c("black", "red", "blue", "black"),
         lty = c(1, 1, 1, NA),   # lty=NA pour les points
         pch = c(NA, NA, NA, 1), # pch=16 pour les points des mesures
         cex = 0.6,
         xpd=TRUE,
         inset = c(-0.25, 0))
}
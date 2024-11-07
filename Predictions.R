#date:03/10/2024, author : Séléna Jean-Mactoux
#Predictions: removing some points at the end and running OP1, then finding worst case scenarios for predictions

source("derive.R")

# ---- Building data base for predictions ----
add_pred_Data = function(Data, nb_points_omitted){
  Data_Pred = Data
  
  nb_pat = length(Data_Pred$Patient_Anonmyized)
  
  Data_Pred$y_pred = rep(NA,nb_pat)
  Data_Pred$parameters_pred = rep(NA,nb_pat)
  Data_Pred$nb_points_omitted = rep(NA,nb_pat)
  
  for (i in 1:nb_pat){
    print(paste("Patient", i))
    if (i!=14){
      result = opti1_pat(i, Data_Pred, nb_points_omitted = nb_points_omitted)
      parameters = result$solution
      
      sol = sol_opti1(i, Data_Pred, parameters)
      
      Data_Pred$y_pred[i] = I(list(sol$y))
      Data_Pred$parameters_pred[i] = I(list(parameters))
      Data_Pred$nb_points_omitted[i] = nb_points_omitted
      save(Data_Pred, file="./Data_processed/Data_Pred.Rda")
    }
  }
  save(Data_Pred, file="./Data_processed/Data_Pred.Rda")
  return(Data_Pred)
}

# ---- Computing resulting functions : solve ode with the estimated parameters ----
compute_odes_op3 = function(pat_nb, Data, bounds, plot=FALSE){
  x1_low = bounds[1]
  p_min = bounds[2:7]
  x1_upp = bounds[8]
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

# ---- Analysis ----
goodness_predictions = function(pat_nb, Data){
  T0 = 10^9
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  
  y_obs = TC_pat/T0
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat)) 
  
  time = Data$time[[pat_nb]]
  
  y = Data$y_pred[[pat_nb]]
  
  y_est = spline(time, y, xout = time_pat)$y
  
  # Only compare the results for the predicted points
  nb_points_omitted = Data$nb_points_omitted[pat_nb]
  y_obs = y_obs[(length(y_obs)-nb_points_omitted +1):length(y_obs)]
  y_est = y_est[(length(y_est)-nb_points_omitted +1):length(y_est)]
  
  # Metrics
  # R²
  ss_total = sum((y_obs - mean(y_est))^2)
  ss_residual = sum((y_obs - y_est)^2)
  r_squared = 1 - (ss_residual / ss_total)
  
  # MAE
  mae = mean(abs(y_obs - y_est))
  
  # MSE
  mse = mean((y_obs - y_est)^2)
  return(list(MAE = mae,MSE = mse, R2 = r_squared))
}


goodness_predictions_analysis = function(Data){
  
  parameters_df = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_1 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_2 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_3 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_4 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_5 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  
  GF_df = data.frame(MAE = numeric(), MSE = numeric(), R2 = numeric())
  
  
  for (i in 1:length(Data$Patient_Anonmyized)){
    if (length(which(is.na(Data$y_pred[[i]])))==0){
      GF_i = goodness_predictions(i, Data)
      parameters_i = data.frame(x1=Data$parameters_pred[[i]][1], sigma=Data$parameters_pred[[i]][2], rho=Data$parameters_pred[[i]][3], eta=Data$parameters_pred[[i]][4], mu=Data$parameters_pred[[i]][5], delta=Data$parameters_pred[[i]][6], alpha=Data$parameters_pred[[i]][7])
  
      GF_df = rbind(GF_df, GF_i)
      parameters_df = rbind(parameters_df, parameters_i)
  
      j = substr(Data$Study_Arm[[i]], 7, 7)
      df = get(paste("parameters_df_", j, sep=""))
      
      df = rbind(df, parameters_i)
      
      assign(paste("parameters_df_", j, sep=""), df)
    }
  }  
  
  Fit_pred = list("GF" = GF_df, "Para_all" = parameters_df, "Para_1" = parameters_df_1, "Para_2" = parameters_df_2, "Para_3" = parameters_df_3, "Para_4" = parameters_df_4,"Para_5" = parameters_df_5)
  save(Fit_pred, file = "./Data_processed/Fit_pred.Rda")
  return(Fit_pred)
}


# ---- Plotting results ----
plot_pred = function(pat_nb, Data){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-Data$nb_points_omitted[pat_nb]+1]-0.01
  
  ymax = max(max(Data$y_opt[[pat_nb]]), max(Data$y_pred[[pat_nb]]), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(Data$y_opt[[pat_nb]]), min(Data$y_pred[[pat_nb]]), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  # plot environment
  par(mar = c(5, 4, 4, 5), mfrow=c(1,1))
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'n', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Prediction results for patient :", pat_nb), ylim=c(ymin, ymax))
  
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
         cex = 0.5, 
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
  par(mar = c(5, 4, 4, 5), mfrow=c(1,1))
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
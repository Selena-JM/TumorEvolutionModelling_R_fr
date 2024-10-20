#date:03/10/2024, author : Séléna Jean-Mactoux
#Optimization pb 1 : Finding optimal parameters for each patient
source("derive.R")

# ---- Optimization pb 1 ----
opti1_pat = function(pat_nb, Data, nb_points_omitted=0, maxeval=500, precision = 10^(-8)){
  n = length(Data$TargetLesionLongDiam_mm[[pat_nb]])
  
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]][1:(n-nb_points_omitted)]
  time_pat = Data$Treatment_Day[[pat_nb]][1:(n-nb_points_omitted)]
  T0 = 10^9
  
  f_minimize_OP1 = function(parameters) {
    #y_pat = data for the patient
    #temps_pat = times of data for the patient, normalised
    #para = parameters for the patient, [sigma, rho, eta, mu, delta, alpha]
    
    # Parameters
    x1 = parameters[1]
    para <- list(#x1 = parameters[1],
      sigma = parameters[2],
      rho = parameters[3],
      eta = parameters[4],
      mu = parameters[5],
      delta = parameters[6],
      alpha = parameters[7])
    
    #Initial conditions
    init = c("X"=x1, "Y"=TC_pat[1]/T0) #y1 = y(tau=tau1) with tau1 = k2*K*T0*t/100, I take y1 = TC_pat/T0
    
    # Normalised time
    time_pat = time_pat / (max(Data$Treatment_Day[[pat_nb]]) - min(Data$Treatment_Day[[pat_nb]]))
    time = seq(from=min(time_pat), to=max(time_pat), by=0.01) #*(max(time_pat)-min(time_pat))
    
    
    #Solve differential equations
    ODE_sol = ode(y=init, times = time, func = derive, parms = para, method="bdf")
    y = ODE_sol[,3]
    
    ## Find the right value for the comparison in cost function
    
    # Finding the closest times of time_pat in time
    # y_est = c()
    # index = 0
    # for (i in 1:length(time_pat)){
    #   index = which.min(abs(time - time_pat[i]))
    #   y_est = c(y_est, y[index])
    # }
    
    # Interpolation
    y_est = spline(time, y, xout = time_pat)$y
    
    #Cost function
    Cost = sum(abs(y_est-TC_pat/T0)^2) #no abs(), and not ^2 in what they wrote but weird
    return(Cost)
  }
  
  
  # Computing of gradient
  gradient_OP1 <- function(parameters) {
    grad <- grad(func = f_minimize_OP1, x = parameters)
    # print(paste("Gradient:", paste(round(grad, 6), collapse=", ")))  # Display the gradient
    return(grad)
  }
  
  # Limits for x1 and parameters
  lower_bounds <- c(10^(-2), rep(c(10^(-2)), 6))
  upper_bounds <- c(10^2, rep(c(10^2), 6)) 
  
  #Starting values for the parameter optimization
  start_para = c("x1" = 1, "sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 1)
  # start_para = c("x1" = 10, "sigma" = 0.1, "rho" = 1, "eta" = 0.01, "mu" = 1, "delta" = 0.1, "alpha" = 0.01)
  # start_para = c("x1" = 100, "sigma" = 100, "rho" = 100, "eta" = 1, "mu" = 10, "delta" = 10, "alpha" = 1)
  
  # Optimization
  result <- nloptr(
    x0 = start_para,
    eval_f = f_minimize_OP1,
    eval_grad_f = gradient_OP1, # Using the defined gradient
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = precision, "maxeval"=maxeval, print_level = 0) #NLOPT_LD_MMA, NLOPT_LD_SLSQP
  )
  
  # print(result)
  return(result)
}


# ---- Computing resulting functions : solve ode with the estimated parameters, given in result ----
sol_opti1 = function(pat_nb, Data, parameters, Y0 = NA){
  
  x1 = parameters[1]
  para <- list(sigma = parameters[2],
               rho = parameters[3],
               eta = parameters[4],
               mu = parameters[5],
               delta = parameters[6],
               alpha = parameters[7])
  
  T0 = 10^9
  
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  
  if (is.na(Y0)){
    Y0 = TC_pat[1]/T0
  }
  
  #Initial conditions
  init = c("X"=x1, "Y"=Y0) #y1 = y(tau=tau1) with tau1 = k2*K*T0*t/100, I take y1 = TC_pat/T0

  # Normalized time
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  time = seq(from=min(time_pat), to=max(time_pat), by=0.01)
  
  #Solve differential equations
  sol = ode(y=init, times = time, func = derive, parms = para, method="bdf")
  
  return(list(time=time, y=sol[,3], parameters=parameters))
}

# ---- Building data base for opti pb 1 ----
add_op1_Data = function(Data_preprocessed){
  Data_bis = Data_preprocessed
  
  nb_pat = length(Data_bis$Patient_Anonmyized)
  Data_bis$time = rep(NA,nb_pat)
  Data_bis$y_opt = rep(NA,nb_pat)
  Data_bis$parameters_opt = rep(NA,nb_pat)
  
  for (i in 1:nb_pat){ 
    print(paste("Patient ",i))
    result = opti1_pat(i, Data_bis)
    parameters = result$solution
    
    sol = sol_opti1(i, Data_bis, parameters)
    
    Data_bis$time[i] = I(list(sol$time))
    Data_bis$y_opt[i] = I(list(sol$y))
    Data_bis$parameters_opt[i] = I(list(parameters))
    save(Data_bis, file="./Data_processed/Data_bis.Rda")
  }
  save(Data_bis, file="./Data_processed/Data_bis.Rda")
  return(Data_bis)
}

# ---- Goodness of fit analysis ----
goodness_fit = function(pat_nb, Data){
    T0 = 10^9
    TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
    y_obs = TC_pat/T0
    time_pat = Data$Treatment_Day[[pat_nb]]
    time_pat = time_pat/(max(time_pat) - min(time_pat)) 
    time = Data$time[[pat_nb]]
    
    y = Data$y_opt[[pat_nb]]
    
    y_est = spline(time, y, xout = time_pat)$y
    
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


goodness_fit_analysis = function(Data){

  parameters_df = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_1 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_2 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_3 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_4 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  parameters_df_5 = data.frame(x1=numeric(), sigma=numeric(), rho=numeric(), eta=numeric(), mu=numeric(), delta=numeric(), alpha=numeric())
  
  GF_df = data.frame(MAE = numeric(), MSE = numeric(), R2 = numeric())
  
  
  for (i in 1:length(Data$Patient_Anonmyized)){
    
    GF_i = goodness_fit(i, Data)
    parameters_i = data.frame(x1=Data$parameters_opt[[i]][1], sigma=Data$parameters_opt[[i]][2], rho=Data$parameters_opt[[i]][3], eta=Data$parameters_opt[[i]][4], mu=Data$parameters_opt[[i]][5], delta=Data$parameters_opt[[i]][6], alpha=Data$parameters_opt[[i]][7])
    
    GF_df = rbind(GF_df, GF_i)
    parameters_df = rbind(parameters_df, parameters_i)
    
    
    j = substr(Data$Study_Arm[[i]], 7, 7)
    df = get(paste("parameters_df_", j, sep=""))
    
    df = rbind(df, parameters_i)
    
    assign(paste("parameters_df_", j, sep=""), df)
  }  
  
  Fit = list("GF" = GF_df, "Para_all" = parameters_df, "Para_1" = parameters_df_1, "Para_2" = parameters_df_2, "Para_3" = parameters_df_3, "Para_4" = parameters_df_4,"Para_5" = parameters_df_5)
  save(Fit, file = "./Data_processed/Fit.Rda")
  return(Fit)
}

# ---- Plotting results ----
plot_op1 = function(pat_nb, Data){
  par(mar = c(5, 4, 4, 5), mfrow=c(1,1))
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  ymax = max(max(Data$y_opt[[pat_nb]]), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(Data$y_opt[[pat_nb]]), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP1 results for patient :", pat_nb), ylim=c(ymin,ymax))
  
  
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  legend("topright", 
         legend = c("y_opt", "observations"), 
         col = c("black", "black"), 
         lty = c(1, NA),   # lty=NA pour les points
         pch = c(NA, 1), # pch=16 pour les points des mesures
         cex = 0.5, 
         xpd=TRUE, 
         inset = c(-0.25, 0)) 
}

boxplot_OP1 = function(parameters_opti1){
  par(mar = c(3, 5, 3, 6))
  boxplot(parameters_opti1$sigma, parameters_opti1$mu, parameters_opti1$delta, parameters_opti1$alpha, parameters_opti1$rho, parameters_opti1$eta, 
          names = c("σ","μ","δ","α","ρ","η"),
          ylim=c(10^(-2), 10^2), 
          log="y",
          col = NA, 
          ylab = "Parameters value (log scale)", 
          border = "blue",  
          medcol = "red",   
          whiskcol = "blue",  
          staplecol = "blue")
  
}
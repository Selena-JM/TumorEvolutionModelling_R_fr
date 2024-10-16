#date:03/10/2024, author : Séléna Jean-Mactoux
#Optimization pb 2 : Parameter identifiability 

source("derive.R")

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
#we use the sum the difference between pmin et pmax as a cost function

opti2_pat = function(pat_nb, Data){
  epsilon = 0.2
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  time = Data$time[[pat_nb]]
  y_opt = Data$y_opt[[pat_nb]]
  parameters_opt = Data$parameters_opt[[pat_nb]]
  
  f_minimize_OP2 = function(bounds) {
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
    print(paste("Patient", i))
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

# ---- Optimization problem 2 bis ----

opti2_bis_pat = function(pat_nb, Data, k){
  epsilon = 0.2
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  time = Data$time[[pat_nb]]
  y_opt = Data$y_opt[[pat_nb]]
  parameters_opt = Data$parameters_opt[[pat_nb]]
  
  f_minimize_OP2_bis = function(p_k) {
    p_k_min = p_k[1]
    p_k_max = p_k[2]
    
    #Cost function
    Cost = -p_k_max+p_k_min #- bc want to maximize
    return(Cost)
  }
  
  eval_jac_g_ineq_bis = function(p_k){
    jacobian(eval_g_ineq_bis, p_k)
  }
  
  eval_g_ineq_bis = function(p_k){
    para = rep(parameters_opt[2:7], 2)
    para[k] = p_k[1]
    para[6+k] = p_k[2]
    
    odes = compute_odes(pat_nb, Data, bounds=para)
    
    y_est_min = odes$y_est_min
    y_est_max = odes$y_est_max
    y_est_opt = spline(time, y_opt, xout = time_pat)$y
    
    ineq_values = c(
      y_est_max - (1 + epsilon) * y_est_opt,   # y_est_max <= (1 + epsilon) * y_est_opt
      (1 - epsilon) * y_est_opt - y_est_max,   # y_est_max >= (1 - epsilon) * y_est_opt
      y_est_min - (1 + epsilon) * y_est_opt,   # y_est_min <= (1 + epsilon) * y_est_opt
      (1 - epsilon) * y_est_opt - y_est_min    # y_est_min >= (1 - epsilon) * y_est_opt
    )
    
    return(ineq_values)
  }
  
  # Computing of gradient
  gradient_OP2_bis = function(p_k) {
    grad = grad(func = f_minimize_OP2_bis, x = p_k)
    # print(paste("Gradient:", paste(round(grad, 6), collapse=", ")))  # Display the gradient
    return(grad)
  }
  
  # Limits for x1 and parameters
  lower_bounds <- c(10^(-2), parameters_opt[k]) #lower bound of p_min is lower bound of all para : 10^-2 and lower bound of p_max is p_opt
  upper_bounds <- c(parameters_opt[k], 10^2) 
  
  #Starting values for the parameter optimization
  start_para = c("p_k_min" = upper_bounds[1], "p_k_max" = lower_bounds[2])
  
  # Optimization
  result <- nloptr(
    x0 = start_para,
    eval_f = f_minimize_OP2_bis,
    eval_grad_f = gradient_OP2_bis, # Using the defined gradient
    eval_g_ineq = eval_g_ineq_bis,
    eval_jac_g_ineq = eval_jac_g_ineq_bis,
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-8, "maxeval"=500, "print_level"=0) #NLOPT_LD_MMA, NLOPT_LD_SLSQP
  )
  
  # print(result)
  return(result)
}

# ---- Building data base for opti pb 2 bis ----
add_op2_bis_Data = function(Data_post_op1){
  Data_ter_bis = Data_post_op1
  
  nb_pat = length(Data_ter_bis$Patient_Anonmyized)
  
  # Data_ter_bis$y_min = rep(NA,nb_pat)
  # Data_ter_bis$y_max = rep(NA,nb_pat)
  # Data_ter_bis$parameters_min = rep(NA,nb_pat)
  # Data_ter_bis$parameters_max = rep(NA,nb_pat)
  
  parameters_min = rep(NA, 6)
  parameters_max = rep(NA, 6)
  for (i in 29:50){
    print(paste("Patient", i))
    # if (i!= 4 & i!=11 & i!=14 & i!=25){
    for (k in 1:6){
      print(paste("Patient", i, ", parameter", k))
      result = opti2_bis_pat(i, Data_ter_bis, k)
      
      p_k = result$solution
      parameters_min[k] = p_k[1]
      parameters_max[k] = p_k[2]
    
      # }
    }
    bounds = c(parameters_min, parameters_max)
    odes = compute_odes(i, Data_ter_bis, bounds)
    
    Data_ter_bis$y_min[i] = I(list(odes$y_min))
    Data_ter_bis$y_max[i] = I(list(odes$y_max))
    Data_ter_bis$parameters_min[i] = I(list(parameters_min))
    Data_ter_bis$parameters_max[i] = I(list(parameters_max))
    save(Data_ter_bis, file="./Data_processed/Data_ter_bis.Rda")
  }
  save(Data_ter_bis, file="./Data_processed/Data_ter_bis.Rda")
  return(Data_ter_bis)
}

# ---- Plotting results ----
plot_op1_op2 = function(pat_nb, Data_ter){
  par(mar = c(5, 4, 4, 5.5))
  
  ymax = max(max(Data_ter$y_min[[pat_nb]]), max(Data_ter$y_max[[pat_nb]]), max(Data_ter$y_opt[[pat_nb]]),max(Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(Data_ter$y_min[[pat_nb]]), min(Data_ter$y_max[[pat_nb]]), min(Data_ter$y_opt[[pat_nb]]),min(Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  plot(Data_ter$time[[pat_nb]], Data_ter$y_opt[[pat_nb]], type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP2 results for patient :", pat_nb), ylim=c(ymin,ymax))
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

plot_op1_op2_bis = function(pat_nb, Data_ter, para){
  par(mar = c(5, 4, 4, 5.5))
  
  parameters_opt = Data_ter$parameters_opt[[pat_nb]]
  parameters = rep(parameters_opt[2:7], 2)
  parameters[para] = Data_ter$parameters_min[[pat_nb]][para]
  parameters[6+para] = Data_ter$parameters_max[[pat_nb]][para]
  
  odes = compute_odes(pat_nb, Data_ter, bounds=parameters)
  
  y_min = odes$y_min
  y_max = odes$y_max
  
  ymax = max(max(y_min), max(y_max), max(Data_ter$y_opt[[pat_nb]]),max(Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(y_min), min(y_max), min(Data_ter$y_opt[[pat_nb]]),min(Data_ter$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  plot(Data_ter$time[[pat_nb]], Data_ter$y_opt[[pat_nb]], type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP2 results for patient :", pat_nb), ylim=c(ymin,ymax))
  lines(Data_ter$time[[pat_nb]], y_min, type = 'l', col = 'red')
  lines(Data_ter$time[[pat_nb]], y_max, type = 'l', col = 'blue')
  
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

plot_range_para = function(pat_nb, Data){
  p_opt = Data$parameters_opt[[pat_nb]]
  sigma = p_opt[1]
  rho = p_opt[2]
  eta = p_opt[3]
  mu = p_opt[4]
  delta = p_opt[5]
  alpha = p_opt[6]
  
  p_min = Data$parameters_min[[pat_nb]]
  low_value=c(p_min[1], p_min[4], p_min[5], p_min[6], p_min[2], p_min[3])
  
  p_max = Data$parameters_max[[pat_nb]]
  high_value=c(p_max[1], p_max[4], p_max[5], p_max[6], p_max[2], p_max[3])
  
  print(Data$parameters_opt[[pat_nb]])
  print(Data$parameters_min[[pat_nb]])
  print(Data$parameters_max[[pat_nb]])
  
  data = data.frame(
    Parameter = factor(c("σ", "μ", "δ", "α", "ρ", "η"), 
                       levels = c("σ", "μ", "δ", "α", "ρ", "η")), 
    y = c(sigma, mu, delta, alpha, rho, eta),
    
    low_Value = low_value,
    high_Value = high_value,  
    Color = brewer.pal(6, "Set2") #c("blue", "green", "lightblue", "orange", "yellow", "red")
    )
  
  # Créer le graphique
  ggplot(data, aes(x = Parameter, y = y)) +
    geom_point(aes(color = Color), size = 2, show.legend = FALSE) +  
    geom_segment(aes(x = Parameter, xend = Parameter, 
                     y = y, yend = low_Value, color=Color), 
                 arrow = arrow(length = unit(0.2, "cm"), ends = "last", type="open"), 
                 size = 0.8,
                 show.legend = FALSE) +  
    geom_segment(aes(x = Parameter, xend = Parameter, 
                     y = y, yend = high_Value, color=Color), 
                 arrow = arrow(length = unit(0.2, "cm"), ends = "last", type="open"),
                 size = 0.8,
                 show.legend = FALSE) + 
      
    scale_color_manual(values = data$Color) +  
    scale_y_log10(limits = c(0.01, 100)) +  
    labs(title = paste("Patient ", pat_nb), x = "", y = "") +  
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), 
          axis.text.x = element_text(size = 14),  # Taille des lettres grecques
          plot.title = element_text(hjust = 0.5, size = 16))  # Centrer le titre
}


plot_histo = function(Data){
  # nb_pat = length(Data$Patient_Anonmyized)
  
  nb_pat = 50
  
  range = 10^seq(from=-2, to=2, by=0.1)
  L = length(range)
  D_histo = matrix(0, nrow = 6, ncol = L)
  
  for (i in 1:nb_pat){
    p_min = Data$parameters_min[[i]]
    p_max = Data$parameters_max[[i]]
    intervals = list(c(p_min[1], p_max[1]), c(p_min[2], p_max[2]), c(p_min[3], p_max[3]), c(p_min[4], p_max[4]), c(p_min[5], p_max[5]), c(p_min[6], p_max[6]))
    
    belongings = sapply(range, function(p) {
      which(sapply(intervals, function(interval) {
        p >= interval[1] && p <= interval[2]
      }))
    })
    
    #count the number of points in results belonging to interval 1, 2 etc
    for (p in 1:L){
      for (k in 1:6){
        if (length(belongings[[p]][belongings[[p]]==k]) != 0){
          D_histo[k,p] = D_histo[k,p] + 1
        }
      }
    }
  }

  par(mfrow=c(3,2), mar=c(2,2,2,2))
  title = c("σ", "ρ", "η", "μ", "δ", "α")
  Color = brewer.pal(6, "Set2")
  for (k in 1:6){
    plot(range, D_histo[k,]/nb_pat, log='x', type='l', main=title[k], col=Color[k])
  }
}

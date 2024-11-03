#date:21/10/24, author : Séléna Jean-Mactoux

# ---- Finding parameter intervals for sensitivity analysis ----
#finding the parameter intervals with other parameters set to optimal values
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
  lower_bounds <- c(10^(-2), parameters_opt[2:7][k]) #lower bound of p_min is lower bound of all para : 10^-2 and lower bound of p_max is p_opt
  upper_bounds <- c(parameters_opt[2:7][k], 10^2) 
  
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
  
  return(result)
}

# ---- Building data base intervals for sensitivity analysis ----
add_op2_bis_Data = function(Data_post_op1){
  Data_ter_bis = Data_post_op1
  
  nb_pat = length(Data_ter_bis$Patient_Anonmyized)
  
  Data_ter_bis$y_min = rep(NA,nb_pat)
  Data_ter_bis$y_max = rep(NA,nb_pat)
  Data_ter_bis$parameters_min = rep(NA,nb_pat)
  Data_ter_bis$parameters_max = rep(NA,nb_pat)
  
  for (i in 1:nb_pat){
    print(paste("Patient", i))
    parameters_min = rep(NA, 6)
    parameters_max = rep(NA, 6)
    for (k in 1:6){
      print(paste("Patient", i, ", parameter", k))
      result = opti2_bis_pat(i, Data_ter_bis, k)
      
      p_k = result$solution
      parameters_min[k] = p_k[1]
      parameters_max[k] = p_k[2]
      
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


# ---- Analysis ----
sensitivity_analysis_multidimensional = function(pat_nb, Data_ter_bis, nb_obs, N = 100, bounds="OP2"){
  #don't do the analysis on x1 because take the bounds found by OP2
  if (bounds=="OP2"){
    param_bounds = data.frame(
      p1 = c(Data_ter_bis$parameters_min[[pat_nb]][1], Data_ter_bis$parameters_max[[pat_nb]][1]),   
      p2 = c(Data_ter_bis$parameters_min[[pat_nb]][2], Data_ter_bis$parameters_max[[pat_nb]][2]),   
      p3 = c(Data_ter_bis$parameters_min[[pat_nb]][3], Data_ter_bis$parameters_max[[pat_nb]][3]),
      p4 = c(Data_ter_bis$parameters_min[[pat_nb]][4], Data_ter_bis$parameters_max[[pat_nb]][4]),
      p5 = c(Data_ter_bis$parameters_min[[pat_nb]][5], Data_ter_bis$parameters_max[[pat_nb]][5]),
      p6 = c(Data_ter_bis$parameters_min[[pat_nb]][6], Data_ter_bis$parameters_max[[pat_nb]][6])
    )
  }
  else {
    param_bounds = data.frame(
      p1 = c(10^(-2), 10^2),   
      p2 = c(10^(-2), 10^2),   
      p3 = c(10^(-2), 10^2),
      p4 = c(10^(-2), 10^2),
      p5 = c(10^(-2), 10^2),
      p6 = c(10^(-2), 10^2)
    )
  }
  param_bounds = t(param_bounds) #want the parameters on the rows

  model = function(parameters){
    # parameters contains samples in [0,1], so need to adjust the samples to fit the interval of each parameter
    
    scaled_params = t(apply(parameters, 1, function(x) {
      10^(x * (log10(param_bounds[, 2]) - log10(param_bounds[, 1])) + log10(param_bounds[, 1]))
    }))

    TC_pat = Data_ter_bis$TargetLesionLongDiam_mm[[pat_nb]]
    time_pat = Data_ter_bis$Treatment_Day[[pat_nb]]
    time_pat = time_pat / (max(Data_ter_bis$Treatment_Day[[pat_nb]]) - min(Data_ter_bis$Treatment_Day[[pat_nb]]))

    T0 = 10^9

    # output = matrix(NA, nrow=N, ncol=1)
    output = numeric(dim(scaled_params)[1])
    for (i in 1:dim(scaled_params)[1]){
      para = c(Data_ter_bis$parameters_opt[[pat_nb]][1], scaled_params[i, ])
      sol = sol_opti1(pat_nb=pat_nb, Data=Data_ter_bis, parameters=para)
      y = sol$y
      time = sol$time

      y_est = spline(time, y, xout = time_pat)$y

      #Cost function
      Cost = sum(abs(y_est-TC_pat/T0)^2)
      # Cost = y_est[nb_obs]
      
      output[i] = Cost
    }

    return(output)
  }

  # Generate samples
  # X1 = randtoolbox::sobol(n = N, dim = 6)  # 6 parameters so 6 dimensions
  # X2 = randtoolbox::sobol(n = N, dim = 6)
  
  X1 = data.frame(matrix(runif(6*N, 0, 1), ncol = 6))
  X2 = data.frame(matrix(runif(6*N, 0, 1), ncol = 6))
  
  sobol_design = sobolSalt(model = model, X1 = as.data.frame(X1), X2 = as.data.frame(X2), nboot = 100, scheme = "B")
  return(sobol_design)
}


# ---- Plot ----
sensitivity_plot_unidimensional = function(pat_nb, Data_ter_bis, N = 100, bounds="OP2_bis"){
  #don't do the analysis on x1 because take the bounds found by OP2
  if (bounds=="OP2_bis"){
    param_bounds = data.frame(
      p1 = c(Data_ter_bis$parameters_min[[pat_nb]][1], Data_ter_bis$parameters_max[[pat_nb]][1]),   
      p2 = c(Data_ter_bis$parameters_min[[pat_nb]][2], Data_ter_bis$parameters_max[[pat_nb]][2]),   
      p3 = c(Data_ter_bis$parameters_min[[pat_nb]][3], Data_ter_bis$parameters_max[[pat_nb]][3]),
      p4 = c(Data_ter_bis$parameters_min[[pat_nb]][4], Data_ter_bis$parameters_max[[pat_nb]][4]),
      p5 = c(Data_ter_bis$parameters_min[[pat_nb]][5], Data_ter_bis$parameters_max[[pat_nb]][5]),
      p6 = c(Data_ter_bis$parameters_min[[pat_nb]][6], Data_ter_bis$parameters_max[[pat_nb]][6])
    )
  }
  else if (bounds=="OP2_ter"){
    param_bounds = data.frame(
      p1 = c(Data_ter_bis$parameters_min[[pat_nb]][1,1], Data_ter_bis$parameters_max[[pat_nb]][1,1]),   
      p2 = c(Data_ter_bis$parameters_min[[pat_nb]][2,2], Data_ter_bis$parameters_max[[pat_nb]][2,2]),   
      p3 = c(Data_ter_bis$parameters_min[[pat_nb]][3,3], Data_ter_bis$parameters_max[[pat_nb]][3,3]),
      p4 = c(Data_ter_bis$parameters_min[[pat_nb]][4,4], Data_ter_bis$parameters_max[[pat_nb]][4,4]),
      p5 = c(Data_ter_bis$parameters_min[[pat_nb]][5,5], Data_ter_bis$parameters_max[[pat_nb]][5,5]),
      p6 = c(Data_ter_bis$parameters_min[[pat_nb]][6,6], Data_ter_bis$parameters_max[[pat_nb]][6,6])
    )
  }
  else {
    param_bounds = data.frame(
      p1 = c(10^(-2), 10^2),   
      p2 = c(10^(-2), 10^2),   
      p3 = c(10^(-2), 10^2),
      p4 = c(10^(-2), 10^2),
      p5 = c(10^(-2), 10^2),
      p6 = c(10^(-2), 10^2)
    )
  }
  param_bounds = t(param_bounds) #want the parameters on the rows
  
  N_tot=6*N
  
  param_samples = matrix(NA, nrow = N_tot, ncol = 6)
  
  # Generate samples
  for (i in 1:6) {
    param_samples[((i-1)*N+1):(i*N), ] = rep(1,6)
    # param_samples[((i-1)*N+1):(i*N), i] = randtoolbox::sobol(n = N, dim = 1)
    param_samples[((i-1)*N+1):(i*N), i] = seq(from=0, to=1, length.out = N)
  }
  
  model = function(parameters){
    # parameters contains samples in [0,1], so need to adjust the samples to fit the interval of each parameter
    
    scaled_params = matrix(NA, nrow=N_tot, ncol=6)
    for (i in 1:6) {
      scaled_params[((i-1)*N+1):(i*N), ] = parameters[((i-1)*N+1):(i*N), ] * matrix(rep(Data_ter_bis$parameters_opt[[pat_nb]][2:7], each=N), N, 6)
      scaled_params[((i-1)*N+1):(i*N), i] = 10^(parameters[((i-1)*N+1):(i*N), i] * (log10(param_bounds[i, 2]) - log10(param_bounds[i, 1])) + log10(param_bounds[i, 1]))
    }

    output = matrix(NA, nrow = N_tot, ncol = length(sol_opti1(pat_nb = pat_nb, Data = Data_ter_bis, parameters = c(Data_ter_bis$parameters_opt[[pat_nb]][1], scaled_params[1, ]))$y))
    for (i in 1:N_tot){
      para = c(Data_ter_bis$parameters_opt[[pat_nb]][1], scaled_params[i, ])
      sol = sol_opti1(pat_nb=pat_nb, Data=Data_ter_bis, parameters=para)$y
      output[i,] = sol
    }
    
    return(output)
  }
  
  
  output = model(param_samples)
  
  #plot
  par(mfrow=c(3,2), mar=c(2,2,2,2), oma=c(0, 0, 2, 0))
  letters = c("σ", "ρ", "η", "μ", "δ", "α")
  
  for (i in 1:6){
    beg = ((i-1)*N+1)
    end = (i*N)
    ymax = max(output[beg:end, ])
    ymin = min(output[beg:end, ])
    
    plot(Data_ter_bis$time[[pat_nb]], output[beg,], type='l', ylim=c(ymin,ymax), xlab="Normalised time",
         ylab="Nb of tumor cells / 10^9", main=paste(letters[i]," in [", round(param_bounds[i,1],3), ",", round(param_bounds[i,2],3), "]", sep=""))
    
    for (j in (beg+1):end){
      lines(Data_ter_bis$time[[pat_nb]],output[j,], type='l')
    }
    
    # ttle = paste("Unidimensional sensitivity analysis, patient", pat_nb)
    ttle = paste("Patient", pat_nb)
    mtext(ttle, outer=TRUE, cex=1.5)
  }
}

sensitivity_plot_histo = function(sobol_design_multi, pat_nb){
  results_inter = sobol_design_multi$S2
  
  results_inter = data.frame(interactions = rownames(results_inter), results_inter, row.names = NULL)
  
  colnames(results_inter)[1] = "interactions"
  
  # create the histogram
  ggplot(results_inter, aes(x = reorder(interactions, original), y = original)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_errorbar(aes(ymin = min..c.i., ymax = max..c.i.), width = 0.2) +
    labs(title = paste("Indices de Sobol de Deuxième Ordre, patient ", pat_nb, sep=""),
         x = "Interactions entre Paramètres",
         y = "Valeur de l'Indice") +
    coord_flip() +  
    theme_minimal()
  
}

sensitivity_plot_heatmap = function(sobol_design_multi, pat_nb){
  parameter_names = title = c("σ", "ρ", "η", "μ", "δ", "α")
  interaction_matrix = matrix(0, nrow = length(parameter_names), ncol = length(parameter_names),
                               dimnames = list(parameter_names, parameter_names))
  results_inter = sobol_design_multi$S2
  
  # Fill the matrix
  index = 1
  for (i in 1:5) {
    for (j in (i+1):6){
      interaction_matrix[i,j] = results_inter$original[index]
      interaction_matrix[j,i] = results_inter$original[index]  # Pour rendre la matrice symétrique
      index = index + 1}
  }
  
  heatmap_data = as.data.frame(as.table(interaction_matrix))
  colnames(heatmap_data) = c("Param1", "Param2", "Value")
  
  # Create heatmap
  ggplot(heatmap_data, aes(Param1, Param2)) +
    geom_tile(aes(fill = Value), color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", 
                         midpoint = 0, limit = c(-0.15, 0.15), 
                         name = "Interaction Index") +
    theme_minimal() +
    labs(title = paste("Heatmap des indices d'interaction de Sobol, patient ", pat_nb, sep=""),
         x = "Paramètres",
         y = "Paramètres")
  
}

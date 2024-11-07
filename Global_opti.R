# date:07/10/2024, author:Séléna Jean-Mactoux
source("derive.R")

# ---- Global parameters authentification ----
interval_analysis = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)

  borne_up_min = sapply(1:6, function(i) {
    max(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  })
  borne_low_min = sapply(1:6, function(i) {
    min(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  })
  
  borne_low_max = sapply(1:6, function(i) {
    min(sapply(Data$parameters_max, function(x) x[i]), na.rm=TRUE)
  })
  
  
  borne_up_max = sapply(1:6, function(i) {
    max(sapply(Data$parameters_max, function(x) x[i]), na.rm=TRUE)
  })
  
  
  #The parameters are distributed across the whole spectrum of possible values
  #Verify that max is always >= min
  
  l_bool = rep(NA, nb_pat)
  for (i in 1:nb_pat){
    l_bool[i] = length(which(Data$parameters_max[[i]] < Data$parameters_min[[i]]))
  }
  print(paste("Number of patients with max < min :", length(which(l_bool>0))))
  
  
  return(list(borne_low_min = borne_low_min, borne_up_min = borne_up_min, borne_low_max = borne_low_max, borne_up_max = borne_up_max))
}


# ---- Clustering ----
#We could try to identify clusters of patients for which the parameters are approximately the same

#can use the mean between the min and max for each parameter for clustering
clustering_extremums = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)
  Data_clustering = matrix(NA, nrow = nb_pat, ncol=20)
  for (pat in 1:nb_pat){
    Data_clustering[pat,1] = pat
    Data_clustering[pat,2:7] = Data$parameters_min[[pat]]
    Data_clustering[pat,8:14] = Data$parameters_opt[[pat]]
    Data_clustering[pat,15:20] = Data$parameters_max[[pat]]
  }
  
  #scale parameters for the clustering
  Data_clustering[,2:20] = scale(Data_clustering[,2:20], center = TRUE, scale = TRUE)
  
  return(Data_clustering)
}


#clustering with the 8 parameters (x1, parameters, Y0)
clustering_para = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)
  Data_clustering = matrix(NA, nrow = nb_pat, ncol=9)
  for (pat in 1:nb_pat){
    Data_clustering[pat,1] = pat
    Data_clustering[pat,2:8] = Data$parameters_opt[[pat]]
    Data_clustering[pat,9] = Data$TargetLesionLongDiam_mm[[pat]][1]/10^9
  }

  #scale parameters for the clustering
  Data_clustering[,2:9] = scale(Data_clustering[,2:9], center = TRUE, scale = TRUE)
  return(Data_clustering)
}


#clustering with the curves
clustering_curves = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)
  time_cluster = seq(from=0.05, to=0.95, by=0.01) 
  
  Data_clustering = matrix(NA, nrow = nb_pat, ncol=length(time_cluster)+1)
  for (pat in 1:nb_pat){
    if (length(which(is.na(Data$time[[pat]]))) == 0){
      Data_clustering[pat,1] = pat
      Data_clustering[pat,2:(length(time_cluster)+1)] = spline(Data$time[[pat]], Data$y_opt[[pat]], xout = time_cluster)$y
    }
  }
  
  #scale parameters for the clustering
  Data_clustering[,2:(length(time_cluster)+1)] = scale(Data_clustering[,2:(length(time_cluster)+1)], center = TRUE, scale = TRUE)
  
  return(Data_clustering)
}

# ---- Cluster analysis ----
# analyzing the means, median and std of each cluster, to get the parameters 
# for the cluster curve
cluster_analysis = function(cluster_hierarchical, nb_clusters, Data_clustering, Data_ter_bis){
  cluster_means = matrix(NA, nrow=nb_clusters, ncol=7)
  cluster_std = matrix(NA, nrow=nb_clusters, ncol=7)
  cluster_median = matrix(NA, nrow=nb_clusters, ncol=7)
  cluster_Y0 = rep(NA, nb_clusters)
  
  for (i in 1:nb_clusters){
    indices = Data_clustering[unique(which(cluster_hierarchical==i)),1]
    parameters = matrix(NA, nrow = length(indices), ncol = 7) #contains the parameters for each patients
    lY0 = rep(NA, length(indices))
    for (j in 1:length(indices)){
      parameters[j,] = Data_ter_bis$parameters_opt[[indices[j]]]
      lY0[j] = Data_ter_bis$TargetLesionLongDiam_mm[[indices[j]]][1]/10^9
    }
    
    cluster_means[i,] = apply(parameters, 2, function(col) mean(col, na.rm=TRUE))
    cluster_median[i,] = apply(parameters, 2, function(col) median(col, na.rm=TRUE))
    cluster_std[i,] = apply(parameters, 2, function(col) sd(col, na.rm=TRUE))
    cluster_Y0[i] = mean(lY0)
  }
  return(list(cluster_means=cluster_means, cluster_std=cluster_std, cluster_median=cluster_median, cluster_Y0=cluster_Y0))
}

#Optimizing the parameters of the cluster curve
global_opti = function(cluster_hierarchical, cluster_nb, Data_clustering, Data, maxeval=500, precision = 10^(-8)){
  T0 = 10^9
  indices = Data_clustering[unique(which(cluster_hierarchical==cluster_nb)),1]
  print(indices)
  
  f_minimize_global = function(parameters) {
    # Parameters
    x1 = parameters[1]
    para <- list(
      sigma = parameters[2],
      rho = parameters[3],
      eta = parameters[4],
      mu = parameters[5],
      delta = parameters[6],
      alpha = parameters[7])
    # Y0 = parameters[8]
    
    y_err = rep(NA, length(indices))
    
    for (pat in 1:length(indices)){
      #Initial conditions
      Y0 = Data$TargetLesionLongDiam_mm[[indices[pat]]][1]/T0
      init = c("X"=x1, "Y"=Y0)
      
      time_pat = Data$Treatment_Day[[indices[pat]]]
      time_pat = time_pat / (max(Data$Treatment_Day[[indices[pat]]]) - min(Data$Treatment_Day[[indices[pat]]]))
      
      time = seq(from=min(time_pat), to=max(time_pat), by=0.01) #*(max(time_pat)-min(time_pat))
      
      #Solve differential equations
      ODE_sol = ode(y=init, times = time, func = derive, parms = para, method="bdf")
      y_cluster = ODE_sol[,3]
      
      y_cluster_est = spline(time, y_cluster, xout = time_pat)$y
      
      # y_pat = spline(Data$time[[indices[pat]]], Data$y_opt[[indices[pat]]], xout = time_pat)$y
      y_pat = Data$TargetLesionLongDiam_mm[[indices[pat]]]/T0
      y_err[pat] = sum((y_pat - y_cluster_est)^2)
    }
    
    #Cost function
    Cost = sum(y_err)
    return(Cost)
  }
  
  
  # Computing of gradient
  gradient_global <- function(parameters) {
    grad <- grad(func = f_minimize_global, x = parameters)
    return(grad)
  }
  
  # Limits for x1 and parameters
  lower_bounds <- c(10^(-2), rep(c(10^(-2)), 6))
  upper_bounds <- c(10^2, rep(c(10^2), 6)) 
  
  #Starting values for the parameter optimization
  start_para = c("x1" = 1, "sigma" = 1, "rho" = 1, "eta" = 1, "mu" = 1, "delta" = 1, "alpha" = 1)# , "Y0"=1
  
  # Optimization
  result <- nloptr(
    x0 = start_para,
    eval_f = f_minimize_global,
    eval_grad_f = gradient_global, 
    lb = lower_bounds,
    ub = upper_bounds,
    opts = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = precision, "maxeval"=maxeval, print_level = 1) 
  )
  
  # print(result)
  return(result)
}
# ---- Plotting results ----
plot_cluster_curve = function(pat_nb, parameters, Data, Y0=NA){
  
  sol = sol_opti1(pat_nb, Data, parameters, Y0)
  time = sol$time
  y = sol$y
  
  par(mar = c(5, 4, 4, 5))
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  ymax = max(max(y), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(y), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
  plot(time, y, type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Global opti results for patient :", pat_nb), ylim=c(ymin,ymax))
  
  
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  legend("topright", 
         legend = c("y_global_opt", "observations"), 
         col = c("black", "black"), 
         lty = c(1, NA),   # lty=NA pour les points
         pch = c(NA, 1), # pch=16 pour les points des mesures
         cex = 0.6, 
         xpd=TRUE, 
         inset = c(-0.25, 0)) 
}

plot_clusters = function(cluster_nb, cluster_hierarchical, Data_clustering, Data, parameters, ylim=TRUE){
  pat_cluster = Data_clustering[unique(which(cluster_hierarchical==cluster_nb)),1]
  nb_pat_cluster = length(pat_cluster)
  
  pat_nb = pat_cluster[1]
  
  sol = sol_opti1(pat_nb, Data, parameters)
  time = sol$time
  y = sol$y - Data$TargetLesionLongDiam_mm[[pat_nb]][1]/10^9
  
  par(mar = c(5, 4, 4, 5.5))
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  ymax=max(y)
  ymin=min(y)
  if(ylim==TRUE){
    for (i in 1:2){
      if (max(Data$y_opt[[pat_cluster[i]]]- Data$TargetLesionLongDiam_mm[[pat_cluster[i]]][1]/10^9) > ymax){
        ymax=max(Data$y_opt[[pat_cluster[i]]]- Data$TargetLesionLongDiam_mm[[pat_cluster[i]]][1]/10^9)
      }
      if (min(Data$y_opt[[pat_cluster[i]]] - Data$TargetLesionLongDiam_mm[[pat_cluster[i]]][1]/10^9) < ymin){
        ymin=min(Data$y_opt[[pat_cluster[i]]]- Data$TargetLesionLongDiam_mm[[pat_cluster[i]]][1]/10^9)
      }
    }
  }
  
  plot(time, y, type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Global opti results for cluster :", cluster), ylim=c(ymin,ymax),lwd = 1.5)
  
  colors = rainbow(length(pat_cluster))
  
  lgd = "Cluster"
  lwd = 1.5
  col = "black"
  # for (i in 1:nb_pat_cluster){
  for (i in 1:2){
    pat_nb = pat_cluster[i]
    lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]] - Data$TargetLesionLongDiam_mm[[pat_nb]][1]/10^9, col=colors[i])
    lgd = c(lgd, paste("Opti", pat_nb))
    col = c(col, colors[i])
    lwd = c(lwd,1)
  }

  legend("topright", 
         legend = lgd, 
         col = col,
         xpd = TRUE, 
         lwd = lwd,
         cex = 0.7, 
         inset = c(-0.25, 0)) 
}

# date:07/10/2024, author:Séléna Jean-Mactoux

# ---- Global parameters authentification ----
global_opti = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)
  
  # borne_min = sapply(1:nb_pat, function(i) {
  #   max(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  # })
  borne_up_min = sapply(1:6, function(i) {
    max(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  })
  borne_low_min = sapply(1:6, function(i) {
    min(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  })
  
  # borne_max = sapply(1:nb_pat, function(i) {
  #   min(sapply(Data$parameters_max, function(x) x[i]), na.rm=TRUE)
  # })
  
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
clustering = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)
  Data_clustering = matrix(NA, nrow = nb_pat, ncol=6)
  for (pat in 1:nb_pat){
    for (para in 1:6){
      Data_clustering[pat,para] = (Data$parameters_min[[pat]][para]+Data$parameters_max[[pat]][para])/2
    }
  }
  return(Data_clustering)
}


plot_global = function(pat_nb, kmeans_result, Data){
  cluster_pat = kmeans_result$cluster[pat_nb]
  x1 = Data$parameters_opt[[pat_nb]][1]
  parameters = c(x1, kmeans_result$centers[cluster_pat, ])
  
  sol = sol_opti1(pat_nb, Data, parameters)
  time = sol$time
  y = sol$y
  
  par(mar = c(5, 4, 4, 5.5))
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  ymax = max(max(y), max(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  ymin = min(min(time), min(Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9))
  
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

plot_clusters = function(cluster, kmeans_result, Data, ylim=TRUE){
  pat_cluster = which(kmeans_result$cluster == cluster)
  nb_pat_cluster = length(pat_cluster)
  
  pat_nb = pat_cluster[1]
  x1 = Data$parameters_opt[[pat_nb]][1]
  parameters = c(x1, kmeans_result$centers[cluster, ])
  
  sol = sol_opti1(pat_nb, Data, parameters)
  time = sol$time
  y = sol$y
  
  par(mar = c(5, 4, 4, 5.5))
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  ymax=max(y)
  ymin=min(y)
  if(ylim==TRUE){
    for (i in 1:length(pat_cluster)){
      if (max(Data$y_opt[[pat_cluster[i]]]) > ymax){
        ymax=max(Data$y_opt[[pat_cluster[i]]])
      }
      if (min(Data$y_opt[[pat_cluster[i]]]) < ymin){
        ymin=min(Data$y_opt[[pat_cluster[i]]])
      }
    }
  }
  
  plot(time, y, type = 'l', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Global opti results for cluster :", cluster), ylim=c(ymin,ymax),lwd = 1.5)
  
  colors = rainbow(length(pat_cluster))
  
  lgd = "Cluster"
  lwd = 1.5
  col = "black"
  for (i in 1:nb_pat_cluster){
    pat_nb = pat_cluster[i]
    lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], col=colors[i])
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

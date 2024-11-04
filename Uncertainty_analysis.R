#date 30/10/24, author Séléna Jean-Mactoux
#Uncertainty analysis on the patients selected for the report

# ---- Sample creation ----
sample_creation = function(pat_nb, Data, N, min=-2, max=2){
  nb_points = length(Data$TargetLesionLongDiam_mm[[pat_nb]])
  sample_pat = matrix(NA, ncol = nb_points, nrow = N)

  # converting LD in mm to number of tumor cells
  vol_tc = 8*10^(-6) #mm^3/TC, volume of a tumor cell
  prop_tc = 3/4 # proportion of TC in the lesion volume
  
  for (i in 1:N){
    sample_pat[i,] = as.numeric(Data$TargetLesionLongDiam_mm[[pat_nb]]) + runif(nb_points,min=min, max=max)

    LD = sample_pat[i,] 
    vol_lesion = (4/3)*pi*(LD/2)^3 #mm^3, spherical shape of the lesion : LD is the long diameter used
    nb_tc = round(vol_lesion*prop_tc/vol_tc) #tc occupy 3/4 of the lesion volume
    sample_pat[i,] = nb_tc 
  }
  return(sample_pat)
}

# ---- Uncertainty propagation ----
compute_uncertainty_OP1 = function(pat_nb, Data, sample, nb_points_omitted=0){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat / (max(Data$Treatment_Day[[pat_nb]]) - min(Data$Treatment_Day[[pat_nb]]))
  time = seq(from=min(time_pat), to=max(time_pat), by=0.01) 
  
  T0 = 10^9
  results_uncertainty_propagation = matrix(NA, nrow=dim(sample)[1], ncol=length(time))
  
  for (i in 1:dim(sample)[1]){
    print(paste("Sample", i))
    TC_sample = sample[i,]
    
    result = opti1_pat(pat_nb, Data, TC_pat = TC_sample, nb_points_omitted=nb_points_omitted)
    parameters = result$solution

    sol = sol_opti1(pat_nb=pat_nb, Data=Data, parameters=parameters, TC_pat = TC_sample)
    y = sol$y

    results_uncertainty_propagation[i,] = y
  }
  return(results_uncertainty_propagation)
}

compute_uncertainty_OP3 = function(pat_nb, Data, sample, results_uncertainty_Pred){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat / (max(Data$Treatment_Day[[pat_nb]]) - min(Data$Treatment_Day[[pat_nb]]))
  time = seq(from=min(time_pat), to=max(time_pat), by=0.01) 
  
  T0 = 10^9
  results_uncertainty_OP3_ylow = matrix(NA, nrow=dim(results_uncertainty_Pred)[1], ncol=length(time))
  results_uncertainty_OP3_yupp = matrix(NA, nrow=dim(results_uncertainty_Pred)[1], ncol=length(time))
  
  
  for (i in 1:dim(sample)[1]){
    print(paste("Sample", i))
    TC_sample = sample[i,]
    
    result = opti3_pat(pat_nb, Data, TC_pat = TC_sample, y_pred=results_uncertainty_Pred[i,])
    
    bounds = result$solution
    parameters_low = bounds[1:7]
    parameters_upp = bounds[8:14]
    
    odes = compute_odes_op3(pat_nb, Data, bounds, plot=FALSE)
    
    results_uncertainty_OP3_ylow[i,] = odes$y_low
    results_uncertainty_OP3_yupp[i,] = odes$y_upp
  }
  return(list(results_uncertainty_OP3_ylow=results_uncertainty_OP3_ylow, results_uncertainty_OP3_yupp=results_uncertainty_OP3_yupp))
}

# ---- Goodness of fit analysis ----
goodness_fit_uncertainty = function(pat_nb, Data, results_OP1_uncertainty){
  T0 = 10^9
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  y_obs = TC_pat/T0
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat)) 
  
  time = Data$time[[pat_nb]]
  
  GF = data.frame(MAE = numeric(), MSE = numeric(), R2 = numeric())
  
  for (i in 1:nrow(results_OP1_uncertainty)){
    y_opt = results_OP1_uncertainty[i,]
    
    y_est = spline(time, y_opt, xout = time_pat)$y
    
    # Metrics
    # R²
    ss_total = sum((y_obs - mean(y_est))^2)
    ss_residual = sum((y_obs - y_est)^2)
    r_squared = 1 - (ss_residual / ss_total)
    
    # MAE
    mae = mean(abs(y_obs - y_est))
    
    # MSE
    mse = mean((y_obs - y_est)^2)
    
    GF = rbind(GF, data.frame(MAE = mae, MSE = mse, R2 = r_squared))
    
  }
  return(GF)
}


goodness_fit_analysis_uncertainty = function(Data, patients){
  GF_tot = data.frame(MAE = numeric(), MSE = numeric(), R2 = numeric())
  
  for (i in patients){
    load(file=paste("./Data_processed/Uncertainty_OP1_pat", i, ".Rda", sep=""))
    
    GF_pat = goodness_fit_uncertainty(i, Data, Uncertainty_OP1)
    save(GF_pat, file=paste("./Data_processed/GF_analysis_uncertainty_pat", i, ".Rda", sep=""))
    
    GF_tot = rbind(GF_tot, GF_pat)
  }  
  
  save(GF_tot, file = "./Data_processed/GF_tot_uncertainty.Rda")
  return(GF_tot)
}



# ---- Plots ----
plot_uncertainty_OP1=function(pat_nb, Data, results_uncertainty_propagation){

  ymin = min(results_uncertainty_propagation)
  ymax = max(results_uncertainty_propagation)
  
  par(mfrow=c(1,1), mar=c(5,4,4,5.5))
  
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type='n', ylim=c(ymin,ymax), xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Uncertainty analysis patient ", pat_nb, sep=""), col='red', lwd=1.5)
  
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat / (max(Data$Treatment_Day[[pat_nb]]) - min(Data$Treatment_Day[[pat_nb]]))
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  for (i in 1:nrow(results_uncertainty_propagation)){
    lines(Data$time[[pat_nb]], results_uncertainty_propagation[i,], type='l', col='black')
  }
  
  lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type='l', col="red", lwd=1.5)

  legend("topright", 
         legend = c("y_opt", "observations", "Sample curves"), 
         col = c("red", "black", "black"), 
         lty = c(1, NA, 1),   # lty=NA pour les points
         pch = c(NA, 1, NA), # pch=16 pour les points des mesures
         lwd = c(1.5, NA, 1),
         cex = 0.5, 
         xpd=TRUE, 
         inset = c(-0.25, 0)) 
}


plot_uncertainty_pred = function(pat_nb, Data, results_uncertainty_propagation, nb_points_omitted=2){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-nb_points_omitted+1]-0.01
  
  ymin = min(results_uncertainty_propagation)
  ymax = max(results_uncertainty_propagation)
  
  # plot environment
  par(mar = c(5, 4, 4, 5.5), mfrow=c(1,1))
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'n', col = 'red', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("Prediction results for patient :", pat_nb), ylim=c(ymin, ymax), lwd=1.5)
  
  points(time_pat, Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)
  
  # Rectangle
  rect(xleft = x_col, ybottom = par("usr")[3], xright = par("usr")[2], 
       ytop = par("usr")[4], col = rgb(0.9, 0.8, 0.8, 0.5), border = NA) 
  
  for (i in 1:nrow(results_uncertainty_propagation)){
    lines(Data$time[[pat_nb]], results_uncertainty_propagation[i,], type='l', col='black')
  }
  
  lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'l', col = 'red', lwd = 1.5)
  lines(Data$time[[pat_nb]], Data$y_pred[[pat_nb]], type = 'l', col = 'lightblue', lwd = 1.5)
  
  legend("topright", 
         legend = c("y_opt", "y_pred", "observations", "Sample curves"), 
         col = c("red", "lightblue", "black", "black"), 
         lty = c(1, 1, NA, 1),   # lty=NA pour les points
         pch = c(NA, NA, 1, NA), # pch=16 pour les points des mesures
         lwd = c(1.5, 1.5, NA, 1),
         cex = 0.5, 
         xpd=TRUE, 
         inset = c(-0.25, 0)) 
}

plot_uncertainty_pred_op3 = function(pat_nb, Data, results_uncertainty_OP3_ylow, results_uncertainty_OP3_yupp){
  time_pat = Data$Treatment_Day[[pat_nb]]
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  
  sorted_time_pat = sort(Data$Treatment_Day[[pat_nb]])
  sorted_time_pat = sorted_time_pat/(max(sorted_time_pat) - min(sorted_time_pat))
  x_col = sorted_time_pat[length(sorted_time_pat)-Data$nb_points_omitted[pat_nb]+1]-0.01
  
  ymin = min(max(min(results_uncertainty_OP3_ylow), min((Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)*0.6)), min(Data$y_pred[[pat_nb]]), min(Data$y_opt[[pat_nb]]))
  ymax = max(min(max(results_uncertainty_OP3_yupp), max((Data$TargetLesionLongDiam_mm[[pat_nb]]/10^9)*1.4)), max(Data$y_pred[[pat_nb]]), max(Data$y_opt[[pat_nb]]))
  
  # ymin = min(min(results_uncertainty_OP3_ylow), min(Data$y_pred[[pat_nb]]), min(Data$y_opt[[pat_nb]]))
  # ymax = max(max(results_uncertainty_OP3_yupp), max(Data$y_pred[[pat_nb]]), max(Data$y_opt[[pat_nb]]))
  
  # plot environment
  par(mar = c(5, 4, 4, 5), mfrow=c(1,1))
  plot(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'n', col = 'black', xlab="Normalised time",
       ylab="Nb of tumor cells / 10^9", main=paste("OP3 results for patient :", pat_nb), ylim=c(ymin, ymax))
  
  # Rectangle
  rect(xleft = x_col, ybottom = par("usr")[3], xright = par("usr")[2], 
       ytop = par("usr")[4], col = rgb(0.9, 0.8, 0.8, 0.5), border = NA)
  
  # Plot 
  index_max = which.max(results_uncertainty_OP3_yupp[,ncol(results_uncertainty_OP3_yupp)])
  y_upp=results_uncertainty_OP3_yupp[index_max,]
  
  
  # mask = results_uncertainty_OP3_yupp[-index_max,]
  # index_max_bis = which.max(mask[,ncol(results_uncertainty_OP3_yupp)])
  # mask_bis = mask[-index_max_bis,]
  # index_max_ter = which.max(mask_bis[,ncol(results_uncertainty_OP3_yupp)])
  # 
  # y_upp=mask[index_max_ter,]
  
  
  index_min =  which.min(results_uncertainty_OP3_ylow[,ncol(results_uncertainty_OP3_ylow)])
  y_low = results_uncertainty_OP3_ylow[index_min,]
  
  # min_global = min(results_uncertainty_OP3_ylow[, ncol(results_uncertainty_OP3_ylow)])
  # mask = results_uncertainty_OP3_ylow[, ncol(results_uncertainty_OP3_ylow)] != min_global
  # index_min_bis = which(mask)[which.min(results_uncertainty_OP3_ylow[mask, ncol(results_uncertainty_OP3_ylow)])]
  # 
  # y_low=results_uncertainty_OP3_ylow[index_min_bis,]
  
  
  
  lines(Data$time[[pat_nb]], Data$y_opt[[pat_nb]], type = 'l', col = 'black')
  lines(Data$time[[pat_nb]], Data$y_pred[[pat_nb]], type = 'l', col = 'black', lty = 4)
  for(i in 1:nrow(results_uncertainty_OP3_ylow)){
    y_low = results_uncertainty_OP3_ylow[i,]
    y_upp = results_uncertainty_OP3_yupp[i,]
    lines(Data$time[[pat_nb]], y_low, col = 'red')
    lines(Data$time[[pat_nb]], y_upp, col = 'blue')
  }
  
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
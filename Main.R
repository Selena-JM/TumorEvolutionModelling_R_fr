# date : 28/09/2024, author : Séléna Jean-Mactoux
# Reproduction of the paper : 
# "Can the Kuznetsov Model Replicate and Predict Cancer Growth in Humans?" 
# by Mohammad El Wajeh, Falco Jung, Dominik Bongartz, Chrysoula Dimitra Kappatou, 
# Narmin Ghaffari Laleh, Alexander Mitsos, Jakob Nikolas Kather


rm(list=ls())

create_data_frames = FALSE

# ---- Libraries ----
library("readxl")
library("deSolve")
library("nloptr")
library("numDeriv")
library("ggplot2")
library("plotly")
library("stats")
library("scales")
library("RColorBrewer")
library("sensitivity")
library("randtoolbox")

# ---- Dependencies ---- 
source("Preprocessing.R")
source("Opti_1.R")
source("Opti_2.R")
source("Opti_3.R")
source("Predictions.R")
source("derive.R")
source("Global_opti.R")
source("Uncertainty_analysis.R")
source("Sensitivity_analysis.R")

# ---- Data preprocessing and processing : building data bases ----

if (create_data_frames == TRUE){
  #### Preprocessing ####
  preprocessing()
  load(file = "./Data_processed/Data.Rda")
  load(file = "./Data_processed/Data_converted.Rda")
  
  #### Optimization problem 1 : finding optimal parameters for each patient ####
  print("OP1")
  Data_OP1 = add_op1_Data(Data_converted)
  
  #goodness of fit analysis for optimal curves (results of OP1)
  Fit = goodness_fit_analysis(Data_OP1)
  
  #### Optimization problem 2 : Parameter identifiability ####
  print("OP2 Method 1")
  Data_OP2_1 = add_op2_Data_1(Data_OP1)
  
  print("OP2 Method 2")
  Data_OP2_2 = add_op2_Data_2(Data_OP1)
  
  #### Predictions ####
  print("Predictions")
  Data_Pred = add_pred_Data(Data_OP2_1, 2)
  
  #### Optimization problem 3 : Parameter identifiability ####
  print("OP3")
  Data_OP3 = add_op3_Data(Data_Pred)
  
  #Evaluation of predictions and uncertainty intervals
  Fit_OP3 = goodness_prediction_intervals_analysis(Data_OP3)
  
  #### Local unidimensional sensitivity analysis ####
  print("OP2 sensitivity analysis")
  Data_OP2_sensi = add_op2_Data_sensi(Data_OP1)

} else{
  files = list.files(path = "./Data_processed/")
  
  for (file in files) {
    load(file=paste("./Data_processed/", file, sep=""))
  }
}


# ---- Plots ----
#### DB OP1 : Data_OP1 ####
pat_nb = 2
plot_op1(pat_nb, Data_OP1)

# visualize parameters over all patients
parameters = Fit$Para_all
boxplot_OP1(parameters)

# Goodness of fit analysis 
GF = Fit$GF
print(summary(GF$MAE, rm.na=True))
print(summary(sqrt(GF$MSE), rm.na=True))
print(summary(GF$R2, rm.na=True))

#### DB OP2 : Data_OP2_2 ####
pat_nb = 2
plot_range_para(pat_nb, Data_OP2_1, method="1")
plot_range_para(pat_nb, Data_OP2_2, method="2")

# Plot curves
plot_op1_op2_method2(pat_nb, Data_OP2_2)

#plot parameters histogram
parameters = Fit$Para_all
plot_histo(Data_OP2_1, parameters, method="1")
plot_histo(Data_OP2_2, parameters, method="2")

#### Predictions : Data_Pred ####
pat_nb = 2
plot_pred(pat_nb, Data_Pred)

#### DB OP3 : Data_OP3 ####
pat_nb = 2
plot_pred_op3(pat_nb, Data_OP3)

# Analysis
print(summary(Fit_OP3$MAE, rm.na=True))
print(summary(sqrt(Fit_OP3$MSE), rm.na=True))
print(summary(Fit_OP3$R2, rm.na=True))
print(summary(Fit_OP3$nb_points_in_interval, rm.na=True))

# ---- Uncertainty analysis ----
#### Proportion of tumor size ####
for (pat_nb in c(2, 11, 92, 169, 251, 252)){
  print(2/as.numeric(Data$TargetLesionLongDiam_mm[[pat_nb]][1]))
}

#### OP1 ####
if (create_data_frames==TRUE){
  for (pat_nb in c(2, 11, 92, 169, 251, 252)){
    sample = sample_creation(pat_nb, Data, 10)
    
    Uncertainty_OP1 = compute_uncertainty_OP1(pat_nb, Data, sample)
    save(Uncertainty_OP1, file=paste("./Data_processed/Uncertainty_OP1_pat", pat_nb, ".Rda", sep=""))
  }
}

#plot results
pat_nb = 2
load(file=paste("./Data_processed/Uncertainty_OP1_pat", pat_nb, ".Rda", sep=""))
plot_uncertainty_OP1(pat_nb, Data_OP1, Uncertainty_OP1)

#### Goodness of fit analysis ####
GF_tot = goodness_fit_analysis_uncertainty(Data_OP1, c(2,11,92,169,251,252))
print(summary(GF_tot$MAE, rm.na=True))
print(summary(sqrt(GF_tot$MSE), rm.na=True))
print(summary(GF_tot$R2, rm.na=True))

load(file=paste("./Data_processed/GF_analysis_uncertainty_pat", 251, ".Rda", sep=""))
print(summary(GF_pat$MAE, rm.na=True))
print(summary(sqrt(GF_pat$MSE), rm.na=True))
print(summary(GF_pat$R2, rm.na=True))
Opti_251 = goodness_fit(251, Data_OP1)
print(Opti_251)

load(file=paste("./Data_processed/GF_analysis_uncertainty_pat", 252, ".Rda", sep=""))
print(summary(GF_pat$MAE, rm.na=True))
print(summary(sqrt(GF_pat$MSE), rm.na=True))
print(summary(GF_pat$R2, rm.na=True))
Opti_252 = goodness_fit(252, Data_OP1)
print(Opti_252)

#### Predictions ####
if (create_data_frames==TRUE){
  for (pat_nb in c(2,11,92,169,251,252)){
    sample = sample_creation(pat_nb, Data, 10)
    save(sample, file=paste("./Data_processed/sample_Pred_pat", pat_nb, ".Rda", sep=""))
    
    Uncertainty_Pred = compute_uncertainty_OP1(pat_nb, Data, sample, nb_points_omitted = 2)
    save(Uncertainty_Pred, file=paste("./Data_processed/Uncertainty_Pred_pat", pat_nb, ".Rda", sep=""))
  }
}

#plot results
pat_nb=2
load(file=paste("./Data_processed/Uncertainty_Pred_pat", pat_nb, ".Rda", sep=""))
plot_uncertainty_pred(pat_nb, Data_Pred, Uncertainty_Pred)

#### Worse predictions ####
if (create_data_frames==TRUE){
  for (pat_nb in c(2,11,92,169,251,252)){
    load(file=paste("./Data_processed/sample_Pred_pat", pat_nb, ".Rda", sep=""))
    load(file=paste("./Data_processed/Uncertainty_Pred_pat", pat_nb, ".Rda", sep=""))
    
    if (pat_nb == 2){ #patient 2 has convergence issues with last sample
      Uncertainty_OP3 = compute_uncertainty_OP3(pat_nb, Data_Pred, sample[-nrow(sample),], Uncertainty_Pred[-nrow(Uncertainty_Pred),])
      save(Uncertainty_OP3, file=paste("./Data_processed/Uncertainty_OP3_pat", pat_nb, ".Rda", sep=""))
    }
    else{
      Uncertainty_OP3 = compute_uncertainty_OP3(pat_nb, Data_Pred, sample, Uncertainty_Pred)
      save(Uncertainty_OP3, file=paste("./Data_processed/Uncertainty_OP3_pat", pat_nb, ".Rda", sep=""))
    }
  }
}

#plot results
pat_nb=2
load(file=paste("./Data_processed/Uncertainty_OP3_pat", pat_nb, ".Rda", sep=""))
plot_uncertainty_pred_op3(pat_nb, Data_Pred, Uncertainty_OP3$results_uncertainty_OP3_ylow, Uncertainty_OP3$results_uncertainty_OP3_yupp)


# ---- Sensitivity analysis ----
#### Unidimensional ####
pat_nb = 2
N=1000
sensitivity_plot_unidimensional(pat_nb, Data_OP2_sensi, N=N)

#### Multidimensional ####
pat_nb=2
N = 1000

sobol_design_multi = sensitivity_analysis_multidimensional(pat_nb, Data_OP2_sensi, N=N, bounds="OP2")

print(sobol_design_multi)

#plot results
par(mfrow=c(1,1))
plot(sobol_design_multi, choice=1)

sensitivity_plot_histo(sobol_design_multi, pat_nb)
sensitivity_plot_heatmap(sobol_design_multi,2)


# ---- Global optimization ----
database = Data_OP2_sensi

#visualizing the max of the min and the min of the max parameters
bornes = interval_analysis(database)
print(bornes) #no intersection between the parameter intervals of all the patients

#trying hierarchical clustering
Data_clustering = clustering_extremums(database)
distance_matrix = dist(Data_clustering[,-1], method = "euclidean")

hc_result = hclust(distance_matrix, method = "ward.D")
par(mfrow=c(1,1))
plot(hc_result, main="Patient Dendrogramme", cex = 0.5, xlab="")

#cutting the tree
nb_clusters = 10
cluster_hierarchical = cutree(hc_result, k = nb_clusters)

#### Means of parameters ####
analysis = cluster_analysis(cluster_hierarchical, nb_clusters, Data_clustering, database)
cluster_means = analysis$cluster_means
cluster_std = analysis$cluster_std
cluster_median = analysis$cluster_median
cluster_Y0 = analysis$cluster_Y0

print(round(cluster_means,2))
print(round(cluster_std,2))

pat_nb = 1
cluster = cluster_hierarchical[pat_nb]

#indices of the patients in the same cluster
indices = Data_clustering[unique(which(cluster_hierarchical==cluster)),1]
print(indices)

#plot cluster curve with data points of patient pat_nb
plot_cluster_curve(pat_nb, cluster_means[cluster,],database)

#### Optimization of parameters ####
pat_nb = 1
cluster = cluster_hierarchical[pat_nb]

#indices of the patients in the same cluster
indices_clustering = unique(which(cluster_hierarchical==cluster))
print(indices_clustering)

#function to optimize parameters of the cluster curve
result_global_opti = global_opti(cluster_hierarchical, cluster, Data_clustering, database)

#results
parameters = result_global_opti$solution
plot_cluster_curve(pat_nb, parameters[1:7], database)


#### Parameter comparison between methods #### 
print(round(parameters,2))
print(round(cluster_means[cluster,], 2))
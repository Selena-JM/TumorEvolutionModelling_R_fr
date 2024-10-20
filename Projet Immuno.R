# date : 28/09/2024, author : Séléna Jean-Mactoux
# Reproduction of the paper : 
# "Can the Kuznetsov Model Replicate and Predict Cancer Growth in Humans?" 
# by Mohammad El Wajeh, Falco Jung, Dominik Bongartz, Chrysoula Dimitra Kappatou, 
# Narmin Ghaffari Laleh, Alexander Mitsos, Jakob Nikolas Kather

rm(list=ls())


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

# ---- Importing Data ----
# Data imported from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009822#pcbi.1009822.s006
options(digits = 19)
Study1 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study1") #peut importer les ID en text pour mieux gérer le fait que c'est des grands nombres
Study2 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study2")
Study3 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study3")
Study4 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study4")
Study5 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study5")

Data_studies = rbind(Study1, Study2, Study3, Study4, Study5)


# ---- Data pre-processing ----
# removing null inputs
Data_raw = data.frame(Patient_Anonmyized = numeric(), Treatment_Day = numeric(), TargetLesionLongDiam_mm = character(), Study_Arm=character())
for (i in 1:length(Data_studies$Patient_Anonmyized)){
  if (Data_studies$TargetLesionLongDiam_mm[i] != "NOT EVALUABLE" & Data_studies$TargetLesionLongDiam_mm[i] != "TOO SMALL TO MEASURE" & !is.na(Data_studies$TargetLesionLongDiam_mm[i])){
    Data_raw = rbind(Data_raw, Data_studies[i,])
  }
}

# patients with at least 6 measurements different from each other 
nb_patients_tot = 0
Data = data.frame(Patient_Anonmyized = numeric(), Treatment_Day = list(), TargetLesionLongDiam_mm = list(), Study_Arm=character())
i = 1
L_tot = length(Data_raw$Patient_Anonmyized)
bool_verif = TRUE
pat_elim = 0

while (i <= L_tot){
  id_pat = Data_raw$Patient_Anonmyized[i]
  nb_patients_tot = nb_patients_tot + 1
  
  nb_rec = 1 #nb of measurements for patient starting at line i
  #list_0 = c()
  while ((i + nb_rec) < L_tot & Data_raw$Patient_Anonmyized[i + nb_rec] == id_pat){
    # if (Data_raw$TargetLesionLongDiam_mm[i + nb_rec] == 0) {
    #   list_0 = c(list_0, i+nb_rec)
    # }
    nb_rec = nb_rec + 1
  }
  
  index_end = i+nb_rec-1
  l_target = unique(Data_raw$TargetLesionLongDiam_mm[i:index_end])
  nb_rec_unique = length(l_target)
  
  if (nb_rec_unique >= 6){ #(nb_rec - length(list_0)) >= 6
    # if (length(list_0) > 0) {
    #   index_end = min(i+nb_rec-1, list_0[1]) #prend en compte le premier 0 mais pas les autres
    # }
    # else {
    #   index_end = i+nb_rec-1
    # }

    indices = which(!duplicated(Data_raw$TargetLesionLongDiam_mm[i:index_end]))
    l_treatment = Data_raw$Treatment_Day[i:index_end]
    l_treatment = l_treatment[indices]

    Patient = data.frame(Patient_Anonmyized = id_pat, Treatment_Day = I(list(l_treatment)), TargetLesionLongDiam_mm = I(list(l_target)), Study_Arm=Data_raw$Study_Arm[i])
    Data = rbind(Data, Patient)

    if (length(unique(Data_raw$Patient_Anonmyized[i:(i+nb_rec-1)])) > 1) {
      bool_verif = FALSE
    }
  }
  
  else if (nb_rec >= 6){
    pat_elim = pat_elim+1
  }
  
  i = i + nb_rec
}

# Verifications

print(paste("Number patients before : ", nb_patients_tot, "/1472", sep = ""))
print(paste("Number patients after : ", length(Data$Patient_Anonmyized),  "/210", sep = ""))
print(paste("No mixing patients : ", bool_verif))
print(paste("Number of patients with too many duplicates : ", pat_elim))


# converting LD in mm to number of tumor cells
vol_tc = 8*10^(-6) #mm^3/TC, volume of a tumor cell
prop_tc = 3/4 # proportion of TC in the lesion volume

for (i in 1:length(Data$Patient_Anonmyized)){
  LD = as.numeric(Data$TargetLesionLongDiam_mm[[i]]) #this column was characters because of the "NOT EVALUABLE3 etc
  vol_lesion = (4/3)*pi*(LD/2)^3 #mm^3, spherical shape of the lesion : LD is the long diameter used
  nb_tc = round(vol_lesion*prop_tc/vol_tc) #tc occupy 3/4 of the lesion volume
  Data$TargetLesionLongDiam_mm[[i]] = nb_tc 
}

save(Data, file="./Data_processed/Data.Rda")

# To know the number of patients in each study arm
nb_pat_study=rep(NA,14)
names = c("Study_1_Arm_1", "Study_1_Arm_2", "Study_1_Arm_3", 
          "Study_2_Arm_1", "Study_2_Arm_2", 
          "Study_3_Arm_1", "Study_3_Arm_2", "Study_3_Arm_3", "Study_3_Arm_4", "Study_3_Arm_5", "Study_3_Arm_6", 
          "Study_4_Arm_1", "Study_4_Arm_2", 
          "Study_5_Arm_1")
for (i in 1:length(nb_pat_study)){
  nb_pat_study[i] = length(Data$Study_Arm[Data$Study_Arm == names[i]])
}

# ---- Data processing : building data bases ----
#### Optimization problem 1 : finding optimal parameters for each patient ####
source("Opti_1.R")
Data_bis = add_op1_Data(Data)
Fit = goodness_fit_analysis(Data_bis)

#### Optimization problem 2 : Parameter identifiability ####
source("Opti_2.R")
Data_ter = add_op2_Data(Data_bis)
Data_ter_bis = add_op2_bis_Data(Data_bis)

#### Predictions ####
source("Opti_3.R")
Data_qua = add_pred_Data(Data_ter, 2)

#### Optimization problem 3 : Parameter identifiability ####
Data_5 = add_op3_Data(Data_qua)
# Data_5 = add_op3_Data(Data_5)


# ---- Tests ----
#### DB OP1 : Data_bis ####
pat_nb = 92
plot_op1(pat_nb, Data_bis)

# visualize parameters over all patients
parameters = Fit$Para_5
boxplot_OP1(parameters)

# Goodness of fit analysis 
GF = Fit$GF
print(summary(GF$MAE, rm.na=True))
print(summary(sqrt(GF$MSE), rm.na=True))
print(summary(GF$R2, rm.na=True))

#### DB OP2 : Data_ter ####
pat_nb = 252
plot_range_para(pat_nb, Data_ter_bis)

# Plot curves
plot_op1_op2(pat_nb, Data_ter)
plot_op1_op2_bis(pat_nb, Data_ter_bis, para=1)

#plot paramters histogram
parameters = Fit$Para_all
plot_histo(Data_ter_bis, parameters)

#### Predictions : Data_qua ####
pat_nb = 2
plot_pred(pat_nb, Data_qua)
# plot_pred_test(pat_nb, Data_ter, nb_points_omitted=2, maxeval=500)


#### DB OP3 : Data_5 ####
pat_nb = 252
# plot_pred_op3_test(pat_nb, Data_qua)
plot_pred_op3(pat_nb, Data_5)


# ---- Global optimization ----
source("Global_opti.R")
database = Data_ter

#visualizing the max of the min and the min of the max parameters
bornes = global_opti(database)
print(bornes)

#creating data for clustering -> takes the optimal parameters as factors in the
# clustering because best way to take into  consideration the form of the curves
# for each patient
Data_clustering = clustering_extremums(database)

#### hierarchical clustering ####
Data_clustering = clustering_para(database)

distance_matrix = dist(Data_clustering[,-1], method = "euclidean")

hc_result = hclust(distance_matrix, method = "ward.D2")
plot(hc_result, main="Patient Dendrogramme")

nb_clusters = 10
cluster_hierarchical = cutree(hc_result, k = nb_clusters)

#finding the centers of each cluster
analysis = cluster_analysis(cluster_hierarchical, nb_clusters, Data_clustering, database)
cluster_means = analysis$cluster_means
cluster_std = analysis$cluster_std
cluster_median = analysis$cluster_median
cluster_Y0 = analysis$cluster_Y0


pat_nb_clustering = 1
cluster = cluster_hierarchical[[pat_nb_clustering]]
pat_nb = Data_clustering[pat_nb_clustering,1]

indices = Data_clustering[unique(which(cluster_hierarchical==cluster)),1]
print(indices)

plot_cluster_curve(pat_nb, cluster_means[cluster,],database)

result_global_opti = global_opti(cluster_hierarchical, cluster, Data_clustering, database)

parameters = result_global_opti$solution
plot_cluster_curve(pat_nb, parameters[1:7], parameters[8], database)


#### clustering kmeans ####
Data_clustering = clustering_extremums(database)

nb_clusters=20

kmeans_result = kmeans(Data_clustering[,-1], centers = nb_clusters, nstart = 25)
print(kmeans_result)

analysis = cluster_analysis(kmeans_result$cluster, nb_clusters, Data_clustering, database)
cluster_means = analysis$cluster_means
cluster_std = analysis$cluster_std
cluster_median = analysis$cluster_median
cluster_Y0 = analysis$cluster_Y0

pat_nb_clustering = 45
cluster = kmeans_result$cluster[[pat_nb_clustering]]
pat_nb = Data_clustering[pat_nb_clustering,1]

indices = Data_clustering[unique(which(kmeans_result$cluster==cluster)),1]
print(indices)

plot_cluster_curve(pat_nb, cluster_means[cluster,], database)
# plot_cluster_curve(pat_nb, kmeans_result$centers[1:7], database)

# plot_clusters(5, kmeans_result, database)


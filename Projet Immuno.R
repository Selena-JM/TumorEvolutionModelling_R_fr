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
library("plotly")

# ---- Importing Data ----
# Data imported from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009822#pcbi.1009822.s006
options(digits = 19)
Study1 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study1") #peut importer les ID en text pour mieux gérer le fait que c'est des grands nombres
Study2 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study2")
Study3 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study3")
Study4 = read_excel("Data.xlsx", col_types = c("guess", "guess", "guess", "guess"), sheet="Study4")

Data_studies = rbind(Study1, Study2, Study3, Study4)


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


# ---- Data processing : building data bases ----
#### Model definition : differential equations ####
derive = function(t, var, para){
  with(as.list(c(var, para)),{
    
    X = var[1]
    Y = var[2]
    
    deriX = sigma + (rho*X*Y)/(eta+Y) - delta*X - mu*X*Y
    deriY = alpha*Y - (10^(-2))*X*Y # E0/T0 = (10^7)/(10^9)
    return(list(c(deriX, deriY)))
  })
}

#### Optimization problem 1 : finding optimal parameters for each patient ####
source("Opti_1.R")
Data_bis = add_op1_Data(Data)
Fit = goodness_fit_analysis(Data_bis)

#### Goodness of fit analysis ####
GF = Fit[[1]]
print(summary(GF$MAE, rm.na=True))
print(summary(sqrt(GF$MSE), rm.na=True))
print(summary(GF$R2, rm.na=True))

#### Optimization problem 2 : Parameter identifiability ####
source("Opti_2.R")
Data_ter = add_op2_Data(Data_bis)

#### Predictions ####
source("Opti_3.R")
Data_qua = add_pred_Data(Data_ter, 2)

#### Optimization problem 3 : Parameter identifiability ####
Data_5 = add_op3_Data(Data_qua)
# Data_5 = add_op3_Data(Data_5)


# ---- Tests ----
#### DB OP1 : Data_bis ####
pat_nb = 11
plot_op1(pat_nb, Data_bis)
boxplot_OP1(pat_nb, parameters_df)

#### DB OP2 : Data_ter ####
pat_nb = 9
plot_op1_op2(pat_nb, Data_ter)

#### Predictions : Data_qua ####
pat_nb = 9
plot_pred(pat_nb, Data_qua)
# plot_pred_test(pat_nb, Data_ter, nb_points_omitted=2, maxeval=500)


#### DB OP3 : Data_5 ####
pat_nb = 9
plot_pred_op3_test(pat_nb, Data_qua)
# plot_pred_op3(pat_nb, Data_5)


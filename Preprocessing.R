# date:07/10/2024, author:Séléna Jean-Mactoux

preprocessing = function(){
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
    
    nb_rec = 1 
    
    while ((i + nb_rec) < L_tot & Data_raw$Patient_Anonmyized[i + nb_rec] == id_pat){
      nb_rec = nb_rec + 1
    }
    
    index_end = i+nb_rec-1
    l_target = unique(Data_raw$TargetLesionLongDiam_mm[i:index_end])
    nb_rec_unique = length(l_target)
    
    if (nb_rec_unique >= 6){ 
      
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
  
  
  save(Data, file="./Data_processed/Data.Rda")
  
  # converting LD in mm to number of tumor cells
  vol_tc = 8*10^(-6) #mm^3/TC, volume of a tumor cell
  prop_tc = 3/4 # proportion of TC in the lesion volume
  
  
  Data_converted = Data
  for (i in 1:length(Data_converted$Patient_Anonmyized)){
    LD = as.numeric(Data_converted$TargetLesionLongDiam_mm[[i]]) #this column was characters because of the "NOT EVALUABLE3 etc
    vol_lesion = (4/3)*pi*(LD/2)^3 #mm^3, spherical shape of the lesion : LD is the long diameter used
    nb_tc = round(vol_lesion*prop_tc/vol_tc) #tc occupy 3/4 of the lesion volume
    Data_converted$TargetLesionLongDiam_mm[[i]] = nb_tc 
  }
  
  
  save(Data_converted, file="./Data_processed/Data_converted.Rda")
  
  # To know the number of patients in each study arm
  nb_pat_study=rep(NA,14)
  names = c("Study_1_Arm_1", "Study_1_Arm_2", "Study_1_Arm_3", 
            "Study_2_Arm_1", "Study_2_Arm_2", 
            "Study_3_Arm_1", "Study_3_Arm_2", "Study_3_Arm_3", "Study_3_Arm_4", "Study_3_Arm_5", "Study_3_Arm_6", 
            "Study_4_Arm_1", "Study_4_Arm_2", 
            "Study_5_Arm_1")
  for (i in 1:length(nb_pat_study)){
    nb_pat_study[i] = length(Data_converted$Study_Arm[Data_converted$Study_Arm == names[i]])
  }
}
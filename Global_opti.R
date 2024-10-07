# date:07/10/2024, author:Séléna Jean-Mactoux

# ---- Global parameters authentification ----
global_opti = function(Data){
  nb_pat = length(Data$Patient_Anonmyized)
  
  # borne_min = sapply(1:nb_pat, function(i) {
  #   max(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  # })
  borne_min = sapply(1:6, function(i) {
    max(sapply(Data$parameters_min, function(x) x[i]), na.rm=TRUE)
  })
  
  # borne_max = sapply(1:nb_pat, function(i) {
  #   min(sapply(Data$parameters_max, function(x) x[i]), na.rm=TRUE)
  # })
  
  borne_max = sapply(1:6, function(i) {
    min(sapply(Data$parameters_max, function(x) x[i]), na.rm=TRUE)
  })
  
  return(list(borne_min = borne_min, borne_max = borne_max))
}
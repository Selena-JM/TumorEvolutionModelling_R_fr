# date : 30/09/2024, author : Séléna Jean-Mactoux
#Brouillon, choses qui peuvent être utiles mais pas dans le code principal

#A faire : 
# - Vérifications de la partie pre processing
# x Trouver meilleur moyen de faire comparaison, parce que la je prend le point le plus proche mais il faudrait une interpolation plutot
# x Améliorer les plots
# x Réorganiser le code en fonctions intelligible
# x Caractériser performance de OP1
# x Trouver pourquoi moi je dois utiliser moindre carrés pour que ça marche et pas eux 
# - Trouver le meilleur algo possible d'opti, voir ce qui se rapproche de ce qu'ils ont utilisé
# - Disclosure ?
# - Faire une base de donnée pour les paramètres optimaux pour chaque patient et résultats op2
# x Passer à OP2
# x Inspecter pourquoi erreur de ode sur certains patients
# - Inspecter pourquoi GROSSES erreur de modèle sur certains patients et les enlever ?
# - Inspecter l'ordination des temps (pour patient 5 c'est pas dans l'ordre) et est ce que ça influe les mauvais résultats
# - Essayer de comprendre pourquoi les prédictions sont si mauvaises (y compris sur les valeurs prises en compte dans l'optimisation) alors que quand on enlève 0 points c'est très bien



# ---- Fonctions pour visualiser le coût en fonction de 2 paramètres (les autres sont fixés) ----
pat_nb = 2
f_minimize_OP1 = function(parameters) {
  #y_pat = data for the patient
  #temps_pat = times of data for the patient, normalised
  #para = parameters for the patient, [sigma, rho, eta, mu, delta, alpha]
  
  # Parameters
  x1 = parameters[1]
  para <- list(#x1 = parameters[1],
    sigma = parameters[2],
    rho = parameters[3],
    eta = parameters[4],
    mu = parameters[5],
    delta = parameters[6],
    alpha = parameters[7])
  T0 = 10^9
  
  TC_pat = Data$TargetLesionLongDiam_mm[[pat_nb]]
  time_pat = Data$Treatment_Day[[pat_nb]]
  
  #Initial conditions
  init = c("X"=x1, "Y"=TC_pat[1]/T0) #y1 = y(tau=tau1) with tau1 = k2*K*T0*t/100, I take y1 = TC_pat/T0
  
  # Normalised time
  time_pat = time_pat/(max(time_pat) - min(time_pat))
  time = seq(from=min(time_pat), to=max(time_pat), by=0.01*(max(time_pat)-min(time_pat)))
  
  
  #Solve differential equations
  y = ode(y=init, times = time, func = derive, parms = para, method="lsoda")
  
  
  #Finding the closest times of time_pat in time
  y_est = c()
  index = 0
  for (i in 1:length(time_pat)){
    index = which.min(abs(time - time_pat[i]))
    y_est = c(y_est, y[index])
  }
  
  #Cost function
  Cost = sum(abs(y_est-TC_pat/T0)^2) #no abs() in what they wrote but weird
  return(Cost)
}

##
# Fonction pour créer une grille de valeurs des paramètres
evaluate_cost_surface <- function(param1_range, param2_range, fixed_parameters, pat_nb) {
  grid <- expand.grid(param1_range, param2_range)
  cost_values <- matrix(0, nrow = length(param1_range), ncol = length(param2_range))
  
  for (i in seq_along(param1_range)) {
    for (j in seq_along(param2_range)) {
      # Fixer sigma et rho avec les valeurs de la grille
      parameters <- c(fixed_parameters[1], param1_range[i], param2_range[j], 
                      fixed_parameters[4], fixed_parameters[5], fixed_parameters[6], fixed_parameters[7])
      
      # Évaluer la fonction de coût pour ces paramètres
      cost_values[i, j] <- f_minimize_OP1(parameters)
    }
  }
  return(list(cost = cost_values, param1 = param1_range, param2 = param2_range))
}

# Choisir les intervalles des deux paramètres que vous souhaitez explorer
param1_range <- seq(0.01, 10, length.out = 50)  # Par exemple pour sigma
param2_range <- seq(0.01, 10, length.out = 50)  # Par exemple pour rho

# Les autres paramètres restent fixes
fixed_parameters <- c(1, 1, 1, 1, 1, 1, 1)  # x1, sigma, rho, eta, mu, delta, alpha

# Évaluer la surface de coût
cost_data <- evaluate_cost_surface(param1_range, param2_range, fixed_parameters, pat_nb = 1)

# Tracé de la surface avec plotly
plot_ly(x = cost_data$param1, y = cost_data$param2, z = cost_data$cost, type = "surface") %>%
  layout(scene = list(xaxis = list(title = 'sigma'),
                      yaxis = list(title = 'rho'),
                      zaxis = list(title = 'Coût')))
##


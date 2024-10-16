#date:16/10/2024, author: Séléna Jean-Mactoux
#Definition of model ODEs
derive = function(t, var, para){
  with(as.list(c(var, para)),{
    
    X = var[1]
    Y = var[2]
    
    deriX = sigma + (rho*X*Y)/(eta+Y) - delta*X - mu*X*Y
    deriY = alpha*Y - (10^(-2))*X*Y # E0/T0 = (10^7)/(10^9)
    return(list(c(deriX, deriY)))
  })
}
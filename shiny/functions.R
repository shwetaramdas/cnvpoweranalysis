##functions.R

#power calculator

##Function 1: Poisson with 100% purity
poisson_pure = function(alpha,W,l,D,N,L){
#  print("Here!")
  n = sqrt(L/W)
  numerator = W*D*abs(N-2)*sqrt(L/W)/(2*l)
  denominator = sqrt(N*W*D/(2*l))
  
  t = numerator / denominator
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  print(power)
  return(power)
}



##Function 2: Poisson with < 100% purity: INCOMPLETE
poisson_impurity = function(alpha,W,l,D,N,L,F){
  n = sqrt(L/W)
  numerator = F*N*sqrt(D*L)*(l)^(-1/2)
  denominator = 2
  
  t = numerator / denominator
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(power)
}


##Function 3: Negative binomial with 100% purity
negativebinomial = function(alpha,W,l,D,N,L,phi){
  n = sqrt(L/W)
  numerator = sqrt(D*L)*abs(N-2)
  denominator = sqrt(4*N*l*(1 + W*D*N*phi/(2*l)))
  
  t = numerator / denominator
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  print(power)
  return(power)
}

negativebinomial_purity = function(alpha,W,l,D,N,L,phi,F){
  n = sqrt(L/W)
  numerator = F*abs(N-2)*sqrt(L*D/l)
  denominator = 2*sqrt(1+W*D*phi/(l))
  
  t = numerator / denominator
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(power)
}

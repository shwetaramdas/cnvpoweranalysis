 #power calculator

##Function 1: Poisson with 100% purity
poisson_pure = function(alpha,W,l,D,N,L){
  n = L/W
  numerator = W*D*abs(N-2)*sqrt(L/W)/(2*l)
  denominator = sqrt(N*W*D/2/l)
  
  t = numerator / denominator
  print(t)
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(power)  
}


##Function 2: Poisson with < 100% purity: INCOMPLETE
poisson_impurity = function(alpha,W,l,D,N,L,F){
  n = sqrt(L/W)
  numerator = F*abs(N-2)*D*sqrt(L*W)
  
  #simulating variance 
  V = var(F*rpois(100000,lambda = (D*N*W)/(2*l)) + (1-F)*rpois(100000,lambda =D*W/l))
  denominator = 2*l*sqrt(V)
  
  t = numerator / denominator
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(power)
}

##Function 3: Negative binomial with 100% purity

##Function 3: Negative binomial with 100% purity
negativebinomial = function(alpha,W,l,D,N,L,muphi){
	#INPUT = diploid read depth
	D = D/2
	numerator = sqrt(D*L)*abs(N-2)
	denominator = sqrt(N*l*(1 + muphi))
	t = numerator / denominator

	##power
	df = round(L/W,0) - 1
	C = qt(1 - alpha,df=df)
	power = 1 - pt(C-t,df) + pt(-C-t,df)

	return(power)
}

##Function 4: Negative binomial with < 100% purity
negativebinomial_purity = function(alpha,W,l,D,N,L,muphi,F){
	#input is diploid read depth
	D = D/2
	n = sqrt(L/W)
	phi = muphi / (N*D/l)
	
	num = F*abs(N-2)*(D*W/l)*sqrt(L/W)
	V = var(F*rnbinom(1000000,size=1/phi, mu=(N*D*W/l)) + (1-F)*rnbinom(1000000,size=1/phi,mu=2*D*W/l))
	t = num / sqrt(V)

	##power
	df = round(L/W,0) - 1
	C = qt(1 - alpha,df=df)
	power = 1 - pt(C-t,df) + pt(-C-t,df)

	return(power)
}

#simulate variance of weighted sum of poissons
#var(F*rpois(10000,lambda = D*N*W/2/l) + (1-F)*rpois(10000,lambda =D*W/l))

##functions.R

#power calculator

##Function 1: Poisson with 100% purity
##Function 1: Poisson with 100% purity
poisson_pure = function(alpha,W,l,D,N,L){
  n = L/W
  numerator = W*D*abs(N-2)/l*sqrt(L/W)
  denominator = sqrt(N*W*D/l)
  
  t = numerator / denominator
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha/2,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)

  return(sprintf("%.2f",power))  
}

##Function 2: Poisson with < 100% purity: INCOMPLETE
poisson_impurity1 = function(alpha,W,l,D,N,L,F){
  n = sqrt(L/W)
  numerator = F*abs(N-2)*sqrt(D*L/l)
  denominator = sqrt(F*(N-2) + 2)
  t = numerator / denominator
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha/2,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(sprintf("%.2f",power))
}


##Function 3: Negative binomial with 100% purity
negativebinomial = function(alpha,W,l,D,N,L,muphi){
	#INPUT = diploid read depth
	numerator = sqrt(D*L)*abs(N-2)
	denominator = sqrt(N*l*(1 + muphi))
	t = numerator / denominator

	##power
	df = round(L/W,0) - 1
	C = qt(1 - alpha/2,df=df)
	power = 1 - pt(C-t,df) + pt(-C-t,df)

	return(sprintf("%.2f",power))
}


##Function 4: Negative binomial with < 100% purity

##Function 4: Negative binomial with < 100% purity
negativebinomial_purity = function(alpha,W,l,D,N,L,muphi,F){
  #input is diploid read depth
  n = sqrt(L/W)
  phi = muphi / (N*D/l)
  
  #	num = F*abs(N-2)*(D*W/l)*sqrt(L/W)
  #	V = var(F*rnbinom(1000000,size=1/phi, mu=(N*D*W/l)) + (1-F)*rnbinom(1000000,size=1/phi,mu=2*D*W/l))
  
  num = F*abs(N-2)*sqrt(D*L/l)
  den = sqrt(F*(N-2)+2) * sqrt(1+muphi)
  t = num / den
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha/2,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(sprintf("%.2f",power))
}



##Function 2: Poisson with < 100% purity: INCOMPLETE
poisson_highimpurity = function(alpha,W,l,D,N,L,F){
  n = sqrt(L/W)
  t = F * abs(N-2) * sqrt(D*L/(2*l))
  
  ##power
  df = round(L/W,0) - 1
  C = qt(1 - alpha/2,df=df)
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  
  return(sprintf("%.2f",power))
}


plot_power <- function(input){
  ps = c()
  if(input$tovary == 'None'){
    print('Select a parameter to vary to generate a power curve.')
  }
  else if((input$tovary=='Read Depth') && (input$tovary2=='None')){
    plot_readdepth(input);
  }else if((input$tovary == 'Read Depth')&&(input$tovary2 != 'None')){
    plot_double(input);
  }else if((input$tovary=='Read length') && (input$tovary2 == 'None')){
    plot_readlength(input);
  }else if((input$tovary=='Window size') && (input$tovary2 == 'None')){
    plot_windowsize(input);
  }else if((input$tovary=='Length of CNV') && (input$tovary2 == 'None')){
    plot_cnvlength(input);
  }else if((input$tovary=='Read length') && (input$tovary2 != 'None')){
    plot_double(input);
  }else if((input$tovary != "None") && (input$tovary2 != 'None')){
    plot_double(input);
  }
}

# functions.R

# D: haploid sequencing depth (for a diploid dataset of 30X total coverage, the haploid coverage is 15)
# L: the length of the CNA or CNV
# W: the size of the window
# N: the ploidy of the CNA or CNV (normal diploid regions have N=2)
# l: average sequencing read length
# F: the portion of the sample that contains the CNA or CNV
# alpha: the significance level
# theta: the variance inflation factor

# Function 0: Universal Applicable CPV Power Analysis
# A modefied version of function 4, which can be applied to all situation.
# This function will return 1 values, the power.
# This function has default value for alpha as 0.05, theta as 0, F as 1.
UA_CPA = function (alpha=0.05,W,l,D,N,L,theta=0,F=1){
  # calculate t score {
  n = L/W  # the sample size of the CNV
  numerator = F*abs(N-2)*sqrt(D*L)
  denominator = sqrt(F*(N-2)+2)*sqrt(l*(1+theta))
  t = numerator / denominator
  # }
  
  # calculate power {
  df = round(n,0) - 1  # define the degrees of freedom, df
  C = qt(1 - alpha/2,df=df)  # define the critical value, C
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  return(power)
  # }
}

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
poisson_impurity = function(alpha,W,l,D,N,L,F){
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

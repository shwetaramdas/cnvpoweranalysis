# power calculator

# D: haploid sequencing depth (for a diploid dataset of 30X total coverage, the haploid coverage is 15)
# L: the length of the CNA or CNV
# W: the size of the window
# N: the ploidy of the CNA or CNV (normal diploid regions have N=2)
# l: average sequencing read length
# F: the portion of the sample that contains the CNA or CNV
# alpha: the significance level
# theta: the variance inflation factor

# Function 0: Universal Applicable CPV Power Analysis
# A modefied version of function 4, which can be apply to all situation.
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

# Here are some examples:

# example: power = 0.3377715
UA_CPA (0.05,1000,100,0.1,3,10000)
# example: power = 0.9782693
UA_CPA (0.05,1000,100,0.1,3,50000)
# example: power = 0.9998702
UA_CPA (0.05,1000,100,0.1,3,100000)

# example: power = 0.2323759
UA_CPA (0.05,1000,100,0.1,3,10000,0.5)
# example: power = 0.9041381
UA_CPA (0.05,1000,100,0.1,3,50000,0.5)
# example: power = 0.9962494
UA_CPA (0.05,1000,100,0.1,3,100000,0.5)

# example: power = 0.124211
UA_CPA (0.05,1000,100,0.1,3,10000,0,0.5)
# example: power = 0.5891682
UA_CPA (0.05,1000,100,0.1,3,50000,0,0.5)
# example: power = 0.8792025
UA_CPA (0.05,1000,100,0.1,3,100000,0,0.5)


# example: power = 0.06612417
UA_CPA (0.05,1000,100,1,3,10000,1,0.1)
# example: power = 0.183031
UA_CPA (0.05,1000,100,1,3,50000,1,0.1)
# example: power = 0.3303402
UA_CPA (0.05,1000,100,1,3,100000,1,0.1)

# Function 1: Poisson sampling with 100% purity and no variance inflation
# refer to the section 4.1 of the paper
# This function will return 2 values, the t score and the power.
poisson_pure = function(alpha,W,l,D,N,L){
  # calculate t score {
  n = L/W  # the sample size of the CNV
  numerator = abs(N-2)*sqrt(D*L/l)
  denominator = sqrt(N)
  t = numerator / denominator  # equation 3b in the section 4.1 of the paper
  print(t)
  # }
  
  # calculate power {
  df = round(n,0) - 1  # define the degrees of freedom, df
  C = qt(1 - alpha/2,df=df)  # define the critical value, C
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  return(power)  
  # }
}

# example: t = 1.825742, power = 0.3377715
poisson_pure (0.05,1000,100,0.1,3,10000)
# example: t = 4.082483, power = 0.9782693
poisson_pure (0.05,1000,100,0.1,3,50000)
# example: t = 5.773503, power = 0.9998702
poisson_pure (0.05,1000,100,0.1,3,100000)

# Function 2: Poisson sampling with less than 100% purity and no variance inflation
# This function will return 2 values, the t score and the power.
# refer to the section 5.1 of the paper
poisson_impurity = function(alpha,W,l,D,N,L,F){
  # calculate t score {
  n = L/W  # the sample size of the CNV
  numerator = F*abs(N-2)*sqrt(D*L)
  denominator = sqrt(F*(N-2)+2)*sqrt(l)
  t = numerator / denominator  # equation 7 in the section 5.1 of the paper
  print(t)
  # }
  
  # calculate power {
  df = round(n,0) - 1  # define the degrees of freedom, df
  C = qt(1 - alpha/2,df=df)  # define the critical value, C
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  return(power)  
  # }
}

# example: t = 1, power = 0.124211
poisson_impurity (0.05,1000,100,0.1,3,10000,0.5)
# example: t = 2.236068, power = 0.5891682
poisson_impurity (0.05,1000,100,0.1,3,50000,0.5)
# example: t = 3.162278, power = 0.8792025
poisson_impurity (0.05,1000,100,0.1,3,100000,0.5)

# Function 3: Negative Binomial sampling with 100% purity and variance inflation
# refer to the section 4.2 of the paper
# This function will return 2 values, the t score and the power.
negativebinomial = function(alpha,W,l,D,N,L,theta){
  # calculate t score {
  n = L/W  # the sample size of the CNV
	numerator = sqrt(D*L)*abs(N-2)
	denominator = sqrt(N*l*(1 + theta))
	t = numerator / denominator	# equation 4b in the section 4.2 of the paper
	print(t)
	# }

	# calculate power {
	df = round(n,0) - 1	# define the degrees of freedom, df
	C = qt(1 - alpha/2,df=df)	# define the critical value, C
	power = 1 - pt(C-t,df) + pt(-C-t,df)
  return(power)
	# }
}

# example: t = 1.490712, power = 0.2323759
negativebinomial (0.05,1000,100,0.1,3,10000,0.5)
# example: t = 3.333333, power = 0.9041381
negativebinomial (0.05,1000,100,0.1,3,50000,0.5)
# example: t = 4.714045, power = 0.9962494
negativebinomial (0.05,1000,100,0.1,3,100000,0.5)

# Function 4: Negative binomial sampling with less than 100% purity and variance inflation
# refer to the section 5.2 of the paper
# This function will return 2 values, the t score and the power.
negativebinomial_purity = function(alpha,W,l,D,N,L,theta,F){
  # calculate t score {
  n = L/W  # the sample size of the CNV
  numerator = F*abs(N-2)*sqrt(D*L)
  denominator = sqrt(F*(N-2)+2)*sqrt(l*(1+theta))
  t = numerator / denominator  # equation 9 in the section 5.2 of the paper
  print(t)
  # }
  
  # calculate power {
  df = round(n,0) - 1  # define the degrees of freedom, df
  C = qt(1 - alpha/2,df=df)  # define the critical value, C
  power = 1 - pt(C-t,df) + pt(-C-t,df)
  return(power)  
  # }
}

# example: t = 2.236068, power = 0.4906241
negativebinomial_purity (0.05,1000,100,1,3,10000,1,0.5)
# example: t = 5, power = 0.9978252
negativebinomial_purity (0.05,1000,100,1,3,50000,1,0.5)
# example: t = 7.071068, power = 0.9999991
negativebinomial_purity (0.05,1000,100,1,3,100000,1,0.5)

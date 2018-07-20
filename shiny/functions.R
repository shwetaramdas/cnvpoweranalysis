# functions.R

# D: haploid sequencing depth (for a diploid dataset of 30X total coverage, the haploid coverage is 15)
# L: the length of the CNA or CNV
# W: the size of the window
# N: the ploidy of the CNA or CNV (normal diploid regions have N=2)
# l: average sequencing read length
# F: the portion of the sample that contains the CNA or CNV
# alpha: the significance level
# theta: the variance inflation factor

# Function 1: Universal Applicable CPV Power Analysis
# This function will return one value, the power.
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

# Function 2: Universal Aplicable CPV Power Analysis Plotter
# This function will return the power curve plot.
UA_CPAP = function(input){
  # Rename the input variables
  i_a=as.numeric(input$alpha)
  i_W=as.numeric(input$W)
  i_l=as.numeric(input$l)
  i_D=as.numeric(input$D)
  i_N=as.numeric(input$N)
  i_L=as.numeric(input$L)
  i_t=as.numeric(input$phi)
  i_F=as.numeric(input$F)
  # Create a list of x-varialbes which are used to generate the power curve
  x_values = seq(input$min1, input$max1, length.out = 1000)
  # Do nothing if None is selected for both
  if ((input$tovary=="None") && (input$tovary2=='None')){
  } else if (input$tovary2=='None') { # Plot a one-curve graph if only one variable is selete
    if (input$tovary=='Read Depth'){
      y_values = UA_CPA (i_a, i_W, i_l, x_values, i_N, i_L, i_t, i_F)
    } else if (input$tovary=='Read length'){
      y_values = UA_CPA (i_a, i_W, x_values, i_D, i_N, i_L, i_t, i_F)
    } else if (input$tovary=='Window size'){
      y_values = UA_CPA (i_a, x_values, i_l, i_D, i_N, i_L, i_t, i_F)
    } else if (input$tovary=='Length of CNV'){
      y_values = UA_CPA (i_a, i_W, i_l, i_D, i_N, x_values, i_t, i_F)
    } else if (input$tovary=='Sample purity'){
      y_values = UA_CPA (i_a, i_W, i_l, i_D, i_N, i_L, i_t, x_values)
    }
    df_1 <<- data.frame(a=x_values,b=y_values) # a global variable df_1 is generated so that user can download the raw data
    if (input$log_scale) {
      ggplot(df_1, aes(x=a,y=b)) + geom_line(colour="forestgreen") + xlab(input$tovary) + ylab("Power") + scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = "b") # generate the actually power plot with log scale for x-axis
    } else {
      ggplot(df_1, aes(x=a,y=b)) + geom_line(colour="forestgreen") + xlab(input$tovary) + ylab("Power")  # generate the actually power plot
    }
  } else {
    MV_list = c(input$val1, input$val2, input$val3) # a list that contians the 3 values for the second parameter to vary
    y_list = list()
    for (i in 1:3) { # assign the new parameter
      if (input$tovary2=='Read Depth'){
        i_D=MV_list[i]
      } else if (input$tovary2=='Read length'){
        i_l=MV_list[i]
      } else if (input$tovary2=='Window size'){
        i_W=MV_list[i]
      } else if (input$tovary2=='Length of CNV'){
        i_L=MV_list[i]
      } else if (input$tovary2=='Sample purity'){
        i_F=MV_list[i]
      }  
      if (input$tovary=='Read Depth'){ # calculate the power base on the new parameter
        y_list[[i]] = UA_CPA (i_a, i_W, i_l, x_values, i_N, i_L, i_t, i_F)
      } else if (input$tovary=='Read length'){
        y_list[[i]] = UA_CPA (i_a, i_W, x_values, i_D, i_N, i_L, i_t, i_F)
      } else if (input$tovary=='Window size'){
        y_list[[i]] = UA_CPA (i_a, x_values, i_l, i_D, i_N, i_L, i_t, i_F)
      } else if (input$tovary=='Length of CNV'){
        y_list[[i]] = UA_CPA (i_a, i_W, i_l, i_D, i_N, x_values, i_t, i_F)
      } else if (input$tovary=='Sample purity'){
        y_list[[i]] = UA_CPA (i_a, i_W, i_l, i_D, i_N, i_L, i_t, x_values)
      }  
    }
    df_1 <<- data.frame(a=x_values,b=y_list[[1]], c=y_list[[2]], d=y_list[[3]]) # a global variable df_1 is generated so that user can download the raw data
    data_cap = paste(" Green: ",input$tovary2,"=",input$val1,"\n","Blue: ",input$tovary2,"=",input$val2,"\n","Violet:",input$tovary2,"=",input$val3) # caption for the line color
    if (input$log_scale) {
      ggplot(df_1, aes(x=a)) + geom_line(aes(y=b), colour="forestgreen") + geom_line(aes(y=c), colour="blue") + geom_line(aes(y=d), colour="darkviolet") + xlab(input$tovary) + ylab("Power") + labs(subtitle=data_cap) + scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = "b") # generate the actually power plot with log sacle for x-axis
    } else {
      ggplot(df_1, aes(x=a)) + geom_line(aes(y=b), colour="forestgreen") + geom_line(aes(y=c), colour="blue") + geom_line(aes(y=d), colour="darkviolet") + xlab(input$tovary) + ylab("Power") + labs(subtitle=data_cap) # generate the actually power plot
    }
  }
}


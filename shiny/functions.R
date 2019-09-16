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
UA_CPA = function (alpha=0.05,W=1,l,D,N,L,theta=0,F=1,input){
  # Rename the input variables
  i_W=as.numeric(input$W)
  i_l=as.numeric(input$l)
  i_D=as.numeric(input$D)
  i_N=as.numeric(input$N)
  theta_t=theta
  # Define theta
  if (input$dis=="Poisson") {
    theta=0
  } else if (input$dis=="Negative Binomial with theta") {
    theta=as.numeric(input$theta_1)
  } else if (input$dis=="Negative Binomial with phi") {
    theta=as.numeric(input$phi_1)*i_N*i_D*i_W/i_l
  }
  if (input$tovary=='Overdispersion parameter' | input$tovary2=='Overdispersion parameter' | input$tovary3=='Overdispersion parameter' ) {
    theta=theta_t
  }
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
  if (input$dis=="Poisson") {
    i_t=0
  } else if (input$dis=="Negative Binomial with theta") {
    i_t=as.numeric(input$theta_1)
  } else if (input$dis=="Negative Binomial with phi") {
    i_t=as.numeric(input$phi_1)*i_N*i_D*i_W/i_l
  }
  i_F=as.numeric(input$F)
  x_values = seq(input$min1, input$max1, length.out = 1000) # Create a list of x-varialbes which are used to generate the power curve
  g <- ggplot() # 
  if ((input$tovary=="None") && (input$tovary2=='None') && (input$tovary3=='None')){
  } # Do nothing if all 3 variable is set to none
  else if ((input$tovary2=='None') && (input$tovary3=='None')) { # Plot a one-curve graph if only one variable is selete
    if (input$tovary=='Read Depth'){
      y_values = UA_CPA (i_a, i_W, i_l, x_values, i_N, i_L, i_t, i_F, input)
    } else if (input$tovary=='Read length'){
      y_values = UA_CPA (i_a, i_W, x_values, i_D, i_N, i_L, i_t, i_F, input)
    } else if (input$tovary=='Window size'){
      y_values = UA_CPA (i_a, x_values, i_l, i_D, i_N, i_L, i_t, i_F, input)
    } else if (input$tovary=='Length of CNV'){
      y_values = UA_CPA (i_a, i_W, i_l, i_D, i_N, x_values, i_t, i_F, input)
    } else if (input$tovary=='Sample purity'){
      y_values = UA_CPA (i_a, i_W, i_l, i_D, i_N, i_L, i_t, x_values, input)
    } else if (input$tovary=='Overdispersion parameter') {
      y_values = UA_CPA (i_a, i_W, i_l, i_D, i_N, i_L, theta=x_values, i_F, input)
    }
    df_1 <<- data.frame(a=x_values,b=y_values) # a global variable df_1 is generated so that user can download the raw data
    g <- ggplot(df_1, aes(x=a,y=b)) + geom_line(colour="forestgreen") + xlab(input$tovary) + ylab("Power")
  } 
  else if ((input$tovary3=='None')) {
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
      } else if (input$tovary2=='Overdispersion parameter'){
        i_t=MV_list[i]
      }  
      if (input$tovary=='Read Depth'){ # calculate the power base on the new parameter
        i_D = x_values
      } else if (input$tovary=='Read length'){
        i_l = x_values
      } else if (input$tovary=='Window size'){
        i_W = x_values
      } else if (input$tovary=='Length of CNV'){
        i_L = x_values
      } else if (input$tovary=='Sample purity'){
        i_F = x_values
      } else if (input$tovary=='Overdispersion parameter'){
        i_t = x_values
      }
      y_list[[i]] = UA_CPA (i_a, i_W, i_l, i_D, i_N, i_L, i_t, i_F, input)
    }
    df_1 <<- data.frame(a=x_values,b=y_list[[1]], c=y_list[[2]], d=y_list[[3]]) # a global variable df_1 is generated so that user can download the raw data
    data_cap = paste(" Color (",input$tovary2,"): ","Green=",input$val1,"   ","Blue=",input$val2,"   ","Violet=",input$val3,sep = "") # caption for the line color
    g = ggplot(df_1, aes(x=a)) + geom_line(aes(y=b), colour="forestgreen") + geom_line(aes(y=c), colour="blue") + geom_line(aes(y=d), colour="darkviolet") + xlab(input$tovary) + ylab("Power") + labs(subtitle=data_cap)
  }
  else  { # all three variable is seleted
    MV_list = c(input$val1, input$val2, input$val3) # a list that contians the 3 values for the second parameter to vary
    MV_list_3 = c(input$val31, input$val32, input$val33)
    y_list = list()
    line_color = list("forestgreen","blue","darkviolet")
    line_type = list("solid","longdash","dotted")
    g = ggplot(df_1, aes(x=a))
    for (i in 1:3) {
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
      } else if (input$tovary2=='Overdispersion parameter'){
        i_t=MV_list[i]
      }
      for (j in 1:3) {
        if (input$tovary3=='Read Depth'){
          i_D=MV_list_3[j]
        } else if (input$tovary3=='Read length'){
          i_l=MV_list_3[j]
        } else if (input$tovary3=='Window size'){
          i_W=MV_list_3[j]
        } else if (input$tovary3=='Length of CNV'){
          i_L=MV_list_3[j]
        } else if (input$tovary3=='Sample purity'){
          i_F=MV_list_3[j]
        } else if (input$tovary3=='Overdispersion parameter'){
          i_t=MV_list_3[j]
        }
        if (input$tovary=='Read Depth'){ # calculate the power base on the new parameter
          i_D = x_values
        } else if (input$tovary=='Read length'){
          i_l = x_values
        } else if (input$tovary=='Window size'){
          i_W = x_values
        } else if (input$tovary=='Length of CNV'){
          i_L = x_values
        } else if (input$tovary=='Sample purity'){
          i_F = x_values
        } else if (input$tovary=='Overdispersion parameter'){
          i_t = x_values
        }
        y_list[[(i-1)*3+j]] = UA_CPA (i_a, i_W, i_l, i_D, i_N, i_L, i_t, i_F, input)
        df_2 = data.frame(a=x_values, b=y_list[[(i-1)*3+j]])
        g = g + geom_line(data = df_2, aes(x= a, y=b), colour=toString(line_color[i]),linetype=toString(line_type[j]))
      }
      data_cap = paste(" Color (",input$tovary2,"): ","Green=",input$val1,"   ","Blue=",input$val2,"   ","Violet=",input$val3,"\n",
                       " Type (",input$tovary3,"): ","Solid=",input$val31,"   ","Dash=",input$val32,"   ","Dot=",input$val33,
                       sep = "") # caption for the line color
      g = g + xlab(input$tovary) + ylab("Power") + labs(subtitle=data_cap)
    }
  }

  if (input$fix_y_range) { # option for fixing y-axis
    g = g + scale_y_continuous(limits = c(0,1), breaks=seq(0, 1, 0.1))
  }
  
  if (input$show_x_at_0 & !input$log_scale) { # option for showing x at 0 for non-log-scale x-axis
    g = g + scale_x_continuous(limits = c(0,NA))
  } 
  
  if (input$log_scale) { # option for using log-scale for x-axis
    g = g + scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = "b") # generate the actually power plot with log sacle for x-axis
  }
  
  g
}


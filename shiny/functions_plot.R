#plot functions
source("functions.R")

#single parameter: read depth
plot_readdepth = function(input){
  Ds = seq(input$min1, input$max1, length.out = 20);
  ps = c();
  for(d in Ds){
    #poisson or negative binomial
    if(input$dis == "pois"){
      if(input$F == 1){ 
        ps = c(ps,poisson_pure(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), d,as.numeric(input$N), as.numeric(input$L)))
      }else{
        ps = c(ps,poisson_impurity(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), d,as.numeric(input$N), as.numeric(input$L),input$F)) 
      }
    }else{
      if(is.na(input$phi)){phi = 1/(input$W*input$D/input$l);}
      else{phi = input$phi;}
      if(input$F == 1){
        ps = c(ps,negativebinomial_purity(input$alpha,input$W,input$l,d,input$N,input$L,phi,input$F))
      }else{
        ps = c(ps, negativebinomial_purity(input$alpha,input$W,input$l,d,input$N,input$L,phi,input$F))
      }
    }  
  }
  toplot = data.frame(a=Ds,b=ps)
  ggplot(toplot, aes(x=a,y=b)) + geom_point() + geom_line() + xlab("Read depth") + ylab("Power") 
  #plot(Ds, ps,type="l",xlab="Read Depth",ylab="Power")
#  points(Ds, ps,pch=19);
}

#single parameter: read length
plot_readlength = function(input){
  ls = seq(input$min1, input$max1, length.out = 20);
  ps = c()
  for(l in ls){
    if(input$dis == "pois"){
      if(input$F == 1){
        ps = c(ps,poisson_pure(as.numeric(input$alpha), as.numeric(input$W), l, as.numeric(input$D),as.numeric(input$N), as.numeric(input$L)))
      }else{
        ps = c(ps,poisson_impurity(input$alpha, input$W, l, input$D,input$N, input$L,input$F))
      }
    }else{
      if(is.na(input$phi)){phi = 1/(input$W*input$D/input$l);}
      else{phi = input$phi;}  
      if(input$F == 1){
            ps = c(ps,negativebinomial(input$alpha, input$W, l, input$D, input$N, phi, input$L))
      }else{
            ps = c(ps,negativebinomial_purity(input$alpha, input$W, l, input$D, input$N, input$L,phi,input$F))
      }
    }
  }
   toplot = data.frame(a=ls,b=ps)
   ggplot(toplot, aes(x=a,y=b)) + geom_point() + geom_line() + xlab("Read length") + ylab("Power") 
#  plot(ls, ps,type="l",xlab="Read Length",ylab="Power")
#  points(ls, ps,pch=19);
}


#single parameter: read length
plot_windowsize = function(input){
  Ws = seq(input$min1, input$max1, length.out = 20);
  ps = c()
  for(W in Ws){
    if(input$dis == "pois"){
      if(input$F == 1){
        ps = c(ps,poisson_pure(input$alpha, W, input$l, input$D, input$N, input$L))
      }else{
        ps = c(ps,poisson_impurity(input$alpha, W, input$l, input$D,input$N, input$L,input$F))
      }
    }else{
      if(is.na(input$phi)){phi = 1/(input$W*input$D/input$l);}
      else{phi = input$phi;}  
      if(input$F == 1){
        ps = c(ps,negativebinomial(input$alpha, W, input$l, input$D, input$N, phi, input$L))
      }else{
        ps = c(ps,negativebinomial_purity(input$alpha, W, input$l, input$D, input$N, input$L,phi,input$F))
      }
    }
  }
  toplot = data.frame(a=Ws,b=ps)
  output$powercurvevalues = toplot
  ggplot(toplot, aes(x=a,y=b)) + geom_point() + geom_line() + xlab("Window size") + ylab("Power") 
  #  plot(ls, ps,type="l",xlab="Read Length",ylab="Power")
  #  points(ls, ps,pch=19);
}

#single parameter: read length
plot_cnvlength = function(input){
  Ls = seq(input$min1, input$max1, length.out = 20);
  ps = c()
  for(L in Ls){
    if(input$dis == "pois"){
      if(input$F == 1){
        ps = c(ps,poisson_pure(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D),as.numeric(input$N), L))
      }else{
        ps = c(ps,poisson_impurity(input$alpha, input$W, as.numeric(input$l), input$D,input$N, L,input$F))
      }
    }else{
      if(is.na(input$phi)){phi = 1/(input$W*input$D/input$l);}
      else{phi = input$phi;}  
      if(input$F == 1){
            ps = c(ps,negativebinomial(input$alpha, input$W, input$l, input$D, input$N, phi, L))
      }else{
            ps = c(ps,negativebinomial_purity(input$alpha, input$W, input$l, input$D, input$N, L,phi,input$F))
      }
    }
  }
   toplot = data.frame(a=Ls,b=ps)
   ggplot(toplot, aes(x=a,y=b)) + geom_point() + geom_line() + xlab("Length of CNV") + ylab("Power")
}

#single parameter: read length
plot_samplepurity = function(input){
  Fs = seq(input$min1, input$max1, length.out = 20);
  ps = c()
  for(F in Fs){
    if(input$dis == "pois"){
        ps = c(ps,poisson_impurity(input$alpha, input$W, as.numeric(input$l), input$D,input$N, input$L,F))
    }else{
      if(is.na(input$phi)){phi = 1/(input$W*input$D/input$l);}
      else{phi = input$phi;}  
        ps = c(ps,negativebinomial_purity(input$alpha, input$W, input$l, input$D, input$N, input$L,phi,F))
    }
  }
   toplot = data.frame(a=Fs,b=ps)
   ggplot(toplot, aes(x=a,y=b)) + geom_point() + geom_line() + xlab("Sample Purity") + ylab("Power")
}



plot_double <- function(input){
  COLOURS = c('black','red','orange','purple','blue','green','dark green')
  V1s = seq(input$min1, input$max1, length.out = 100);
  V2s = c(input$val1, input$val2, input$val3);
  
  V1 = input$tovary
  V2 = input$tovary2
  
  Ps = list()
  for(n1 in 1:length(V2s)){
	print(n1)
    v2 = V2s[n1]
    ps = c()
    for(i in 1:length(V1s)){
      ps = c(ps, getpower(input,V1s[i],v2,V1,V2))
    }
    Ps[[n1]] = ps
  }

  if(V1 == "Length of CNV"){V1s = log10(V1s);}
  toplot = data.frame(a=V1s, b=Ps[[1]],Parameter_2=V2s[1],col='black')
  
  for(i in 2:length(Ps)){
    temp = data.frame(a=V1s,b=Ps[[i]],Parameter_2=V2s[i],col=COLOURS[i])
    toplot = rbind(toplot, temp)
  }
  
  ggplot(toplot,aes(x=a,y=b,colour=as.factor(Parameter_2))) + geom_point() + geom_line() + xlab(V1) + ylab("Power") + labs(color=V2) 
}


getpower <- function(input,v1,v2,V1,V2){
  if(input$dis == 'pois'){
    if((V1=='Read Depth')&&(V2=='Read length')){return(poisson_impurity(input$alpha,input$W,v2,v1,input$N,input$L,input$F));}
    else if((V1=='Read Depth')&&(V2=='Window size')){return(poisson_impurity(input$alpha,v2,input$l,v1,input$N,input$L,input$F));}
    else if((V1=='Read Depth')&&(V2=='Length of CNV')){return(poisson_impurity(input$alpha,input$W,input$l,v1,input$N,v2,input$F));}
    else if((V1=='Read Depth')&&(V2=='Sample purity')){return(poisson_impurity(input$alpha,input$W,input$l,v1,input$N,input$L,v2));}
	
    else if((V1=='Read length')&&(V2=='Read Depth')){return(poisson_impurity(input$alpha,input$W,v1,v2,input$N,input$L,input$F));}
    else if((V1=='Read length')&&(V2=='Length of CNV')){return(poisson_impurity(input$alpha,input$W,v1,input$D,input$N,v2,input$F));}
    else if((V1=='Read length')&&(V2=='Window size')){return(poisson_impurity(input$alpha,v2,v1,input$D,input$N,input$L,input$F));}
    else if((V1=='Read length')&&(V2=='Sample purity')){return(poisson_impurity(input$alpha,input$W,v1,input$D,input$N,input$L,v2));}
	
	else if((V1=='Length of CNV')&&(V2=='Read Depth')){return(poisson_impurity(input$alpha,input$W,input$l,v2,input$N,v1,input$F));}
	else if((V1=='Length of CNV')&&(V2=='Read length')){return(poisson_impurity(input$alpha,input$W,v2,input$D,input$N,v1,input$F));}
	else if((V1=='Length of CNV')&&(V2=='Window size')){return(poisson_impurity(input$alpha,v2,input$l,input$D,input$N,v1,input$F));}
	else if((V1=='Length of CNV')&&(V2=='Sample purity')){return(poisson_impurity(input$alpha,input$W,input$l,input$D,input$N,v1,v2));}

	else if((V1=='Sample purity')&&(V2=='Read Depth')){return(poisson_impurity(input$alpha,input$W,input$l,v2,input$N,input$L,v1));}
	else if((V1=='Sample purity')&&(V2=='Read length')){return(poisson_impurity(input$alpha,input$W,v2,input$D,input$N,input$L,v1));}
	else if((V1=='Sample purity')&&(V2=='Window size')){return(poisson_impurity(input$alpha,v2,input$l,input$D,input$N,input$L,v1));}
	else if((V1=='Sample purity')&&(V2=='Length of CNV')){return(poisson_impurity(input$alpha,input$W,input$l,input$D,input$N,v2,v1));}
	
  }else{
    if(is.na(input$phi)) {
      phi = 1/(input$W*input$D/input$l)
    }else{
      phi = input$phi
    }
#    print(paste(V1,V2))
    if((V1=='Read Depth')&&(V2=='Read length')){return(negativebinomial_purity(input$alpha,input$W,v2,v1,input$N,input$L,phi,input$F));}
    else if((V1=='Read Depth')&&(V2=='Window size')){return(negativebinomial_purity(input$alpha,v2,input$l,v1,input$N,input$L,phi,input$F));}
    else if((V1=='Read Depth')&&(V2=='Length of CNV')){return(negativebinomial_purity(input$alpha,input$W,input$l,v1,input$N,v2,phi,input$F));}
    else if((V1=='Read Depth')&&(V2=='Sample purity')){return(negativebinomial_purity(input$alpha,input$W,input$l,v1,input$N,input$L,phi,v2));}    
	
    else if((V1=='Read length')&&(V2=='Read Depth')){return(negativebinomial_purity(input$alpha,input$W,v1,v2,input$N,input$L,phi,input$F));}
    else if((V1=='Read length')&&(V2=='Length of CNV')){return(negativebinomial_purity(input$alpha,input$W,v1,input$D,input$N,v2,phi,input$F));}
    else if((V1=='Read length')&&(V2=='Window size')){return(negativebinomial_purity(input$alpha,v2,v1,input$D,input$N,input$L,phi,input$F));}
    else if((V1=='Read length')&&(V2=='Sample purity')){return(negativebinomial_purity(input$alpha,input$W,v1,input$D,input$N,input$L,phi,input$F));}
	
	else if((V1=='Length of CNV')&&(V2=='Read Depth')){return(negativebinomial_purity(input$alpha,input$W,input$l,v2,input$N,v1,phi,input$F));}
	else if((V1=='Length of CNV')&&(V2=='Read length')){return(negativebinomial_purity(input$alpha,input$W,v2,input$D,input$N,v1,phi,input$F));}
	else if((V1=='Length of CNV')&&(V2=='Window size')){return(negativebinomial_impurity(input$alpha,v2,input$l,input$D,input$N,v1,phi,input$F));}
	else if((V1=='Length of CNV')&&(V2=='Sample purity')){return(negativebinomial_impurity(input$alpha,input$W,input$l,input$D,input$N,v1,phi,v2));}

	else if((V1=='Sample purity')&&(V2=='Read Depth')){return(negativebinomial_purity(input$alpha,input$W,input$l,v2,input$N,input$L,phi,v1));}
	else if((V1=='Sample purity')&&(V2=='Read length')){return(negativebinomial_purity(input$alpha,input$W,v2,input$D,input$N,input$L,phi,v1));}
	else if((V1=='Sample purity')&&(V2=='Window size')){return(negativebinomial_impurity(input$alpha,v2,input$l,input$D,input$N,input$L,phi,v1));}
	else if((V1=='Sample purity')&&(V2=='Length of CNV')){return(negativebinomial_impurity(input$alpha,input$W,input$l,input$D,input$N,v2,phi,v1));}
	
  }
}

powervalues <- function(input){
  if(input$tovary2 == "None"){
	V1s = seq(input$min1, input$max1, length.out = 20);
	V1 = input$tovary;
	  ps = c()
	  for(i in 1:length(V1s)){
	     ps = c(ps, poisson_impurity(as.numeric(input$alpha), as.numeric(input$W), as.numeric(input$l), as.numeric(input$D), as.numeric(input$N), as.numeric(input$L),as.numeric(input$F)))
	  }
	  toplot = data.frame(a=V1s, b=ps)
	  colnames(toplot)[1] = input$tovary
	  colnames(toplot)[2] = 'Power'
  }else{
	  V1 = input$tovary
	  V2 = input$tovary2
	  V1s = seq(input$min1, input$max1, length.out = 20);
	  V2s = c(input$val1, input$val2, input$val3);
	  
	  Ps = list()
	  toplot = data.frame()
	  for(n1 in 1:length(V2s)){
		v2 = V2s[n1]
		ps = c()
		for(i in 1:length(V1s)){
		  ps = c(ps, getpower(input,V1s[i],v2,V1,V2))
		}
		Ps[[n1]] = ps
		toplot1 = data.frame(a=V1s, b=ps,Parameter_2=rep(v2,length(ps)))
		toplot = rbind(toplot, toplot1)
	  }
	   colnames(toplot)[1] = input$tovary
	   colnames(toplot)[2] = 'Power'
	   colnames(toplot)[3] = input$tovary2
	   
	}
	
	return(toplot)
}
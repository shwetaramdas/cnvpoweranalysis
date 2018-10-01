#experimental parameters
W = 1000
l=100
L=10000
D=30
N=3
mu=N*D/2


#Negative binomial simulations
#What is effect of overdispersion on variance of aggregated read counts, for different numbers of events 
Vs = c()
SIZES= c(mu*mu/0.5, mu*mu/1, mu*mu/2, mu*mu/3)
for(phis in SIZES){
  df2 = data.frame(nrow=1000,ncol=1)
  for(i in 1:100){
    df2 = cbind(df2, rnbinom(1000,size=phis,mu=mu))
  }
  df2 = df2[,-c(1:2)]
  Vs = c(Vs,var(apply(df2,2,mean)))
}
print(Vs)
V1s = Vs


Vs = c()
for(phis in SIZES){
  df2 = data.frame(nrow=10000,ncol=1)
  for(i in 1:1000){
    df2 = cbind(df2, rnbinom(10000,size=phis,mu=mu))
  }
  df2 = df2[,-c(1:2)]
  Vs = c(Vs,var(apply(df2,2,mean)))
}
print(Vs)
V2s = Vs


Vs = c()
for(phis in SIZES){
  df2 = data.frame(nrow=1000,ncol=1)
  for(i in 1:10000){
    df2 = cbind(df2, rnbinom(1000,size=phis,mu=mu))
  }
  df2 = df2[,-c(1:2)]
  Vs = c(Vs,var(apply(df2,2,mean)))
}
print(Vs)

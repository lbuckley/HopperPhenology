#This code should be run after the model is fit to the data using the code in
#"Wrapper R code_salamander application.R".  See that file for more information.

#Calculate the population growth rate, lambda, over the time period of the study
#Do this using the full posterior distribution
#Then summarize using the mean and the 95% credible interval

# ---------------------------------------------
# function to calculate lambda, 
# population growth, for given matrix
# Courtesy of Evan Cooch
# ---------------------------------------------

calc_lam=function(mat) {
  va=eigen(mat)$values;
  rmax=max(Re(va));
  j=which(Re(va)==rmax);
  lambda=Re(va[j]);
  return(lambda)
}

#########################################

#Pull out the relevant parameters from the workspace
chains=as.matrix(StageSample1)

lambda1=chains[,"lambda[1]"]
lambda2=chains[,"lambda[2]"]
lambda3=chains[,"lambda[3]"]

gamma1=chains[,"gamma[1]"]
gamma2=chains[,"gamma[2]"]
gamma3=chains[,"gamma[3]"]

omega1=chains[,"omega[1]"]
omega2=chains[,"omega[2]"]


#Examine the population growth rate under the assumption of a Leslie matrix
pop.growth=rep(0,length(gamma1))

for (i in 1:length(gamma1)) {
a=matrix(0,nrow=3,ncol=3)
a[1,] = c( 0,          0,         gamma1[i])
a[2,] = c( omega1[i],  0,         0)
a[3,] = c( 0,          omega1[i], omega2[i])

pop.growth[i] = calc_lam(a)
}

summary(pop.growth)
hist(pop.growth)
median(pop.growth)
quantile(pop.growth,c(0.025,0.975))

#Examine the population size for the 21 sampled locations in 50 years
pop.size=matrix(0,nrow=length(gamma1), ncol=5)
colnames(pop.size)=c("year 8", "year 18", "year 28", "year 38", "year 48")

for (i in 1:length(gamma1)) {
  junk=array(0,dim=c(4,48,21))
  junk[1,1,]=rpois(21,lambda1[i])
  junk[2,1,]=rpois(21,lambda2[i])
  junk[3,1,]=rpois(21,lambda3[i])
  junk[4,1,]=apply(junk[1:3,1,],2,sum)
  for (t in 2:48) {
    junk[1,t,] = rpois(21,(junk[3,t-1,]*gamma1[i] + gamma2[i]))
    junk[2,t,]=  rbinom(21,junk[1,t-1,],omega1[i]) + rpois(21,gamma2[i])
    junk[3,t,]= (rbinom(21,junk[2,t-1,],omega1[i]) + rbinom(21,junk[3,t-1,],omega2[i]) 
                 + rpois(21,gamma3[i]) )
    junk[4,t,]=apply(junk[1:3,t,],2,sum)
  }
 a=apply(junk[4,,],1,sum)  
 pop.size[i,]=c(a[8],a[18],a[28],a[38],a[48])
}

apply(pop.size,2,summary)
apply(pop.size,2,quantile,c(0.025,0.5,0.975))



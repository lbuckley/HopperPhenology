#This is the JAGS code file to run the structured Dail-Madsen model
#in Zipkin et al. 2014 Ecology.  See "Wrapper R code_salamander application.R" for 
#more information and instructions on how to run the code.


model {
  
  # Specify the priors for all parameters in the model
  lambda[1] ~ dunif(0, 50)
  lambda[2] ~ dunif(0, 50)
  lambda[3] ~ dunif(0, 50)
  
  gamma[1] ~ dunif(0, 50)
  gamma[2] ~ dunif(0, 50)
  gamma[3] ~ dunif(0, 50)
  
  omega[1] ~ dunif(0, 1)
  omega[2] ~ dunif(0, 1)
  
  p[1] ~ dunif(0, 1)
  p[2] ~ dunif(0, 1)
 
  #Create a loop across all j sites
  for(j in 1:nSites) {
    
    
    #Intitate the model for year 1 - poisson with parameter lambda
    #The abundance matrix N is specified location x year x stage
    N[j,1,1] ~ dpois(lambda[1])
    N[j,1,2] ~ dpois(lambda[2])
    N[j,1,3] ~ dpois(lambda[3])
    
    #Create S and G vectors for year 1, which are not used in the model
    #S and G for year one are set to be consistent with the intial values
    for (r in 1:3) {
    S[j,1,r] ~ dpois(2)  
    G[j,1,r] ~ dpois(20)
    }
    
    #Specify the model for years 2 through nYears 
    for(t in 2:nYears) {
      
      #Estimate survivorship
      S[j,t,1] ~ dbin(omega[1], N[j,t-1,1]) 
      S[j,t,2] ~ dbin(omega[1], N[j,t-1,2])
      S[j,t,3] ~ dbin(omega[2], N[j,t-1,3])
      
      #Estimate recruitment (gamma1) and movement (gamma2 and gamma3)
      G[j,t,1] ~ dpois(gamma[1]*N[j,t-1,3] + gamma[2])
      G[j,t,2] ~ dpois(gamma[2])
      G[j,t,3] ~ dpois(gamma[3])
      
      #Sum all stages to get total N at each site j in each year t
      N[j,t,1] <- G[j,t,1]
      N[j,t,2] <- S[j,t,1] + G[j,t,2]
      N[j,t,3] <- S[j,t,2] + S[j,t,3] + G[j,t,3]
    }
      
      #Loop accross reps to estimate detection probability for all years
      #The data matrix n is specified location x year x stage x rep
    for (t in 1:nYears){
      for(k in 1:nReps){
        n[j,t,1,k] ~ dbin(p[1], (N[j,t,1]+N[j,t,2])) #Stages S1 and S2 are indistinguishable
        n[j,t,2,k] ~ dbin(p[2], N[j,t,3]) #Detection probability is the same for all adults
      }  }
      
    }
  
  
  #sum up the number of individuals in all locations to estimate annual
  #total N for each stage
  for (t in 1:nYears){
    Ntotal[1,t] <- sum(N[,t,1])
    Ntotal[2,t] <- sum(N[,t,2])
    Ntotal[3,t] <- sum(N[,t,3])
  }
    
}
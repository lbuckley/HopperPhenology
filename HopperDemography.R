#This code is associated with Zipkin et al. 2014 Ecology
#The code loads and formats the spotted salamander data and then fits a structured
#Dail-Madsen model using the associated JAGS code "JAGS code_salamander application.R".
#The code produces a summary of the parameter estimates for the model and model diagnostics
#(e.g., traceplots and R-hat statisic)

#The "salamander_leslie matrix.R" code estimates the population growth rate, lambda,
#using the full posterior distribution.

#This code runs with R version 3.2.2


#Load the library rjags
library(rjags)

#Data - desired format: (site, year, stage, rep)
#Read in data files
#Note that the site names have been changed as per the request of the data providers
ads<- read.table("saldata.aggregated_A.csv", header=TRUE,sep=",",na.strings=TRUE)
juv<- read.table("saldata.aggregated_J.csv", header=TRUE,sep=",",na.strings=TRUE)

#Format the data to the necessary structure
#Adults
junads<-ads[,c(2,4,6,8,10,12,14)]
julads<-ads[,c(3,5,7,9,11,13,15)]
#Juveniles
junjuv<-juv[,c(2,4,6,8,10,12,14)]
juljuv<-juv[,c(3,5,7,9,11,13,15)]

nSites=length(ads$Stream) #Total number of locations
nYears=dim(junads)[2]     #Total number of survey years
nObsAges=2                #Number of observable stages
nReps=2                   #Number of replicate surveys 

#Create an observation matrix n
n <- array(NA, dim=c(nSites,nYears,nObsAges,nReps), 
           dimnames=list(ads$Stream,2005:2011,c("Juveniles","Adults"),c("Jun","Jul")) )       

#Fill in the adult data to n observation matrix
for (j in 1:nSites){
  for (t in 1:nYears){
    n[j,t,2,1]<-junads[j,t]
    n[j,t,2,2]<-julads[j,t]
  }
}

#Fill in the juvenile data to n observation matrix
for (j in 1:nSites){
  for (t in 1:nYears){
    n[j,t,1,1]<-junjuv[j,t]
    n[j,t,1,2]<-juljuv[j,t]
  }
}

#Set the intital values to run the JAGS model
lamNew=NA; omegaNew=NA; gammaNew=NA; pNew=NA;
lamNew[1] <- 5 * 10          #Initial population abundance (juveniles)
lamNew[2] <- 5 * 10          #Initial population abundance (adults)
omegaNew[1] <- 0.90          #Survival rate (juveniles)
omegaNew[2] <- 0.90          #Survival rate (adults)
gammaNew[1] <- 5 * 10        #Recruitment rate 
gammaNew[2] <- 5 * 10        #Movement rate (juveniles)
gammaNew[3] <- 5 * 10        #Movement rate (adults)

pNew[1] <- 0.8              # Detection probability
pNew[2] <- 0.8

#Must specify intital values for N,G,S that are consistent with the model
#NNew and GNew are empty arrays
#Fix SNew to some value greater than zero
NNew <- array(NA, dim=c(nSites,nYears,3), dimnames=list(paste("Site",1:nSites),
                                                        paste("Year",1:nYears),c("S1", "S2","S3")) )

GNew <- array(NA, dim=c(nSites,nYears,3), 
              dimnames=list(paste("Site",1:nSites),paste("Year",1:(nYears)),
                            c("S1", "S2","S3")) )     # S = Survivors; G=NewRecruits
SNew <- array(2, dim=c(nSites,nYears,3), 
              dimnames=list(paste("Site",1:nSites),paste("Year",1:(nYears)),
                            c("S1", "S2","S3")) )

#Find the max observed number of individuals for each site
nmax<- apply(n,c(1,2,3),max)

#TransMat not necessary for this code but it does show how the population
#transitions into other stages
TransMat <- array( c(0,1,0,  0,0,1,  0,0,1), dim=c(3,3), 
                   dimnames=list(c("S1_t+1", "S2_t+1","S3_t+1"),
                                 c("S1_t","S2_t","S3_t")))

#Add a big number to nmax for each stage to fill in the NNew matrix
NNew[,1,1] <- nmax[,1,1] +15
NNew[,1,2] <- nmax[,1,1] +15
NNew[,1,3] <- nmax[,1,2] +15

#Very important!!
#Make sure SNew+GNew=NNew or jags will not run!
GNew = NNew-SNew

#Create all the necessary inputs for JAGS
#Bundle data
Dat <- list(nSites=nSites, nYears=nYears, nReps=nReps, n=n)

#Set intial values
InitStage <- function() list(gamma=gammaNew, N=NNew, 
                             omega=omegaNew, G=GNew, p=pNew, S=SNew) 

# Parameters to be monitored
ParsStage <- c("lambda", "gamma", "omega", "p", 
               "Ntotal")

# Sequential - start the adaptive phase of the model 
StageBurnin1 <- jags.model(paste("JAGS code_salamander application.r",sep=""), 
                           Dat, InitStage, n.chains=3, n.adapt=10000)


#Keep 10000/5 samples after the initial burn in of 10000 (above)
Niter=50000
StageSample1 <- coda.samples(StageBurnin1, ParsStage, n.iter=Niter, thin=10)

#Print out a summary of the parameter estimates
summary(StageSample1)

#Graph the model results
plot(StageSample1)

#Check the R-hat statistic to ensure convergence
library(coda)
gelman.diag(StageSample1)

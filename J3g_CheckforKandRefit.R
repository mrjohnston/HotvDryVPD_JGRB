#Checking all the models for k value; re-fitting those that need it!
# --> note that if edf << k', increasing k doesn't seem to do anything... don't refit
# refit if: p<0.05 and edf is within THREE units of k' [I just picked 3 to be ~ conservative]
#NO air models need refit according to this criteria

setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/") #AIR
setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/") #LEAF

gppmods<-list.files(pattern="GPP_D_*")
neemods<-list.files(pattern="NEE_D_*")
lemods<-list.files(pattern="LE_D_*")

for(i in 1:length(gppmods)){
  print(paste0("site (i): ",i))
  gpp<-readRDS(gppmods[i])
  nee<-readRDS(neemods[i])
  le<-readRDS(lemods[i])
  for(e in 1:length(gpp)){
    #print(paste("GPP, mod:",e,"FILE: ",gppmods[i]))
    r <- k.check(gpp[[e]])
    tocheck<-sum(r[,"k'"]-r[,'edf']<3)
    tocheckp<-sum(r[,"p-value"]<0.05)
    if(tocheck>0 & tocheckp>0) {print(paste("CHECK model: ",gppmods[i],"; model: ",e))}
    #print(paste("NEE, mod:",e,"FILE: ",neemods[i]))
    r <- k.check(nee[[e]])
    tocheck<-sum(r[,"k'"]-r[,'edf']<3)
    tocheckp<-sum(r[,"p-value"]<0.05)
    if(tocheck>0 & tocheckp>0) {print(paste("CHECK model: ",neemods[i],"; model: ",e))}
    #print(paste("LE, mod:",e,"FILE: ",lemods[i]))
    r <- k.check(le[[e]])
    tocheck<-sum(r[,"k'"]-r[,'edf']<3)
    tocheckp<-sum(r[,"p-value"]<0.05)
    if(tocheck>0 & tocheckp>0) {print(paste("CHECK model: ",lemods[i],"; model: ",e))}
  }
}

#"CHECK model:  LE_D_Air_Vcp.rds ; model:  6"
setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/") #AIR
lemods<-list.files(pattern="LE_D_*")
le<-readRDS(lemods[6])
gam.check(le[[6]]) #not a k problem

#Note, from Simon Wood presentation:
#"pattern may also be caused by
#mean-variance problems, missing covariates, structural
#infelicities, zero inflation etc. . ."


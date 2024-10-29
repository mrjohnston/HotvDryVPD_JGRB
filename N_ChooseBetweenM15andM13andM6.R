#Only selecting between M15 and M13 AND gpp[[6]] (3-way)
  #M13 is 14th in the list (i.e. gpp[[14]])
  #M15 is last in the list (i.e. gpp[[16]])
      #check the calls to be sure: gpp[[14]]$call
library(mgcv)

choosemod3<-function(lmods){ #list of models
  modsummary<-AIC(lmods[[6]],lmods[[14]],lmods[[16]])
  ord<-order(modsummary$AIC) #first is minimum
  #if min AIC model has AIC at least 2 units lower than the next lowest,
  #choose regardless of degrees of freedom
  if(modsummary$AIC[ord[1]]-modsummary$AIC[ord[2]] < (-2)){ 
    return(ord[1])
    stop("Done. lowest AIC has no competition")
  } 
  else {
    #which models are within 2 AIC units, including self?
    within2<-which(modsummary$AIC[ord[1]]-modsummary$AIC >=(-2))   
    #corresponding degrees of freedom:
    theirDF<-modsummary$df[within2]
    #pick the lowest degrees of freedom (aka fewest parameters)
    bestmod<-within2[which(theirDF==min(theirDF))]
    return(bestmod)
  }
}

#AIR  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/") #DAILY outliers removed
gppmods<-list.files(pattern="GPP_D_*")
neemods<-list.files(pattern="NEE_D_*")
lemods<-list.files(pattern="LE_D_*")

airresults<-data.frame(matrix(NA,nrow=26,ncol=7))
names(airresults)<-c("site","GPPchoice","GPPdevexpl",
                     "NEEchoice","NEEdevexpl",
                     "LEchoice","LEdevexpl")

for(i in 1:length(gppmods)){
  
  site<-substr(gppmods[i],nchar(gppmods[i])-6,nchar(gppmods[i])-4)
  print(paste(site,i))
  airresults[i,1]<-site
  gpp<-readRDS(gppmods[i])
  choose<-choosemod3(gpp)
  
  airresults[i,2]<-choose 
  if(choose==1){ #3-way interaction
    airresults[i,3]<-summary(gpp[[6]])$dev.expl
  }
  if(choose==2){ #no interaction 
    airresults[i,3]<-summary(gpp[[14]])$dev.expl
  }
  if(choose==3){ #2-way interaction
    airresults[i,3]<-summary(gpp[[16]])$dev.expl
  }
  
  nee<-readRDS(neemods[i])
  choose<-choosemod3(nee)
  airresults[i,4]<-choose
  if(choose==1){ #3-way interaction
    airresults[i,5]<-summary(nee[[6]])$dev.expl
  }
  if(choose==2){
    airresults[i,5]<-summary(nee[[14]])$dev.expl
  }
  if(choose==3){
    airresults[i,5]<-summary(nee[[16]])$dev.expl
  }
  
  le<-readRDS(lemods[i])
  choose<-choosemod3(le)
  airresults[i,6]<-choose
  if(choose==1){ #3-way interaction
    airresults[i,7]<-summary(le[[6]])$dev.expl
  }
  if(choose==2){
    airresults[i,7]<-summary(le[[14]])$dev.expl
  }
  if(choose==3){
    airresults[i,7]<-summary(le[[16]])$dev.expl
  }
}

table(airresults[,2])
table(airresults[,4])
table(airresults[,6])
write.csv(airresults,"/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/airresults_ALL_modsummary_with3way.csv",row.names = FALSE)


#LEAF  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/") #DAILY outliers removed
gppmods<-list.files(pattern="GPP_D_*")
neemods<-list.files(pattern="NEE_D_*")
lemods<-list.files(pattern="LE_D_*")

leafresults<-data.frame(matrix(NA,nrow=26,ncol=7))
names(leafresults)<-c("site","GPPchoice","GPPdevexpl",
                     "NEEchoice","NEEdevexpl",
                     "LEchoice","LEdevexpl")

for(i in 1:length(gppmods)){
  
  site<-substr(gppmods[i],nchar(gppmods[i])-6,nchar(gppmods[i])-4)
  print(paste(site,i))
  leafresults[i,1]<-site
  gpp<-readRDS(gppmods[i])
  choose<-choosemod3(gpp)
  
  leafresults[i,2]<-choose #2 means more complicated model
  if(choose==1){
    leafresults[i,3]<-summary(gpp[[6]])$dev.expl
  }
  if(choose==2){
    leafresults[i,3]<-summary(gpp[[14]])$dev.expl
  }
  if(choose==3){
    leafresults[i,3]<-summary(gpp[[16]])$dev.expl
  }
  
  nee<-readRDS(neemods[i])
  choose<-choosemod3(nee)
  leafresults[i,4]<-choose
  if(choose==1){
    leafresults[i,5]<-summary(nee[[6]])$dev.expl
  }
  if(choose==2){
    leafresults[i,5]<-summary(nee[[14]])$dev.expl
  }
  if(choose==3){
    leafresults[i,5]<-summary(nee[[16]])$dev.expl
  }
  
  le<-readRDS(lemods[i])
  choose<-choosemod3(le)
  leafresults[i,6]<-choose
  if(choose==1){
    leafresults[i,7]<-summary(le[[6]])$dev.expl
  }
  if(choose==2){
    leafresults[i,7]<-summary(le[[14]])$dev.expl
  }
  if(choose==3){
    leafresults[i,7]<-summary(le[[16]])$dev.expl
  }
}
write.csv(leafresults,"/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/leafresults_ALL_modsummary_with3way.csv",row.names = FALSE)


table(leafresults[,2])
table(leafresults[,4])
table(leafresults[,6])

########### Alternative: significant or not:
#for(i in 1:length(gppmods)){
#  gpp<-readRDS(gppmods[i])
#  summary(gpp[[16]])
#}

############################## #gam.check plots for chosen models - FOR SUPP
leafresults<-read.csv("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/leafresults_ALL_modsummary_with3way.csv")

setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/") #DAILY outliers removed
gppmods<-list.files(pattern="GPP_D_*")
neemods<-list.files(pattern="NEE_D_*")
lemods<-list.files(pattern="LE_D_*")

png("/Users/doris/Documents/UIowa/Figures/final/Leaf_GPPmod_eval.png",
    width=1600,height=2000,units="px",res=200)
par(mfrow=c(7,4))

for (u in 1:length(gppmods)){
  gpp<-readRDS(gppmods[u])
  r<-which(leafresults$site==substr(gppmods[u],12,14))
  if(leafresults$GPPchoice[r] ==1){
    modchoice<-gpp[[6]]
    inwords<-"3-way (Eqn.10)"
  }
  if(leafresults$GPPchoice[r] ==2){
    modchoice<-gpp[[14]]
    inwords<-"no int (Eqn.8)"
  }
  if(leafresults$GPPchoice[r] ==3){
    modchoice<-gpp[[16]]
    inwords<-"2-way (Eqn.9)"
  }
  
  res<-residuals(modchoice, type= "deviance")
  observed.y <- napredict(modchoice$na.action, modchoice$y)
  #par("plt" = c(.15,.9,.15,.9))
  par("plt" = c(.25,.95,.3,.88))
  plot(fitted(modchoice), observed.y, xlab = "", 
       ylab = "", main="",
       pch=20,col=rgb(0,0,0,.3))
  mtext(paste0(leafresults$site[r],"; ",inwords),side=3,line=0.1,cex=0.8)
  mtext(expression(paste("Fitted GPP"["s"])),side =1,line=1.8,cex=0.6)
  mtext(expression(paste("Observed GPP"["s"])),side=2,line=1.8,cex=0.6)
  
  abline(0,1,col="red")
  par(new=TRUE)
  par("plt" = c(0.25, 0.46, 0.65, 0.88))
  hist(res, xlab = "", main = "",ylab="",yaxt='n',xaxt='n')
  abline(v=0,col="red")
}
dev.off()

png("/Users/doris/Documents/UIowa/Figures/final/Leaf_NEEmod_eval.png",
    width=1600,height=2000,units="px",res=200)
par(mfrow=c(7,4))

for (u in 1:length(neemods)){
  nee<-readRDS(neemods[u])
  r<-which(leafresults$site==substr(neemods[u],12,14))
  if(leafresults$NEEchoice[r] ==1){
    modchoice<-nee[[6]]
    inwords<-"3-way (Eqn.10)"
  }
  if(leafresults$NEEchoice[r] ==2){
    modchoice<-nee[[14]]
    inwords<-"no int (Eqn.8)"
  }
  if(leafresults$NEEchoice[r] ==3){
    modchoice<-nee[[16]]
    inwords<-"2-way (Eqn.9)"
  }
  
  res<-residuals(modchoice, type= "deviance")
  observed.y <- napredict(modchoice$na.action, modchoice$y)
  par("plt" = c(.25,.95,.3,.88))
  plot(fitted(modchoice), observed.y, xlab = "", 
       ylab = "", main="",
       pch=20,col=rgb(0,0,0,.3))
  mtext(paste0(leafresults$site[r],"; ",inwords),side=3,line=0.1,cex=0.8)
  mtext(expression(paste("Fitted NEP"["s"])),side =1,line=1.8,cex=0.6)
  mtext(expression(paste("Observed NEP"["s"])),side=2,line=1.8,cex=0.6)
  
  abline(0,1,col="red")
  par(new=TRUE)
  par("plt" = c(0.25, 0.46, 0.65, 0.88))
  hist(res, xlab = "", main = "",ylab="",yaxt='n',xaxt='n')
  abline(v=0,col="red")
}
dev.off()

png("/Users/doris/Documents/UIowa/Figures/final/Leaf_LEmod_eval.png",
    width=1600,height=2000,units="px",res=200)
par(mfrow=c(7,4))

for (u in 1:length(lemods)){
  le<-readRDS(lemods[u])
  r<-which(leafresults$site==substr(lemods[u],11,13)) #hard coded numbers change bc LE has 2 letters, ugh
  if(leafresults$LEchoice[r] ==1){
    modchoice<-le[[6]]
    inwords<-"3-way (Eqn.10)"
  }
  if(leafresults$LEchoice[r] ==2){
    modchoice<-le[[14]]
    inwords<-"no int (Eqn.8)"
  }
  if(leafresults$LEchoice[r] ==3){
    modchoice<-le[[16]]
    inwords<-"2-way (Eqn.9)"
  }
  
  res<-residuals(modchoice, type= "deviance")
  observed.y <- napredict(modchoice$na.action, modchoice$y)
  par("plt" = c(.25,.95,.3,.88))
  plot(fitted(modchoice), observed.y, xlab = "", 
       ylab = "", main="",
       pch=20,col=rgb(0,0,0,.3))
  mtext(paste0(leafresults$site[r],"; ",inwords),side=3,line=0.1,cex=0.8)
  mtext(expression(paste("Fitted LE"["s"])),side =1,line=1.8,cex=0.6)
  mtext(expression(paste("Observed LE"["s"])),side=2,line=1.8,cex=0.6)
  
  abline(0,1,col="red")
  par(new=TRUE)
  par("plt" = c(0.25, 0.46, 0.65, 0.88))
  hist(res, xlab = "", main = "",ylab="",yaxt='n',xaxt='n')
  abline(v=0,col="red")
}
dev.off()

#########################
airresults<-read.csv("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/airresults_ALL_modsummary_with3way.csv")

setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/") #DAILY outliers removed
gppmods<-list.files(pattern="GPP_D_*")
neemods<-list.files(pattern="NEE_D_*")
lemods<-list.files(pattern="LE_D_*")

png("/Users/doris/Documents/UIowa/Figures/final/Air_GPPmod_eval.png",
    width=1600,height=2000,units="px",res=200)
par(mfrow=c(7,4))

for (u in 1:length(gppmods)){
  gpp<-readRDS(gppmods[u])
  r<-which(airresults$site==substr(gppmods[u],11,13))
  if(airresults$GPPchoice[r] ==1){
    modchoice<-gpp[[6]]
    inwords<-"3-way (Eqn.10)"
  }
  if(airresults$GPPchoice[r] ==2){
    modchoice<-gpp[[14]]
    inwords<-"no int (Eqn.8)"
  }
  if(airresults$GPPchoice[r] ==3){
    modchoice<-gpp[[16]]
    inwords<-"2-way (Eqn.9)"
  }
  
  res<-residuals(modchoice, type= "deviance")
  observed.y <- napredict(modchoice$na.action, modchoice$y)
  #par("plt" = c(.15,.9,.15,.9))
  par("plt" = c(.25,.95,.3,.88))
  plot(fitted(modchoice), observed.y, xlab = "", 
       ylab = "", main="",
       pch=20,col=rgb(0,0,0,.3))
  mtext(paste0(airresults$site[r],"; ",inwords),side=3,line=0.1,cex=0.8)
  mtext(expression(paste("Fitted GPP"["a"])),side =1,line=1.8,cex=0.6)
  mtext(expression(paste("Observed GPP"["a"])),side=2,line=1.8,cex=0.6)
  
  abline(0,1,col="red")
  par(new=TRUE)
  par("plt" = c(0.25, 0.46, 0.65, 0.88))
  hist(res, xlab = "", main = "",ylab="",yaxt='n',xaxt='n')
  abline(v=0,col="red")
}
dev.off()

png("/Users/doris/Documents/UIowa/Figures/final/Air_NEEmod_eval.png",
    width=1600,height=2000,units="px",res=200)
par(mfrow=c(7,4))

for (u in 1:length(neemods)){
  nee<-readRDS(neemods[u])
  r<-which(airresults$site==substr(neemods[u],11,13))
  if(airresults$NEEchoice[r] ==1){
    modchoice<-nee[[6]]
    inwords<-"3-way (Eqn.10)"
  }
  if(airresults$NEEchoice[r] ==2){
    modchoice<-nee[[14]]
    inwords<-"no int (Eqn.8)"
  }
  if(airresults$NEEchoice[r] ==3){
    modchoice<-nee[[16]]
    inwords<-"2-way (Eqn.9)"
  }
  
  res<-residuals(modchoice, type= "deviance")
  observed.y <- napredict(modchoice$na.action, modchoice$y)
  par("plt" = c(.25,.95,.3,.88))
  plot(fitted(modchoice), observed.y, xlab = "", 
       ylab = "", main="",
       pch=20,col=rgb(0,0,0,.3))
  mtext(paste0(airresults$site[r],"; ",inwords),side=3,line=0.1,cex=0.8)
  mtext(expression(paste("Fitted NEP"["a"])),side =1,line=1.8,cex=0.6)
  mtext(expression(paste("Observed NEP"["a"])),side=2,line=1.8,cex=0.6)
  
  abline(0,1,col="red")
  par(new=TRUE)
  par("plt" = c(0.25, 0.46, 0.65, 0.88))
  hist(res, xlab = "", main = "",ylab="",yaxt='n',xaxt='n')
  abline(v=0,col="red")
}
dev.off()

png("/Users/doris/Documents/UIowa/Figures/final/Air_LEmod_eval.png",
    width=1600,height=2000,units="px",res=200)
par(mfrow=c(7,4))

for (u in 1:length(lemods)){
  le<-readRDS(lemods[u])
  r<-which(airresults$site==substr(lemods[u],10,12)) #hard coded numbers change bc LE has 2 letters, ugh
  if(airresults$LEchoice[r] ==1){
    modchoice<-le[[6]]
    inwords<-"3-way (Eqn.10)"
  }
  if(airresults$LEchoice[r] ==2){
    modchoice<-le[[14]]
    inwords<-"no int (Eqn.8)"
  }
  if(airresults$LEchoice[r] ==3){
    modchoice<-le[[16]]
    inwords<-"2-way (Eqn.9)"
  }
  
  res<-residuals(modchoice, type= "deviance")
  observed.y <- napredict(modchoice$na.action, modchoice$y)
  par("plt" = c(.25,.95,.3,.88))
  plot(fitted(modchoice), observed.y, xlab = "", 
       ylab = "", main="",
       pch=20,col=rgb(0,0,0,.3))
  mtext(paste0(airresults$site[r],"; ",inwords),side=3,line=0.1,cex=0.8)
  mtext(expression(paste("Fitted LE"["a"])),side =1,line=1.8,cex=0.6)
  mtext(expression(paste("Observed LE"["a"])),side=2,line=1.8,cex=0.6)
  
  abline(0,1,col="red")
  par(new=TRUE)
  par("plt" = c(0.25, 0.46, 0.65, 0.88))
  hist(res, xlab = "", main = "",ylab="",yaxt='n',xaxt='n')
  abline(v=0,col="red")
}
dev.off()

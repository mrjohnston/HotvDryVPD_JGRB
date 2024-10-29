#Figure 6
library(mgcv)
predictFUN<-function(model,df){
  pred<-predict.gam(object=model,type="terms",se.fit=TRUE,newdata=df)
  cols<-which(colnames(pred[[1]])%in%c("VPDairmod","ti(VPDairmod)",
                                       "ti(Tadj_svp_air,VPDairmod)",
                                       "ti(VPDairmod,Tadj_svp_air)",
                                       "ti(SMAP_rootzoneSM,VPDairmod)",
                                       "ti(VPDairmod,SMAP_rootzoneSM)",
                                       "ti(VPDairmod,Tadj_svp_air,SMAP_rootzoneSM)",
                                       "ti(VPDairmod,SMAP_rootzoneSM,Tadj_svp_air)",
                                       "ti(Tadj_svp_air,SMAP_rootzoneSM,VPDairmod)",
                                       "ti(Tadj_svp_air,VPDairmod,SMAP_rootzoneSM)",
                                       "ti(SMAP_rootzoneSM,VPDairmod,Tadj_svp_air)",
                                       "ti(SMAP_rootzoneSM,Tadj_svp_air,VPDairmod)"))
  result<-rowSums(pred[[1]][,c(cols)]) #+coef(model)[1] 
  sub<-pred[[2]][,c(cols)]
  sub2<-sub^2
  resultSE<-sqrt(rowSums(sub2))#all squared)rowSums(pred[[2]][,c(cols)])
  return(list(result,resultSE))
}
choice<-read.csv("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/airresults_ALL_modsummary_with3way.csv")
#choice = 1 means 3-way interaction; 2 means no interaction; 3 means 2 way interaction 

#[1] Get useful quantities into "choice" 
info<-read.csv("/Volumes/lss_ech2oLab/AMF_CoreSites/FilteredSiteInfo_v5.csv")
infosub<-info[which(substr(info$site,4,6)%in%choice$site),]
substr(infosub$site,4,6)==choice$site #ok, same order

#Add growing season T and growing season VPD to infosub
infosub$GSVPD<-NA
infosub$GSairtemp<-NA
for (q in 1:nrow(infosub)){
  site<-substr(infosub$site,4,6)[q]
  dat<-read.csv(paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/J_DailyFiles/US_",site,"_D_SVP_Tadj.csv"))
  infosub$GSVPD[q]<-mean(dat$VPDairmod)
  infosub$GSairtemp[q]<-mean(dat$Tairmod)
}

choice$MAT<-infosub$MAT
choice$MAP<-infosub$MAP
choice$veg<-infosub$veg
choice$GSairT<-infosub$GSairtemp
choice<-choice[order(choice$GSairT),] #order by growing season airT

#Models (air)
setwd("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/") #DAILY outliers removed

gppmods<-list.files(pattern="GPP_D_*")
neemods<-list.files(pattern="NEE_D_*")
lemods<-list.files(pattern="LE_D_*")

Tlowquant<-.05
Thighquant<-.95
SMquant=0.5

reslow<-reshigh<-res0<-data.frame(matrix(NA,nrow=25,ncol=13))

varset=c("GPP","NEE","LE")
png(filename="/Users/doris/Documents/UIowa/Figures/ThreeDensities_AIR.png",
    width=2400,height=3215,units="px",res=220)
# zones <- matrix(c(1,15,28,
#                   2,16,29,
#                   3,17,30,
#                   4,18,31,
#                   5,19,32,
#                   6,20,33,
#                   7,21,34,
#                   8,22,35,
#                   9,23,36,
#                   10,24,0,
#                   11,25,0,
#                   12,26,0,
#                   13,27,0,
#                   14,0,0),
zones <- matrix(c(1,16,29,
                  2,17,30,
                  3,18,31,
                  4,19,32,
                  5,20,33,
                  6,21,34,
                  7,22,35,
                  8,23,36,
                  9,24,37,
                  10,25,0,
                  11,26,0,
                  12,27,0,
                  13,28,0,
                  14,0,0,
                  15,0,0),
                ncol = 3, byrow = TRUE)
layout(zones, widths=c(1,1,1),
       heights = c(1,1,1,1,1,1,1,1,1,1,1,1,1))
par(mar=c(2,2,1,.5)) #bottom, left, top, right


for(u in 1:length(varset)){
  var=varset[u]
  
  for(i in 1:nrow(choice)){
    site<-choice$site[i]
    
    reslow[i,1]<-site
    reshigh[i,1]<-site
    res0[i,1]<-site
    
    if(var=="GPP"){
      gpp<-readRDS(gppmods[grep(site,gppmods)])
      if(choice$GPPchoice[i]==1){
        mod<-gpp[[6]]
      }
      if(choice$GPPchoice[i]==2){
        mod<-gpp[[14]]
        #next
      }
      if(choice$GPPchoice[i]==3){
        mod<-gpp[[16]]
      }
    }
    
    if(var=="NEE"){
      nee<-readRDS(neemods[grep(site,neemods)])
      if(choice$NEEchoice[i]==1){
        mod<-nee[[6]]
      }
      if(choice$NEEchoice[i]==2){
        mod<-nee[[14]]
        #next
      }
      if(choice$NEEchoice[i]==3){
        mod<-nee[[16]]
      }
    }
    
    if(var=="LE"){
      le<-readRDS(lemods[grep(site,lemods)])
      if(choice$LEchoice[i]==1){
        mod<-le[[6]]
      }
      if(choice$LEchoice[i]==2){
        mod<-le[[14]]
        # next
      }
      if(choice$LEchoice[i]==3){
        mod<-le[[16]]
      }
    }
    
    #ok, got my site & my model.. Now I need the dataframes.
    dat<-read.csv(paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/J_DailyFiles/US_",site,"_D_SVP_Tadj.csv"))
    dat$dat.Year<-as.factor(dat$dat.Year) #UGH, why didn't this work when I tried to save as part of J0??
    
    preddf<-dat[,which(names(dat)%in%c("SMAP_rootzoneSM","dat.Year","VPDairmod","Tadj_svp_air"))]
    preddfT0<-preddfTlow<-preddfThigh<-preddf
    
    preddfT0$Tadj_svp_air<-0
    preddfT0$SMAP_rootzoneSM<-quantile(preddf$SMAP_rootzoneSM,SMquant)
    
    preddfTlow$Tadj_svp_air<-quantile(preddf$Tadj_svp_air,Tlowquant)
    preddfTlow$SMAP_rootzoneSM<-quantile(preddf$SMAP_rootzoneSM,SMquant)
    
    preddfThigh$Tadj_svp_air<-quantile(preddf$Tadj_svp_air,Thighquant)
    preddfThigh$SMAP_rootzoneSM<-quantile(preddf$SMAP_rootzoneSM,SMquant)
    
    #USE MODEL TO PREDICT - only VPD terms (no intercept)
    
    predT0y<-predictFUN(mod,preddfT0)[[1]] 
    predTlow<-predictFUN(mod,preddfTlow)[[1]]
    predThigh<-predictFUN(mod,preddfThigh)[[1]]
    
    denlow<-density(predTlow)
    denhigh<-density(predThigh)
    den0<-density(predT0y)
    
    #PLOT THE RESULTS - only if there's an interaction
    if(max(abs(predTlow/predT0y)!=1)){
      
      if(var=="GPP"){
        reslow[i,2]<-quantile(predTlow,0.25)
        reslow[i,3]<-quantile(predTlow,0.5)
        reslow[i,4]<-quantile(predTlow,0.75)
        reslow[i,5]<-median(denlow$x[which(denlow$y==max(denlow$y))])#just in case multiple maxes
        reshigh[i,2]<-quantile(predThigh,0.25)
        reshigh[i,3]<-quantile(predThigh,0.5)
        reshigh[i,4]<-quantile(predThigh,0.75)
        reshigh[i,5]<-median(denhigh$x[which(denhigh$y==max(denhigh$y))])
        res0[i,2]<-quantile(predT0y,0.25)
        res0[i,3]<-quantile(predT0y,0.5)
        res0[i,4]<-quantile(predT0y,0.75)
        res0[i,5]<-median(den0$x[which(den0$y==max(den0$y))])#just in case multiple maxes
      }
      if(var=="NEE"){
        reslow[i,6]<-quantile(predTlow,0.25)
        reslow[i,7]<-quantile(predTlow,0.5)
        reslow[i,8]<-quantile(predTlow,0.75)
        reslow[i,9]<-median(denlow$x[which(denlow$y==max(denlow$y))])#just in case multiple maxes
        reshigh[i,6]<-quantile(predThigh,0.25)
        reshigh[i,7]<-quantile(predThigh,0.5)
        reshigh[i,8]<-quantile(predThigh,0.75)
        reshigh[i,9]<-median(denhigh$x[which(denhigh$y==max(denhigh$y))])
        res0[i,6]<-quantile(predT0y,0.25)
        res0[i,7]<-quantile(predT0y,0.5)
        res0[i,8]<-quantile(predT0y,0.75)
        res0[i,9]<-median(den0$x[which(den0$y==max(den0$y))])#just in case multiple maxes
      }
      if(var=="LE"){
        reslow[i,10]<-quantile(predTlow,0.25)
        reslow[i,11]<-quantile(predTlow,0.5)
        reslow[i,12]<-quantile(predTlow,0.75)
        reslow[i,13]<-median(denlow$x[which(denlow$y==max(denlow$y))])#just in case multiple maxes
        reshigh[i,10]<-quantile(predThigh,0.25)
        reshigh[i,11]<-quantile(predThigh,0.5)
        reshigh[i,12]<-quantile(predThigh,0.75)
        reshigh[i,13]<-median(denhigh$x[which(denhigh$y==max(denhigh$y))])
        res0[i,10]<-quantile(predT0y,0.25)
        res0[i,11]<-quantile(predT0y,0.5)
        res0[i,12]<-quantile(predT0y,0.75)
        res0[i,13]<-median(den0$x[which(den0$y==max(den0$y))])#just in case multiple maxes
      }
      
      quantilecutoff<-0.1
      plot(denlow,xlim=c(min(denlow$x[which(denlow$y>quantile(denlow$y,quantilecutoff))],
                             denhigh$x[which(denhigh$y>quantile(denhigh$y,quantilecutoff))],
                             den0$x[which(den0$y>quantile(den0$y,quantilecutoff))]),
                         max(denlow$x[which(denlow$y>quantile(denlow$y,quantilecutoff))],
                             denhigh$x[which(denhigh$y>quantile(denhigh$y,quantilecutoff))],
                             den0$x[which(den0$y>quantile(den0$y,quantilecutoff))])),col="blue",main="",
           ylim=c(0,max(denlow$y,denhigh$y,den0$y)),yaxt='n',cex.axis=1.5)
      polygon(denlow,col=rgb(0,0,1,0.2),border=rgb(0,0,1,0.2))
      lines(denhigh,col="red")
      polygon(denhigh,col=rgb(1,0,0,0.2),border=rgb(1,0,0,0.2))
      lines(den0,col="black")
      polygon(den0,col=rgb(0,0,0,0.2),border=rgb(1,0,0,0.2))
      
      if(var%in%c("GPP","LE")){
        mtext(paste0(site,"; ", var,"a"),side=3,line=-1.2,adj=0) 
      }
      if(var=="NEE"){
        mtext(paste0(site,"; NEPa"),side=3,line=-1.2,adj=0) 
      }
    }
  }
} #Close 'u' loop; on to the next variable

dev.off()

options(scipen=999)
names(reshigh)<-names(reslow)<-
  names(res0)<-c("site","GPPQ1","GPPmed","GPPQ3","GPPDenMaxX",
                                 "NEPQ1","NEPmed","NEPQ3","NEPDenMaxX",
                                 "LEQ1","LEmed","LEQ3","LEDenMaxX")
#reslow[reslow == 0] <- NA
#reshigh[reshigh == 0] <- NA

# ###############
# abs(reshigh$GPPQ3-reshigh$GPPQ1)
# abs(reslow$GPPQ3-reslow$GPPQ1)
# abs(res0$GPPQ3-res0$GPPQ1)
# 
# #IQRs
# abs(reshigh$GPPQ3-reshigh$GPPQ1) - abs(res0$GPPQ3-res0$GPPQ1) #4/14 neg
# abs(reslow$GPPQ3-reslow$GPPQ1) - abs(res0$GPPQ3-res0$GPPQ1) #4/14 neg
# 
# abs(reshigh$NEPQ3-reshigh$NEPQ1) - abs(res0$NEPQ3-res0$NEPQ1) #2/13 neg
# abs(reslow$NEPQ3-reslow$NEPQ1) - abs(res0$NEPQ3-res0$NEPQ1) #7/13 neg
# 
# abs(reshigh$LEQ3-reshigh$LEQ1) - abs(res0$LEQ3-res0$LEQ1) #5/9 neg
# abs(reslow$LEQ3-reslow$LEQ1) - abs(res0$LEQ3-res0$LEQ1) #3/9 neg
# 
# 
# summary(reshigh$GPPDenMaxX)
# summary(reslow$GPPDenMaxX)
# summary(res0$GPPDenMaxX)
# 
# abs(reshigh$GPPDenMaxX) - abs(res0$GPPDenMaxX) #2/14 are negative
# 
# abs(reslow$GPPDenMaxX) - abs(res0$GPPDenMaxX) #4/13 are negative
# abs(reshigh$GPPDenMaxX) - abs(reslow$GPPDenMaxX) #5/14 are negative
# 
# abs(reshigh$NEPDenMaxX) - abs(res0$NEPDenMaxX) #2/13 negative
# abs(reslow$NEPDenMaxX) - abs(res0$NEPDenMaxX) #9/13 negative
# abs(reshigh$NEPDenMaxX) - abs(reslow$NEPDenMaxX) #3/13 negative
# 
# abs(reshigh$LEDenMaxX) - abs(res0$LEDenMaxX) #5/9 negative
# abs(reslow$LEDenMaxX) - abs(res0$LEDenMaxX) #1/9 negative
# abs(reshigh$LEDenMaxX) - abs(reslow$LEDenMaxX) #7/9 negative
# ###
# reshigh$GPPDenMaxX - res0$GPPDenMaxX #3/14 are negative
# summary(reshigh$GPPDenMaxX - res0$GPPDenMaxX) 
# 
# reshigh$NEPDenMaxX - res0$NEPDenMaxX #1/13 negative
# summary(reshigh$NEPDenMaxX - res0$NEPDenMaxX) 

###PAPER NUMBERS:
#Most important for air:
reshigh$GPPDenMaxX - res0$GPPDenMaxX
summary(reshigh$GPPDenMaxX - res0$GPPDenMaxX )
summary((reshigh$GPPDenMaxX - res0$GPPDenMaxX)/abs(res0$GPPDenMaxX))

reshigh$NEPDenMaxX - res0$NEPDenMaxX
summary(reshigh$NEPDenMaxX - res0$NEPDenMaxX)
summary((reshigh$NEPDenMaxX - res0$NEPDenMaxX)/abs(res0$NEPDenMaxX))

reslow$GPPDenMaxX - res0$GPPDenMaxX
reslow$NEPDenMaxX - res0$NEPDenMaxX

reslow$LEDenMaxX - res0$LEDenMaxX
summary(reslow$LEDenMaxX - res0$LEDenMaxX)
summary((reslow$LEDenMaxX - res0$LEDenMaxX)/abs(res0$LEDenMaxX))

#redefine to reasily compare to res for leaves
reshigha<-reshigh
reslowa<-reslow
res0a<-res0
#write.csv(reslow,"/Users/doris/Documents/UIowa/reslow_leaf.csv")
#write.csv(reshigh,"/Users/doris/Documents/UIowa/reshigh_leaf.csv")

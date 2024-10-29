#Final (?) version of H1a_makeresFILTER.R + H1a_makeSiteTogetherplots.R

library(viridis)
library(plyr)
setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/H2.fin_Filtered_withRH/")
files<-list.files()

interval=0.1

#Remove extreme VPDs, temperatures, and RHs...
RemoveExtremes<-function(datsub){
  VPDairmin<-quantile(datsub$VPDairmod,0.05,na.rm = TRUE)
  VPDairmax<-quantile(datsub$VPDairmod,0.95,na.rm = TRUE)
  VPDleafmin<-quantile(datsub$VPDsurfEeqnmod,0.05,na.rm = TRUE)
  VPDleafmax<-quantile(datsub$VPDsurfEeqnmod,0.95,na.rm = TRUE)
  
  dq05rH<-quantile(datsub$rH,.05,na.rm = TRUE)
  dq95rH<-quantile(datsub$rH,.95,na.rm = TRUE)
  
  dq05Tair<-quantile(datsub$Tairmod,.05,na.rm = TRUE)
  dq95Tair<-quantile(datsub$Tairmod,.95,na.rm = TRUE)
  dq05Tleaf<-quantile(datsub$TsurfEeqnmod,.05,na.rm = TRUE)
  dq95Tleaf<-quantile(datsub$TsurfEeqnmod,.95,na.rm = TRUE)
  
  datsubair<-datsub[which(datsub$VPDairmod>=VPDairmin & datsub$VPDairmod<=VPDairmax &
                            datsub$rH>=dq05rH & datsub$rH<=dq95rH &
                            datsub$Tairmod>=dq05Tair & datsub$Tairmod<=dq95Tair),]
  datsubleaf<-datsub[which(datsub$VPDsurfEeqnmod>=VPDleafmin & datsub$VPDsurfEeqnmod<=VPDleafmax &
                            datsub$rH>=dq05rH & datsub$rH<=dq95rH &
                            datsub$TsurfEeqnmod>=dq05Tleaf & datsub$TsurfEeqnmod<=dq95Tleaf),]
  return(list(datsubair,datsubleaf))
}

for(i in 1:length(files)){
  print(i)
  dat<-read.csv(files[i])
  
  datsub<-dat[which(dat$growseasfilter=="keep" & dat$growtimefilter=="keep"),]
  dair<-RemoveExtremes(datsub)[[1]]
  if(i!=2){
      dleaf<-RemoveExtremes(datsub)[[2]] #doesn't actually throw an error for Ha1, just has 0 rows
    }
  #Air
    rangeTa<-max(dair$Tairmod)-min(dair$Tairmod) #total range
    rangerHa<-max(dair$rH)-min(dair$rH) #total range
    VPDmina<-min(dair$VPDairmod)
    VPDmaxa<-max(dair$VPDairmod)
    VPDmediansitea<-median(dair$VPDairmod) #median VPD at site
 #Leaf
    if(i!=2){
    rangeTl<-max(dleaf$TsurfEeqnmod)-min(dleaf$TsurfEeqnmod) #total range
    rangerHl<-max(dleaf$rH)-min(dleaf$rH) #total range
    VPDminl<-min(dleaf$VPDsurfEeqnmod)
    VPDmaxl<-max(dleaf$VPDsurfEeqnmod)
    VPDmediansitel<-median(dleaf$VPDsurfEeqnmod) #median VPD at site
}
  #Set up the sequence of VPDs to subset (changes by site)
  starta<-round_any(VPDmina, .1, f = floor) #function from plyr
  enda<-round_any(VPDmaxa, .1, f = ceiling)
  s1a<-s2a<-seq(starta,enda,by=interval) 
  s1a<-s1a[-length(s1a)]
  s2a<-s2a[-1]
  
  if(i!=2){
  startl<-round_any(VPDminl, .1, f = floor) #function from plyr
  endl<-round_any(VPDmaxl, .1, f = ceiling)
  s1l<-s2l<-seq(startl,endl,by=interval) 
  s1l<-s1l[-length(s1l)]
  s2l<-s2l[-1]
  } 
  airres<-leafres<-data.frame(matrix(NA,ncol=14))

  for(j in 1:length(s1a)){
    
      suba<-dair[which(dair$VPDairmod<s2a[j] & dair$VPDairmod>s1a[j]),]
    
      sq05Ta<-quantile(suba$Tairmod,.05,na.rm = TRUE)
      sq50Ta<-quantile(suba$Tairmod,.5,na.rm = TRUE)
      sq95Ta<-quantile(suba$Tairmod,.95,na.rm = TRUE)
      sq05rHa<-quantile(suba$rH,.05,na.rm = TRUE)
      sq50rHa<-quantile(suba$rH,.5,na.rm = TRUE)
      sq95rHa<-quantile(suba$rH,.95,na.rm = TRUE)
      
      srangeT_90a<-sq95Ta - sq05Ta
      srangerH_90a<- sq95rHa - sq05rHa
      
      #Middle 90% of data
      normrangeT_90a<- srangeT_90a / rangeTa
      normrangerH_90a<- srangerH_90a / rangerHa
      airres[j,1]<-substr(files[i],1,6)
      airres[j,2]<-median(dair$VPDairmod)
      airres[j,3]<-(s1a[j]+s2a[j])/2
      airres[j,4]<-srangeT_90a
      airres[j,5]<-srangerH_90a
      airres[j,6]<-normrangeT_90a
      airres[j,7]<-normrangerH_90a
      airres[j,8]<-nrow(suba)
      airres[j,c(9:14)]<-c(sq05Ta,sq50Ta,sq95Ta,sq05rHa,sq50rHa,sq95rHa)
  } #end j1

  if(i!=2){
for(j in 1:length(s1l)){
      subl<-dleaf[which(dleaf$VPDsurfEeqnmod<s2l[j] & dleaf$VPDsurfEeqnmod>s1l[j]),]
      
      sq05Tl<-quantile(subl$TsurfEeqnmod,.05,na.rm = TRUE)
      sq50Tl<-quantile(subl$TsurfEeqnmod,.5,na.rm = TRUE)
      sq95Tl<-quantile(subl$TsurfEeqnmod,.95,na.rm = TRUE)
      sq05rHl<-quantile(subl$rH,.05,na.rm = TRUE)
      sq50rHl<-quantile(subl$rH,.5,na.rm = TRUE)
      sq95rHl<-quantile(subl$rH,.95,na.rm = TRUE)

      srangeT_90l<-sq95Tl - sq05Tl
      srangerH_90l<- sq95rHl - sq05rHl

    #Middle 90% of data
    normrangeT_90l<- srangeT_90l / rangeTl
    normrangerH_90l<- srangerH_90l / rangerHl
    
    leafres[j,1]<-substr(files[i],1,6)
    leafres[j,2]<-median(dleaf$VPDsurfEeqnmod)
    leafres[j,3]<-(s1l[j]+s2l[j])/2
    leafres[j,4]<-srangeT_90l
    leafres[j,5]<-srangerH_90l
    leafres[j,6]<-normrangeT_90l
    leafres[j,7]<-normrangerH_90l
    leafres[j,8]<-nrow(subl)
    leafres[j,c(9:14)]<-c(sq05Tl,sq50Tl,sq95Tl,sq05rHl,sq50rHl,sq95rHl)
}} #end j2
  
  write.csv(airres,paste0("/Users/doris/git/HotvDryVPD/H1aCSVs/",substr(files[i],1,6),"_AIR_HH.csv"),row.names = FALSE)
  
  if(i!=2){
  write.csv(leafres,paste0("/Users/doris/git/HotvDryVPD/H1aCSVs/",substr(files[i],1,6),"_LEAF_HH.csv"),row.names = FALSE)
  }
   } #end site loop (i)
  


  
  
  
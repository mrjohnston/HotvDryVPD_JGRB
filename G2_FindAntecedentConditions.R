library(lutz)

#Get antecedent conditions, pre-filter.
setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/G2.3_CombineSMAPandRPOutandPrecip/") 
files<-list.files()
`%nin%` = Negate(`%in%`)

for (i in 1:length(files)){

  print(files[i])
  dat<-read.csv(files[i])
  site<-substr(files[i],1,6)
  
  #Ensure that the data are continuous
  if(site%nin%c("US_Ha1","US_MMS")){ #HALF-HOURLY
  timeseq<-seq(from=as.POSIXct(dat$posix_utc,format="%Y-%m-%d %H:%M:%S",tz="UTC")[1],
                to=as.POSIXct(dat$posix_utc,format="%Y-%m-%d %H:%M:%S",tz="UTC")[nrow(dat)],
                by="30 min")
  } #UTc tag is v important, otherwise things get mucky muck w/ daylight savings!
  
  if(site%in%c("US_Ha1","US_MMS")){ #HOURLY
    timeseq<-seq(from=as.POSIXct(dat$posix_utc,format="%Y-%m-%d %H:%M:%S",tz="UTC")[1],
               to=as.POSIXct(dat$posix_utc,format="%Y-%m-%d %H:%M:%S",tz="UTC")[nrow(dat)],
               by="1 hour")
  }
  print(nrow(dat)==length(timeseq)) #want this to be TRUE
  
  #Ok, data are continuous. Now add the antecedent conditions required for filter
  #Make precip and swc filter columns
  
  #SWC differences - anytime in the previous HOUR (was: previous timestep, but this standardizes)
  dat$swc_diff_lag1<-NA
  dat$swc_diff_lag1<-c(NA,diff(dat[,grep("SWC",names(dat))],lag=1,differences = 1))
  dat$swc_diff_lag2<-NA
  dat$swc_diff_lag2<-c(NA,NA,diff(dat[,grep("SWC",names(dat))],lag=2,differences = 1))
  
  if(site%in%c("US_Ha1","US_MMS")){ #HOURLY - remove lag2 b/c irrelevant (it would mean 2 hours)
    dat$swc_diff_lag2<-NA
  }
  #So the filter should be: if either of these things is >0.02, omit.

  #Precipitation sums - over the last 24 & 48 hours - this code is really inefficient.
  #Keep the sums NA if there's no data
  dat$precip_sums24hr<-NA
  dat$precip_sums48hr<-NA
  
  for(q in 1:nrow(dat)){
    
    if(dat$Hour[2]-dat$Hour[1]==1){ #HOURLY
      
      if(q>23){
        if(sum(!is.na(dat$P[c(q-23):q]))>0){ #if all NA, sum stays NA
          dat$precip_sums24hr[q]<-sum(dat$P[c(q-23):q],na.rm=TRUE)
        } }
      
      if(q>47){
        if(sum(!is.na(dat$P[c(q-47):q]))>0){
          dat$precip_sums48hr[q]<-sum(dat$P[c(q-47):q],na.rm=TRUE)
        } }
    }
    
    if(dat$Hour[2]-dat$Hour[1]==0.5){ #HALF-HOURLY
      
      if(q>47){
        if(sum(!is.na(dat$P[c(q-47):q]))>0){
          dat$precip_sums24hr[q]<-sum(dat$P[c(q-47):q],na.rm=TRUE)
        } }
      
      if(q>95){
        if(sum(!is.na(dat$P[c(q-95):q]))>0){
          dat$precip_sums48hr[q]<-sum(dat$P[c(q-95):q],na.rm=TRUE)
        } }
      
    }
  }
  #################

  dat$WS_ant<-c(NA,dat$WS[-nrow(dat)])
  
  write.csv(dat,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/G2.4_AntecedentConditions/",site,"_withAnte.csv"))
}

#UGH, again. This script really SHOULD include the soil moisture antecedent.
#But it's in the next script (H_Filter).




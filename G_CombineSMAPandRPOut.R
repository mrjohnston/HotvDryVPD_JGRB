# [1] Combine SMAP extractions with output of REddyProc
# Updated May 2023, no longer calculating WUE, Gs.
# Compare to old script: Shared/lss_ech2oLab/miriam/scripts/CombineRPandSMAP.R

# [2] Add precipitation from the raw file -- I really should have extracted this prior to REddyProc :-(

#Combine SMAP & REddyProc data
library(lubridate)

#Only have to do this once...
#info<-read.csv("/Volumes/lss_ech2oLab/miriam/scripts/info_core.csv")
#info$tz <- lutz::tz_lookup_coords(info$LOCATION_LAT, info$LOCATION_LONG, method="accurate")
#info$utcoffset<-NA
#for(y in 1:nrow(info)){
#  info$utcoffset[y]<-lutz::tz_offset("2018-01-01", info$tz[y])$utc_offset_h
#}
#write.csv(info, "/Volumes/lss_ech2oLab/miriam/scripts/info_core_withTZ.csv")

SMAPfiles<-list.files(path = "/Volumes/lss_ech2oLab/AMF_CoreSites/E2_SMAPOut/ExtractedData/", full.names = TRUE)
RPfiles<-list.files(path = "/Volumes/lss_ech2oLab/AMF_CoreSites/C2_ReddyProcOut/",full.names = TRUE)
info<-read.csv("/Volumes/lss_ech2oLab/miriam/scripts/info_core_withTZ.csv")

for(i in 1:length(RPfiles)){ #NEXT TIME: START WITH MMS
  gc()
  site<-substr(RPfiles[i],nchar(RPfiles[i])-15,nchar(RPfiles[i])-13)
  print(paste0(site, ", ",i,"/",length(RPfiles)))
  RP<-read.csv(RPfiles[i])
  SMAP<-read.csv(SMAPfiles[which(substr(unlist(lapply(strsplit(SMAPfiles,"//"),"[[",2)),4,6)==site)])
  siteinfo<-info[which(substr(info$SITE_ID,4,6)==site),]
  utcoffset <- siteinfo$utcoffset
  
  #Note that the RP times are local standard time, NOT UTC
  #Get UTC (posix format) time in RP file
  RP$min<-0
  RP$min[grep(".5",RP$Hour,fixed=TRUE)]<-30
  RP$Hour_full<-RP$Hour
  RP$Hour_full[grep(".5",RP$Hour,fixed=TRUE)]<-RP$Hour_full[grep(".5",RP$Hour,fixed=TRUE)]-0.5
  RP$posix_local<-as.POSIXct(paste(RP$Year,RP$DoY,RP$Hour_full,RP$min),format = "%Y %j %H %M",tz="UTC")
  #^TZ is always the standard version, not daylight version -- so call it UTC, so there's no daylight savings,
  #and then convert the time accordingly (below)
  RP$posix_utc<-RP$posix_local-(utcoffset*60*60)
  
  #Get UTC (posix format) time in SMAP
  SMAP$posix_utc<-as.POSIXct(paste(SMAP$year,SMAP$month,SMAP$day,SMAP$time_intervalcenter_utc),
                             format = "%Y %m %d %H:%M",tz="UTC")
  
  #Match them!
  RP$SMAP_rootzoneSM <- RP$SMAP_surfaceSM <- RP$SMAP_posix_utc <- NA
  nrowRP<-nrow(RP)
  
  for(j in 1:nrow(RP)){ #Someday, I should use an apply function?
    #for(j in 1:10000){
    
    if(j %% 10000==0) {
      cat(paste0(j, "/" , nrowRP, " --- "))
    }
    
    rowselec<-which(abs(RP$posix_utc[j]-SMAP$posix_utc) == min(abs(RP$posix_utc[j]-SMAP$posix_utc)))
    rowselec<-rowselec[1] #Because some times will be equidistant to 2 SMAP acquisitions.
    
    selec<-SMAP[rowselec, c(which(names(SMAP)=="rootzone_SM"),
                            which(names(SMAP)=="surface_SM"),
                            which(names(SMAP)=="posix_utc"))]
    RP[j,c(which(names(RP)=="SMAP_rootzoneSM"),
           which(names(RP)=="SMAP_surfaceSM"),
           which(names(RP)=="SMAP_posix_utc"))]<-selec
    
  } #end loop through rows of file
  
  RP$SMAP_posix_utc<-as.POSIXct(RP$SMAP_posix_utc,origin="1970-01-01", tz="UTC")
  write.csv(RP,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/G2.2_CombineSMAPandRPOut/US_",site,"_merge.csv"),row.names=F)
  
} #end loop through files

##############################################################################
#[2] Ugh, I actually need precip from the raw files -- add here.
source("/Users/doris/git/HotvDryVPD/date2doy_fxn.R")
library(amerifluxr)

rawfiles<-list.files(path = "/Volumes/lss_ech2oLab/AMF_CoreSites/A2_Raw/", full.names = TRUE) 
Combinefiles<-list.files(path = "/Volumes/lss_ech2oLab/AMF_CoreSites/G2.2_CombineSMAPandRPOut/",full.names = TRUE) 

for(i in 24:length(Combinefiles)){

  #Get the RP and SMAP Combined and AMF raw files:
  site<-substr(Combinefiles[i],nchar(Combinefiles[i])-12,nchar(Combinefiles[i])-10)
  print(paste("REddyProc site ",i,": ",site))
  R<-read.csv(Combinefiles[i])
  print(paste("Reading raw file: " , rawfiles[which(substr(unlist(lapply(strsplit(rawfiles,"//"),"[[",2)),8,10)==site)]))
  raw<- amf_read_base(file = rawfiles[which(substr(unlist(lapply(strsplit(rawfiles,"//"),"[[",2)),8,10)==site)],
                      unzip = TRUE,
                      parse_timestamp = FALSE)
  
  #Get Year, DoY, Hour columns on the raw file to match with RP
  R$tomatch<-paste0(R$Year,".",R$DoY,".",R$Hour)
  raw.time <- as.character(as.matrix(raw$TIMESTAMP_START)) #This is local standard time.
  yr <- as.numeric(substr(raw.time, 1, 4))
  mo <- as.numeric(substr(raw.time, 5, 6))
  dy <- as.numeric(substr(raw.time, 7, 8))
  hr <- as.numeric(substr(raw.time, 9, 10))
  mn <- as.numeric(substr(raw.time, 11, 12))
  doy <- as.numeric(mapply(date2doy, yr=yr, mo=mo, dy=dy))
  raw$Year <- yr
  raw$DoY <- doy
  raw$Hour <- hr + mn/60 + 0.5  #Ending hour of the half-hourly period
  rm(list = c('yr', 'mo', 'dy', 'hr', 'mn', 'doy', 'raw.time'))
  raw$tomatch<-paste0(raw$Year,".",raw$DoY,".",raw$Hour)
  head<-names(raw)
  
  #If it exists, get P for the filter
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  if(site=="BRG"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P_1_1_1', head))] 
  }
  if(site=="Ha1"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P$', head))] 
  }
    if(site=="Ho1"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P_RAIN_1_1_1', head))] 
  }
  if(site=="Ho2"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P_RAIN_1_1_1', head) )] 
    }
  if(site=="Me6"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head) )] 
  }
  if(site=="MMS"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P_1_1_1$', head))] 
  }
  if(site=="NC2"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),grep(pattern ='^P_1_1_1$', head))] 
  }
  if(site=="NC3"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P_1_1_1$', head))] 
  }
  if(site=="NR1"){  
    raw2 <- raw[, c(grep(pattern ='tomatch', head),grep(pattern ='^P_1_1_1$', head))] 
    #duplicated R files here...troubleshooting:
    ##not duplicated in Ameriflux raw
    #testB<-read.csv("/Volumes/lss_ech2oLab/AMF_CoreSites/B2_ReddyProcIn/US_NR1_REddyIN.csv")
    #testB$dd<-paste0(testB$Year,testB$DoY,testB$Hour)
    #sum(duplicated(testB$dd)) #22591 - weird... 
    ##ok, it's because Year+DoY+Hour XXXX 1 10.5 = XXXX 11 0.5 -- but why wasn't this an issue for others?
    #Fixed now for all, re-run from the beginning.
  }
  if(site=="Ro4"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),grep(pattern ='^P$', head))] 
  }
  if(site=="Seg"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),  grep(pattern ='^P$', head))] 
  }
  if(site=="Ses"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P$', head))] 
     }
  if(site=="SRG"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),grep(pattern ='P$', head) )] 
  }
  if(site=="SRM"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='P$', head))] 
  }
  if(site=="Syv"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P_PI_F$', head))] 
  }
  if(site=="Ton"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),grep(pattern ='^P$', head))] 
  }
  if(site=="UMB"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head))] 
  }
  if(site=="UMd"){ #NO PRECIP HERE
    raw2 <- raw[, c(grep(pattern ='tomatch', head))] 
    P<-rep(NA,length(raw2))
    raw2<-cbind(raw2,P)
    raw2<-as.data.frame(raw2)
  }
  if(site=="Var"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head) )] 
  }
  if(site=="Vcm"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head))] 
  }
  if(site=="Vcp"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),   grep(pattern ='^P$', head))] 
      }
  if(site=="Vcs"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head),grep(pattern ='^P$', head))] 
  }
  if(site=="WCr"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P_PI_F$', head))] 
  }
  if(site=="Whs"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head))] 
  }
  if(site=="Wjs"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head))] 
      }
  if(site=="Wkg"){ 
    raw2 <- raw[, c(grep(pattern ='tomatch', head), grep(pattern ='^P$', head))] 
  }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  names(raw2)<-c("tomatch2","P")
  raw2$P[which(raw2$P<0)]<-NA
  
  #Put the ReddyProc file and the raw file together; double check the match
  raw3<-raw2[which(raw2$tomatch2%in%R$tomatch),]
  tog<-cbind(R,raw3)
  if(!all.equal(tog$tomatch,tog$tomatch2)){ #I don't expect this ever to happen
    stop("Dates not matching")
  }
  tog<-tog[,-which(names(tog)%in%c("tomatch","tomatch2"))]
  
  write.csv(tog,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/G2.3_CombineSMAPandRPOutandPrecip/US_",site,"_merge2.csv"),row.names=F)
}


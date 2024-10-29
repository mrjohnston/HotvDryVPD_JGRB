#REddyProc - THIS SHOULD BE RUN ON THE CLUSTER; IT IS WAY TOO SLOW.
    #To Run: qsub REddyProcX.job, where X ={1,2,3,4,5} in: /Shared/lss_ech2oLab/miriam/scripts
    #These call B_REddyProcX.R, where X={1,2,3,4,5} in: /Shared/lss_ech2oLab/miriam/scripts
    #Each of the R scripts loops through 5 or 6 of the 26 files
    #Note: issues with amerifluxr on the Cluster; B_REddyProcX.R reads in "info." 
library("amerifluxr")
library(REddyProc)
library(dplyr)
library(lutz)
library(bigleaf)
library(suncalc)

setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/B2_ReddyProcIn/")
filenames<-list.files("/Volumes/lss_ech2oLab/AMF_CoreSites/B2_ReddyProcIn/",full.names = FALSE)
sites<-substr(filenames,4,6)
all<-amf_site_info()
info<-all[which(sapply(strsplit(all$SITE_ID,"-"),"[[",2)%in%sites),
         which(colnames(all)%in%c("SITE_ID","LOCATION_LAT","LOCATION_LONG"))]

#####################################
#FOR CLUSTER:
#library(REddyProc)
#library(dplyr)
#library(lutz)
#library(bigleaf)
#setwd("/Shared/lss_ech2oLab/AMF_CoreSites/B2_ReddyProcIn/")
#filenames<-list.files("/Shared/lss_ech2oLab/AMF_CoreSites/B2_ReddyProcIn/",full.names = FALSE)
#sites<-substr(filenames,4,6)
#info<-read.csv("/Shared/lss_ech2oLab/miriam/scripts/info_core.csv") 
#####################################

for(i in 1:length(sites)){
#for(i in 15:15){
  print(paste("SITE IS ",sites[i],"########################"))
  idx <- grep(sites[i], info$SITE_ID, ignore.case=TRUE)
  lat <- info$LOCATION_LAT[idx]
  lon <- info$LOCATION_LON[idx]
  tz <- tz_lookup_coords(lat, lon, method = "fast", warn = FALSE)
  utc <- tz_offset("2018-01-01", tz)$utc_offset_h
  
  ec.s<-read.csv(filenames[grep(sites[i],filenames)])
  
  #Convert logical to numeric
  ids<-which(sapply(ec.s,mode)=="logical")
  for(j in ids){
    ec.s[,j]<-as.numeric(ec.s[,j])
  }
  
  #ec.s<-ec.s[,which(colnames(ec.s)%in%c("Year","DoY","Hour","NEE","LE","Rg","SWout","LWout","Tair","rH","Ustar","VPD"))]
  head <- names(ec.s)
  
  #Deal with site-specific issues:
  if(sites[i]=="BRG"){ #no NEE or LE in 2016
    ec.s<-ec.s[-which(ec.s$Year==2016),]
  }
  if(sites[i]=="Ha1" | sites[i]=="MMS"){ #hourly data are on 0.5; add 0.5 to all.
    ec.s$Hour<-ec.s$Hour+0.5
  }
  if(sites[i]=="Ho1"){ #no ustar in 1995
    ec.s<-ec.s[-which(ec.s$Year==1995),]
  }
  if(sites[i]=="Me6"){ #no Tair 2010-2018
    ec.s<-ec.s[-which(ec.s$Year%in%c(2010,2011,2012,2013,2014,2015,2016,2017,2018)),]
  }
  if(sites[i]=="SRM"){ #no Ustar, NEE, LE, Rg for 2003
    ec.s<-ec.s[-which(ec.s$Year==2003),]
  }
  if(sites[i]=="Syv"){ #Missing Rg
    ec.s<-ec.s[-which(ec.s$Year%in%c(2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011)),]
  }
  if(sites[i]=="UMd"){ #missing Rg
    ec.s<-ec.s[-which(ec.s$Year==2007),] 
  }
  if(sites[i]=="Vcm"){ #missing Rg
    ec.s<-ec.s[-which(ec.s$Year==2007),] 
  }
  if(sites[i]=="Whs"){ #missing Tair, RH, VPD
    ec.s<-ec.s[-which(ec.s$Year%in%c(2007,2008)),] 
  }

  # Time step: hourly or half-hourly?
  time_step <- ifelse(ec.s$Hour[2]-ec.s$Hour[1]==1, 24, 48) #48 means half-hourly
  
  z <- 5.5 # Threshold (from Papale 2006)... MRJ: SAME as Papale 2006, I think?
  d1 <- rep(NA, dim(ec.s)[1])
  d2 <- rep(NA, dim(ec.s)[1])
  d1[2:dim(ec.s)[1]] <- diff(ec.s[,grep(pattern ='NEE', head)]) #difference between neighboring data (e.g. element 1 = FC_2 - FC_1) 
  d2[1:dim(ec.s)[1]-1] <- diff(ec.s[,grep(pattern ='NEE', head)]) #offset of ^^
  di <- d1-d2 # a difference of differences
  
  # Nighttime
  Md <- median(di[ec.s$Rg < 20 & !is.na(ec.s$Rg)], na.rm=T)
  MAD <- median(abs(di[ec.s$Rg < 20 & !is.na(ec.s$Rg)] - Md), na.rm=T)
  nighttime.outs <- (di[ec.s$Rg < 20 & !is.na(ec.s$Rg)] < (Md - ((z*MAD)/0.6745))) | (di[ec.s$Rg < 20 & !is.na(ec.s$Rg)] > (Md + ((z*MAD)/0.6745)))
  
  # Daytime
  Md <- median(di[ec.s$Rg >= 20 & !is.na(ec.s$Rg)], na.rm=T)
  MAD <- median(abs(di[ec.s$Rg >= 20 & !is.na(ec.s$Rg)] - Md), na.rm=T)
  daytime.outs <- (di[ec.s$Rg >= 20 & !is.na(ec.s$Rg)] < (Md - ((z*MAD)/0.6745))) | (di[ec.s$Rg >= 20 & !is.na(ec.s$Rg)] > (Md + ((z*MAD)/0.6745)))
  
  # Combine
  outs <- rep(FALSE, length(di))
  outs[ec.s$Rg < 20 & !is.na(ec.s$Rg)] <- nighttime.outs
  outs[ec.s$Rg >= 20 & !is.na(ec.s$Rg)] <- daytime.outs
  
  # Remove outliers from NEE in ec.s
  ec.s$NEE[outs] <- NA
  
  ##### Add Posix time stamp - we're still in local standard...
  ec.sp <- fConvertTimeToPosix(ec.s, 'YDH', Year = 'Year', Day = 'DoY', Hour='Hour')
  
  ##### Change nonsensical Rg<0 --> 0 #(already did this)
  #ec.sp$Rg[ec.sp$Rg < 0] <- 0
  
  ##### Process with REddyProc #####
  #Wutzer et al. 2018
  #Basic and extensible post-processing of eddy covariance flux data with REddyProc
  # Initialize
  eddyC <- sEddyProc$new(
    sites[i], ec.sp, 
    c('NEE','LE','Rg','Tair','rH','VPDair',#'Ustar'),DTS=time_step)
      'TsurfEeqn','TsurfE95','TsurfE98',
      'VPDsurfEeqn','VPDsurfE95','VPDsurfE98','Ustar'),DTS=time_step) 
  #NOTE: ^^ This crashes R if the hours are on the .5
  
  # UStar threshold distribution -- this is a little different from what's in the Wutzler paper.
  #In Appendix B, they did uStarTh <- eddyC$sEstUstarThresholdDistribution( 
  uStarTh <- eddyC$sEstimateUstarScenarios(  #This fails with the SRM file I just downloaded... what's the dfference?
    nSample = 100L, probs = c(0.05, 0.5, 0.95))
  uStarThAgg <- eddyC$sGetEstimatedUstarThresholdDistribution()
  eddyC$useSeaonsalUStarThresholds()
  
  # Gap-filling
  eddyC$sMDSGapFillUStarScens("NEE")
  eddyC$sSetLocationInfo(LatDeg = lat, LongDeg = lon, TimeZoneHour = utc)
  #^^I assume this is after NEE gapfill because it is needed for the partitioning, specifically (could have been before)
  
  eddyC$sMDSGapFill('Tair', FillAll = FALSE, minNWarnRunLength = NA)
  
  if(sum(is.na(ec.s$TsurfEeqn))!=nrow(ec.s)){
    eddyC$sMDSGapFill('TsurfEeqn', FillAll = FALSE, minNWarnRunLength = NA)
  }
  if(sum(is.na(ec.s$TsurfE95))!=nrow(ec.s)){
    eddyC$sMDSGapFill('TsurfE95', FillAll = FALSE, minNWarnRunLength = NA)
  }
  if(sum(is.na(ec.s$TsurfE98))!=nrow(ec.s)){
  eddyC$sMDSGapFill('TsurfE98', FillAll = FALSE, minNWarnRunLength = NA)
  }
  
  eddyC$sMDSGapFill('rH', FillAll = FALSE, minNWarnRunLength = NA) 
  eddyC$sMDSGapFill('Rg', FillAll = FALSE, minNWarnRunLength = NA)
  eddyC$sMDSGapFill('VPDair', FillAll = FALSE, minNWarnRunLength = NA)
  
  if(sum(is.na(ec.s$VPDsurfEeqn))!=nrow(ec.s)){
    eddyC$sMDSGapFill('VPDsurfEeqn', FillAll = FALSE, minNWarnRunLength = NA)
  }
  if(sum(is.na(ec.s$VPDsurfE95))!=nrow(ec.s)){
    eddyC$sMDSGapFill('VPDsurfE95', FillAll = FALSE, minNWarnRunLength = NA)
  }
  if(sum(is.na(ec.s$VPDsurfE98))!=nrow(ec.s)){
    eddyC$sMDSGapFill('VPDsurfE98', FillAll = FALSE, minNWarnRunLength = NA)
  }
  eddyC$sMDSGapFill('LE', FillAll = FALSE, minNWarnRunLength = NA)

  # Partition fluxes (nighttime method)
  eddyC$sMRFluxPartitionUStarScens() #also a little different than appendix B in Wutzler...
  
  # Export filled data and add to original file
  FilledEddyData.F <- eddyC$sExportResults()
  FilledEddyData.ToAdd <- FilledEddyData.F [,which(names(FilledEddyData.F)%in%c('NEE_uStar_f', 'NEE_uStar_fqc', 
                                                            'Rg_f', 'Rg_fqc','Tair_f', 'Tair_fqc',
                                                            'TsurfEeqn_f', 'TsurfEeqn_fqc','TsurfE95_f', 'TsurfE95_fqc',
                                                            'TsurfE98_f','TsurfE98_fqc', 'Tair_fqc',  'rH_f', 'rH_fqc',
                                                            'VPDair_f', 'VPDair_fqc','VPDsurfEeqn_f', 'VPDsurfEeqn_fqc',
                                                            'VPDsurfE95_f', 'VPDsurfE95_fqc','VPDsurfE98_f', 'VPDsurfE98_fqc',
                                                            'Reco_uStar', 'GPP_uStar_f', 'GPP_uStar_fqc',
                                                            'LE_f', 'LE_fqc'))]
  names(FilledEddyData.ToAdd)<-gsub("_uStar","",names(FilledEddyData.ToAdd))
  names(FilledEddyData.ToAdd)[which(names(FilledEddyData.ToAdd)=="Reco")]<-"Reco_f" #for consistency
  #names(FilledEddyData.ToAdd) -> headerToAdd
  #headerToAdd[c(1:27)] <- c('NEE_f', 'NEE_fqc',
  #                           'Rg_f', 'Rg_fqc', 'Tair_f', 'Tair_fqc','TsurfEeqn_f', 'TsurfEeqn_fqc',
  #                          'TsurfE95_f', 'TsurfE95_fqc','TsurfE98_f', 'TsurfE98_fqc','rH_f', 'rH_fqc',
  #                          'VPDair_f', 'VPDair_fqc','VPDsurfEeqn_f', 'VPDsurfEeqn_fqc',
  #                          'VPDsurfE95_f', 'VPDsurfE95_fqc','VPDsurfE98_f', 'VPDsurfE98_fqc',
  #                           'Reco_f', 'GPP_f', 'GPP_fqc','LE_f', 'LE_fqc') #changed because now there's rH in there!
  #names(FilledEddyData.ToAdd) <- headerToAdd
  
  `%ni%` <- Negate(`%in%`)
  ec.f <- cbind(ec.sp[,names(ec.sp) %ni% c('DateTime')], FilledEddyData.ToAdd)
  #ec.f <- ec.f[, c('Year','DoY','Hour','Ustar','NEE_f', 'NEE_fqc',
  #                 'Rg_f', 'Rg_fqc', 'Tair_f', 'Tair_fqc', 
  #                 'TsurfEeqn_f','TsurfEeqn_fqc', 'TsurfE95_f','TsurfE95_fqc', 'TsurfE98_f','TsurfE98_fqc',
  #                 'rH_f', 'rH_fqc', 'SWout', 'LWout','WS', 'SWC','albedo','emiss',
  #                 'VPDair_f', 'VPDair_fqc','VPDsurfEeqn_f', 'VPDsurfEeqn_fqc',
  #                 'VPDsurfE95_f', 'VPDsurfE95_fqc','VPDsurfE98_f', 'VPDsurfE98_fqc',
  #                 'Reco_f', 'GPP_f', 'GPP_fqc','LE_f', 'LE_fqc')]
  
  if(sites[i]=="Ha1" | sites[i]=="MMS"){ #put hourly data back.
    ec.f$Hour<-ec.f$Hour-0.5
  }
  
  rm(list = c('ec.s', 'ec.sp','eddyC', 'uStarSuffixes', 'uStarThAgg', 'uStarThAnnual', 'FilledEddyData.F', 'FilledEddyData.ToAdd', 'uStarTh'))
  gc()
  
  #WRITE THE REDDY PROC FILE
  fnout <- paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/C2_ReddyProcOut/US_",sites[i],"_REddyOUT.csv")
  #fnout <- paste0("/Shared/lss_ech2oLab/AMF_CoreSites/C2_ReddyProcOut/US_",sites[i],"_REddyOUT.csv")
  write.table(ec.f, file=fnout, row.names = FALSE, sep = ',') #Writes the REddyProc File
  
}
  
#STEP 1: 
# - Download Ameriflux files
# - Standardize columns in prep for REddyProc

#install.packages("amerifluxr")
library("amerifluxr")
library(bigleaf) #PPFD.to.Rg function
library(REddyProc)
source("/Users/doris/git/HotvDryVPD/date2doy_fxn.R")

# #ONLY HAVE TO RUN THIS PART ONCE:
#  amf_download_base(user_id = "mrjohnst",
#                    user_email = "mjohnston@g.harvard.edu",
#                    site_id = c("US-BRG","US-Ha1","US-Ho1","US-Ho2",
#                                "US-Me6","US-MMS","US-NC2","US-NC3",
#                                "US-NR1","US-Ro4","US-Seg","US-Ses",
#                                "US-SRG","US-SRM","US-Syv","US-Ton",
#                                "US-UMB","US-UMd","US-Var","US-Vcm",
#                                "US-Vcp","US-Vcs","US-WCr","US-Whs",
#                                "US-Wjs","US-Wkg"),
#                    data_product = "BASE-BADM",
#                    data_policy = "CCBY4.0",
#                    agree_policy = TRUE,
#                    intended_use = "synthesis",
#                    intended_use_text = "Research on VPD effects on plant function",
#                    verbose = TRUE,
#                    out_dir = "/Volumes/lss_ech2oLab/AMF_CoreSites/A2_Raw/")


files<-list.files("/Volumes/lss_ech2oLab/AMF_CoreSites/A2_Raw",full.names = TRUE)
filenames<-list.files("/Volumes/lss_ech2oLab/AMF_CoreSites/A2_Raw",full.names = FALSE)

for(y in 1:length(files)){
  print(filenames[y])
    #Read the file
  ec <- amf_read_base(file = files[y],
                      unzip = TRUE,
                      parse_timestamp = FALSE)
  site<-substr(filenames[y],8,10)
  print(site)
#  print(names(ec)[grep("SW",names(ec))])
#  print(names(ec)[grep("LW",names(ec))])
#  print(names(ec)[grep("WS",names(ec))])
#  print(names(ec)[grep("SWC",names(ec))])
#  print("###############################################")
#}

  #Organize times
  ec.time <- as.character(as.matrix(ec$TIMESTAMP_START)) #This is local standard time.
  yr <- as.numeric(substr(ec.time, 1, 4))
  mo <- as.numeric(substr(ec.time, 5, 6))
  dy <- as.numeric(substr(ec.time, 7, 8))
  hr <- as.numeric(substr(ec.time, 9, 10))
  mn <- as.numeric(substr(ec.time, 11, 12))
  dt <- ISOdate(yr, mo, dy, hr, mn)
  doy <- as.numeric(mapply(date2doy, yr=yr, mo=mo, dy=dy))
  ec$yr <- yr
  ec$doy <- doy
  ec$hr <- hr + mn/60 + 0.5  #Ending hour of the half-hourly period
  rm(list = c('yr', 'mo', 'dy', 'hr', 'mn', 'doy', 'ec.time'))
  head <- names(ec)

  #Select standard columns for REddyProc
  #Note that I *really* should have included "P" (precip) here, for the filter...
  #ugh. This is currently done in G_CombiineSMAPandRPOut.R
  
  if(site=="BRG"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT_1_1_1$', head),
                   grep(pattern ='^LW_OUT_1_1_1$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR', head))]
  }
  if(site=="Ha1"){
    ec$SW_IN <- PPFD.to.Rg(ec[,grep(pattern ='^PPFD_IN_1_1_1.*', head)]) #assign a SW_IN var from PPFD, if it doesn't exist
    head <- names(ec)
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR', head))]
    ec.s<-cbind(ec.s[,1:6],as.numeric(rep(NA,nrow(ec.s))),as.numeric(rep(NA,nrow(ec.s))),ec.s[,7:9],
                as.numeric(rep(NA,nrow(ec.s))),ec.s[,10]) #empty cols for SWout, LWout, SWC
  }
  if(site=="Ho1"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT_1_1_1$', head),
                   grep(pattern ='^LW_OUT_1_1_1$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_2_1_1$', head), #1_1_1 absent
                   grep(pattern ='^USTAR_1_1_1$', head))]
  }
  if(site=="Ho2"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT_1_1_1$', head),
                   grep(pattern ='^LW_OUT_1_1_1$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^USTAR_1_1_1$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }

  if(site=="Me6"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="MMS"){
    ec.s <- ec[,  c(grep(pattern ='^yr.*', head),
                    grep(pattern ='^doy.*', head),
                    grep(pattern ='^hr.*', head),
                    grep(pattern ='^FC_1_1_1$', head),
                    grep(pattern ='^LE_1_1_1$', head),
                    grep(pattern ='^SW_IN_1_1_1$', head),
                    grep(pattern ='^SW_OUT_1_1_1$', head),
                    grep(pattern ='^LW_OUT_1_1_1$', head),
                    grep(pattern ='^TA_1_1_1$', head),
                    grep(pattern ='^RH_1_1_1$', head),
                    grep(pattern ='^WS_1_1_1$', head),
                    grep(pattern ='^SWC_1_1_1$', head),
                    grep(pattern ='^USTAR.*', head))]
  }
  if(site=="NC2"){
    head <- names(ec)
    ec.s <- ec[,  c(grep(pattern ='^yr.*', head),
                    grep(pattern ='^doy.*', head),
                    grep(pattern ='^hr.*', head),
                    grep(pattern ='^FC_1_1_1$', head),
                    grep(pattern ='^LE_1_1_1$', head),
                    grep(pattern ='^SW_IN_1_1_1$', head),
                    grep(pattern ='^SW_OUT_1_1_1$', head),
                    grep(pattern ='^LW_OUT_1_1_1$', head),
                    grep(pattern ='^TA_1_1_1$', head),
                    grep(pattern ='^RH_1_1_1$', head),
                    grep(pattern ='^WS_1_1_1$', head),
                    grep(pattern ='^SWC_1_1_1$', head),
                    grep(pattern ='^USTAR.*', head))]
  }
  if(site=="NC3"){
    ec.s <- ec[,  c(grep(pattern ='^yr.*', head),
                    grep(pattern ='^doy.*', head),
                    grep(pattern ='^hr.*', head),
                    grep(pattern ='^FC_1_1_1$', head),
                    grep(pattern ='^LE_1_1_1$', head),
                    grep(pattern ='^SW_IN_1_1_1$', head),
                    grep(pattern ='^SW_OUT_1_1_1$', head),
                    grep(pattern ='^LW_OUT_1_1_1$', head),
                    grep(pattern ='^TA_1_1_1$', head),
                    grep(pattern ='^RH_1_1_1$', head),
                    grep(pattern ='^WS_1_1_1$', head),
                    grep(pattern ='^SWC_1_1_1$', head),
                    grep(pattern ='^USTAR.*', head))]
  }
  if(site=="NR1"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT_1_1_1$', head),
                   grep(pattern ='^LW_OUT_1_1_1$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR_1_1_1$', head))]
  }
  if(site=="Ro4"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^SWC$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="Seg"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }
  if(site=="Ses"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }
  if(site=="SRG"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="SRM"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_PI_1_1_A$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="Syv"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT_1_1_1$', head),
                   grep(pattern ='^LW_OUT_1_1_1$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR_1_1_1$', head))]
  }
  if(site=="Ton"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT_1_1_1$', head),
                   grep(pattern ='^LW_OUT_1_1_1$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_PI_1_1_A$', head),
                   grep(pattern ='^USTAR_1_1_1*', head))]
  }
  if(site=="UMB"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="UMd"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="Var"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN_1_1_1$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^SWC_PI_1_1_A$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="Vcm"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }
  if(site=="Vcp"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }
  if(site=="Vcs"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }
  if(site=="WCr"){
    ec$avp<-fCalcAVPfromVMFandPress(VMF=ec$H2O_1_1_1/1000, Press=ec$PA_1_1_1*10) 
    ec$RH<-fCalcRHfromAVPandTair(AVP = ec$avp, Tair=ec$TA_1_1_1)
    head <- names(ec)
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC_1_1_1$', head),
                   grep(pattern ='^LE_1_1_1$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR_1_1_1', head))]
  }
  if(site=="Whs"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  if(site=="Wjs"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA$', head),
                   grep(pattern ='^RH$', head),
                   grep(pattern ='^WS$', head),
                   grep(pattern ='^USTAR$', head))]
    ec.s<-cbind(ec.s[,1:11],as.numeric(rep(NA,nrow(ec.s))),ec.s[,12]) #empty col for SWC
  }
  if(site=="Wkg"){
    ec.s <- ec[, c(grep(pattern ='^yr.*', head),
                   grep(pattern ='^doy.*', head),
                   grep(pattern ='^hr.*', head),
                   grep(pattern ='^FC$', head),
                   grep(pattern ='^LE$', head),
                   grep(pattern ='^SW_IN$', head),
                   grep(pattern ='^SW_OUT$', head),
                   grep(pattern ='^LW_OUT$', head),
                   grep(pattern ='^TA_1_1_1$', head),
                   grep(pattern ='^RH_1_1_1$', head),
                   grep(pattern ='^WS_1_1_1$', head),
                   grep(pattern ='^SWC_1_1_1$', head),
                   grep(pattern ='^USTAR$', head))]
  }
  #Standardize naming
      names(ec.s) <- c("Year","DoY","Hour","NEE","LE","Rg",
                       "SWout","LWout","Tair","rH","WS","SWC","Ustar")
      ec.s$NEE[which(ec.s$NEE<c(-999))]<-NA
      ec.s$LE[which(ec.s$LE<c(-999))]<-NA
      ec.s$Rg[which(ec.s$Rg<0)]<-0 
      ec.s$SWout[which(ec.s$SWout<0)]<-NA
      ec.s$LWout[which(ec.s$LWout<0)]<-NA 
      ec.s$Tair[which(ec.s$Tair<c(-30))]<-NA 
      ec.s$Tair[which(ec.s$Tair>60)]<-NA    #changed to match Tsurf - was 50
      ec.s$rH[which(ec.s$rH>100)]<-NA
      ec.s$rH[which(ec.s$rH<0)]<-NA 
      ec.s$WS[which(ec.s$WS<0)]<-NA 
      ec.s$SWC[which(ec.s$SWC<0)]<-NA 
 
    #Calculations     
      ec.s$albedo<- ec.s$SWout/ec.s$Rg
      ec.s$albedo[which(ec.s$Rg==0)]<-NA
      ec.s$albedo[which(ec.s$albedo>1)]<-NA
      ec.s$emiss<-0.99-(0.16*ec.s$albedo)
        SB =0.0000000567 
      ec.s$TsurfEeqn<-((ec.s$LWout/(SB*ec.s$emiss))^0.25)-273.15 #with emissivity= -0.16*albedo +0.99
        ec.s$TsurfEeqn[which(ec.s$TsurfEeqn<c(-30))]<-NA 
        ec.s$TsurfEeqn[which(ec.s$TsurfEeqn>60)]<-NA 
      ec.s$TsurfE95<-((ec.s$LWout/(SB*0.95))^0.25)-273.15 #with emissivity= 0.95
        ec.s$TsurfE95[which(ec.s$TsurfE95<c(-30))]<-NA 
        ec.s$TsurfE95[which(ec.s$TsurfE95>60)]<-NA
      ec.s$TsurfE98<-((ec.s$LWout/(SB*0.98))^0.25)-273.15 #with emissivity= 0.98
        ec.s$TsurfE98[which(ec.s$TsurfE98<c(-30))]<-NA 
        ec.s$TsurfE98[which(ec.s$TsurfE98>60)]<-NA
      #VPD  
      ec.s$VPDair <- fCalcVPDfromRHandTair(ec.s$rH, ec.s$Tair) #calculate the vpd, if we don't have it.
        #Note: the warning here ^^ is because there are missing values in RH and Temp.
      ec.s$VPDsurfEeqn <- fCalcVPDfromRHandTair(ec.s$rH, ec.s$TsurfEeqn) #calculate the vpd, if we don't have it.
      ec.s$VPDsurfE95 <- fCalcVPDfromRHandTair(ec.s$rH, ec.s$TsurfE95) #calculate the vpd, if we don't have it.
      ec.s$VPDsurfE98 <- fCalcVPDfromRHandTair(ec.s$rH, ec.s$TsurfE98) #calculate the vpd, if we don't have it.
      
      #Write file
      write.csv(ec.s,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/B2_ReddyProcIn/US_",site,"_REddyIN.csv"))

}


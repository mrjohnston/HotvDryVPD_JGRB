#Adapted/cleaned from FilterAMFCore_v3.R; MRJ 6/7/23
#There are several kinds of filters here:
  #[0] Based on the quality of the fill (VPD, GPP, LE, NEE, all the air temp stuff)
  #[1] Sanity values & 3sd filtes
  #[2] Wet filter
  #[3] Wind filter
  #[4] Growing season and time of day filter
#I'm going to write out a new spreadsheet with filters applied, for simplicity.
library(zoo)
library(prospectr)
library(suncalc)
library(lutz)
info<-read.csv("/Volumes/lss_ech2oLab/AMF_CoreSites/FilteredSiteInfo_v3.csv")

setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/G2.4_AntecedentConditions/")
`%nin%` = Negate(`%in%`)

files<-list.files()
wetfilterlog<-data.frame(matrix(nrow=26,ncol=7))
names(wetfilterlog)<-c("site","ntot","lossP24","lossP48","lossSWClag1","lossSWClag2","lossSMAP")

for(i in 1:length(files)){
  dat<-read.csv(files[i])
  site<-substr(files[i],1,6)
  print(paste0(site, "; ",i, " / ", length(files)))
  
  wetfilterlog[i,1]<-site
  wetfilterlog[i,2]<-nrow(dat)
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
#[0] Fill quality
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  fillquality<-function(dataframe,filledvar,fillquality){
    #filledvar and fillquality are column names (character)
    #Get the relevant column numers
    filledvarnum<-which(names(dataframe)==filledvar)
    fillqualitynum<-which(names(dataframe)==fillquality)
    #initialize newvar
    newvar<-rep(NA,length=nrow(dataframe))
    #fill with high-quality results only
    newvar[which(dataframe[,fillquality]%in%c(0,1))]<-dataframe[,filledvar][which(dataframe[,fillquality]%in%c(0,1))]
    return(newvar)
  }
  
  NEEmod<-fillquality(dataframe=dat, filledvar = "NEE_f", fillquality = "NEE_fqc")
  GPPmod<-fillquality(dataframe=dat, filledvar = "GPP_f", fillquality = "GPP_fqc")
  LEmod<-fillquality(dataframe=dat, filledvar = "LE_f", fillquality = "LE_fqc")
  
  VPDairmod<-fillquality(dataframe=dat, filledvar = "VPDair_f", fillquality = "VPDair_fqc")
  Tairmod<-fillquality(dataframe=dat, filledvar = "Tair_f", fillquality = "Tair_fqc")
  rH<-fillquality(dataframe=dat, filledvar = "rH_f", fillquality = "rH_fqc")
  
  #Note: Ha1 doesn't report LWout, so doesn't have the surface stuff...
  if(site!="US_Ha1"){
  VPDsurfEeqnmod<-fillquality(dataframe=dat, filledvar = "VPDsurfEeqn_f", fillquality = "VPDsurfEeqn_fqc")
  VPDsurfEe98mod<-fillquality(dataframe=dat, filledvar = "VPDsurfE98_f", fillquality = "VPDsurfE98_fqc")
  VPDsurfEe95mod<-fillquality(dataframe=dat, filledvar = "VPDsurfE95_f", fillquality = "VPDsurfE95_fqc")
  TsurfEeqnmod<-fillquality(dataframe=dat, filledvar = "TsurfEeqn_f", fillquality = "TsurfEeqn_fqc")
  TsurfEe98mod<-fillquality(dataframe=dat, filledvar = "TsurfE98_f", fillquality = "TsurfE98_fqc")
  TsurfEe95mod<-fillquality(dataframe=dat, filledvar = "TsurfE95_f", fillquality = "TsurfE95_fqc")
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#[1] Sanity checks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  #NEEmod #no check here
  GPPmod[c(which(GPPmod<0))]<-NA
  LEmod[c(which(LEmod<0))]<-NA
  VPDairmod[c(which(VPDairmod<0))]<-NA
  Tairmod[c(which(Tairmod>50), which(Tairmod<c(-30)))]<-NA
  rH[c(which(rH>100), which(rH<0))]<-NA
  
  if(site!="US_Ha1"){
  VPDsurfEeqnmod[c(which(VPDsurfEeqnmod<0))]<-NA
  VPDsurfEe98mod[c(which(VPDsurfEe98mod<0))]<-NA
  VPDsurfEe95mod[c(which(VPDsurfEe95mod<0))]<-NA
  TsurfEeqnmod[c(which(TsurfEeqnmod>50), which(TsurfEeqnmod<c(-30)))]<-NA
  TsurfEe98mod[c(which(TsurfEe98mod>50), which(TsurfEe98mod<c(-30)))]<-NA
  TsurfEe95mod[c(which(TsurfEe95mod>50), which(TsurfEe95mod<c(-30)))]<-NA
  }
  
  #Remove points >< 3 SD from the mean
  VPDairmod[c(which(VPDairmod < (mean(VPDairmod,na.rm=T)-(3*sd(VPDairmod,na.rm=T)))),
              which(VPDairmod > (mean(VPDairmod,na.rm=T)+(3*sd(VPDairmod,na.rm=T)))))]<-NA 
  VPDsurfEeqnmod[c(which(VPDsurfEeqnmod < (mean(VPDsurfEeqnmod,na.rm=T)-(3*sd(VPDsurfEeqnmod,na.rm=T)))),
                   which(VPDsurfEeqnmod > (mean(VPDsurfEeqnmod,na.rm=T)+(3*sd(VPDsurfEeqnmod,na.rm=T)))))]<-NA
  VPDsurfEe98mod[c(which(VPDsurfEe98mod < (mean(VPDsurfEe98mod,na.rm=T)-(3*sd(VPDsurfEe98mod,na.rm=T)))),
                   which(VPDsurfEe98mod > (mean(VPDsurfEe98mod,na.rm=T)+(3*sd(VPDsurfEe98mod,na.rm=T)))))]<-NA
  VPDsurfEe95mod[c(which(VPDsurfEe95mod < (mean(VPDsurfEe95mod,na.rm=T)-(3*sd(VPDsurfEe95mod,na.rm=T)))),
                   which(VPDsurfEe95mod > (mean(VPDsurfEe95mod,na.rm=T)+(3*sd(VPDsurfEe95mod,na.rm=T)))))]<-NA
  
  Tairmod[c(which(Tairmod < (mean(Tairmod,na.rm=T)-(3*sd(Tairmod,na.rm=T)))),
            which(Tairmod > (mean(Tairmod,na.rm=T)+(3*sd(Tairmod,na.rm=T)))))]<-NA 
  TsurfEeqnmod[c(which(TsurfEeqnmod < (mean(TsurfEeqnmod,na.rm=T)-(3*sd(TsurfEeqnmod,na.rm=T)))),
                 which(TsurfEeqnmod > (mean(TsurfEeqnmod,na.rm=T)+(3*sd(TsurfEeqnmod,na.rm=T)))))]<-NA
  TsurfEe98mod[c(which(TsurfEe98mod < (mean(TsurfEe98mod,na.rm=T)-(3*sd(TsurfEe98mod,na.rm=T)))),
                 which(TsurfEe98mod > (mean(TsurfEe98mod,na.rm=T)+(3*sd(TsurfEe98mod,na.rm=T)))))]<-NA
  TsurfEe95mod[c(which(TsurfEe95mod < (mean(TsurfEe95mod,na.rm=T)-(3*sd(TsurfEe95mod,na.rm=T)))),
                 which(TsurfEe95mod > (mean(TsurfEe95mod,na.rm=T)+(3*sd(TsurfEe95mod,na.rm=T)))))]<-NA
  
  rH[c(which(rH < (mean(rH,na.rm=T)-(3*sd(rH,na.rm=T)))),
            which(rH > (mean(rH,na.rm=T)+(3*sd(rH,na.rm=T)))))]<-NA 
  
  #Make sure that SMAP is NA when no measurements.
  if(site%in%c("US_Ha1","US_MMS")) { #hourly
    firstrowSMAP<-min(which(dat$SMAP_rootzoneSM==unique(dat$SMAP_rootzoneSM)[2]))-3  }
  if(site%nin%c("US_Ha1","US_MMS")) { #half-hourly 
    firstrowSMAP<-min(which(dat$SMAP_rootzoneSM==unique(dat$SMAP_rootzoneSM)[2]))-6  } #aka 6 timesteps before the second unique value
  #sprintf("%.10f",unique(dat$SMAP_rootzoneSM)[2])
  #sprintf("%.10f",unique(dat$SMAP_rootzoneSM)[1])
  if(firstrowSMAP<1){
    firstrowSMAP=1}
  SMAP_rootzoneSM<-dat$SMAP_rootzoneSM
  if(firstrowSMAP!=1){
    SMAP_rootzoneSM[1:c(firstrowSMAP-1)]<-NA
  }
  SMAP_rootzoneSM[c(which(SMAP_rootzoneSM < (mean(SMAP_rootzoneSM,na.rm=T)-(3*sd(SMAP_rootzoneSM,na.rm=T)))),
                    which(SMAP_rootzoneSM > (mean(SMAP_rootzoneSM,na.rm=T)+(3*sd(SMAP_rootzoneSM,na.rm=T)))))]<-NA
  
  NEEmod[c(which(NEEmod < (mean(NEEmod,na.rm=T)-(3*sd(NEEmod,na.rm=T)))),
           which(NEEmod > (mean(NEEmod,na.rm=T)+(3*sd(NEEmod,na.rm=T)))))]<-NA
  GPPmod[c(which(GPPmod < (mean(GPPmod,na.rm=T)-(3*sd(GPPmod,na.rm=T)))),
           which(GPPmod > (mean(GPPmod,na.rm=T)+(3*sd(GPPmod,na.rm=T)))))]<-NA 
  LEmod[c(which(LEmod < (mean(LEmod,na.rm=T)-(3*sd(LEmod,na.rm=T)))),
          which(LEmod > (mean(LEmod,na.rm=T)+(3*sd(LEmod,na.rm=T)))))]<-NA 
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# [2] Windfilter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  windfilter<-vector(mode="character",length=nrow(dat))
  windfilter[which(dat$WS>1)]<-"keep"
  windfilter[which(dat$WS<1)]<-"toss"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# [3] Wetfilter, with notes about how much data loss in each step
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  wetfilter<-rep("keep",length=nrow(dat))

  wetfilter[which(dat$precip_sums24hr>0)]<-"toss" #if there's any rain in the past 24 hr
  wetfilterlog[i,3]<-table(wetfilter)[2]
  wetfilter[which(dat$precip_sums48hr>5)]<-"toss" #if there's >5mm rain in the past 48 hr
  wetfilterlog[i,4]<-table(wetfilter)[2]
  
  wetfilter[which(dat$swc_diff_lag1>0.02)]<-"toss"
  wetfilterlog[i,5]<-table(wetfilter)[2]
  wetfilter[which(dat$swc_diff_lag2>0.02)]<-"toss" #This is an hour for the HH sites;
  wetfilterlog[i,6]<-table(wetfilter)[2]
  #swc_diff_lag2 = NA for the hourly sites
  
  #Make "wetfilter" also = 1 when SMAP_surfaceSM_diff > 0.02
  #Note that there has to be a contingency if there isn't variable SMAP data for the site (i.e. site ends before March 2015)
 #This snippet really should have been in the antecedent code...
  if(length(unique(dat$SMAP_posix_utc))>1){
    startSMdiff<-min(which(dat$SMAP_posix_utc==unique(dat$SMAP_posix_utc)[2])) #Row with the second unique date
    smapreps<-length(which(dat$SMAP_posix_utc==unique(dat$SMAP_posix_utc)[2])) #should be 6 for half-hourly data; 3 for hourly.
    
    #Get lagged data, by smapreps
    dat$SMAP_surf_diff<-NA
    dat$SMAP_surf_diff<-c(rep(NA,smapreps),diff(dat[,grep("SMAP_surfaceSM",names(dat))],lag=smapreps,differences = 1))
    
    #Fix the first smapreps NA values...
    if(startSMdiff<smapreps){ #this is when the record starts after SMAP
      firstbit<-c(rep(NA,startSMdiff-1), #rows where I don't have a baseline # for the difference - should actually be NA
                  rep(dat$SMAP_surfaceSM[startSMdiff]-dat$SMAP_surfaceSM[startSMdiff-1],smapreps-(startSMdiff-1))) #value for the next NAs
      dat$SMAP_surf_diff[1:length(firstbit)]<-firstbit
    }
    
    if(startSMdiff>smapreps){ #this happens when there's time in the record before SMAP
      firstbit<-c(rep(NA,startSMdiff-1)) #rows where I don't have a baseline # for the difference - should be NA
      dat$SMAP_surf_diff[1:length(firstbit)]<-firstbit
    }
  }
    wetfilter[which(dat$SMAP_surf_diff>0.02)]<-"toss"
    wetfilterlog[i,7]<-table(wetfilter)[2]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# [4] Growing time filter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  #Consider any time with sun angle at 0.5 radians or above 
  #(see email "solar angle vs. potential PPFD" 9/23/22)
  sunaltthresh<-0.50 #include HH where the sun altitude is >0.5 radians
  datefun <- function(DF){ #This is going to be used to calc solar angle
    #https://stackoverflow.com/questions/66483922/how-to-combine-
    #year-day-of-the-year-and-hour-columns-to-datetime-in-r
    d <- with(DF, paste(Year, DoY))
    d <- as.Date(d, "%Y %j")
    hm <- DF[["Hour"]]*60
    d <- paste(d, paste(hm %/% 60, hm %% 60, 0, sep = ":"))
    d <- as.POSIXct(d, "%Y-%m-%d %H:%M:%S",tz="UTC")
    return(d)
  }
  idx<-which(info$site==substr(files[i],1,6))
  tz <- tz_lookup_coords(lat=info$lat[idx],
                         lon=info$lon[idx],
                         method = "fast", warn = FALSE)
  d_utc<-datefun(dat) - (tz_offset("2018-01-01", tz)$utc_offset_h * 60 * 60)
  altiDF<-getSunlightPosition(date = d_utc, 
                              lat = info$lat[idx],
                              lon = info$lon[idx], keep = "altitude")
  sunaltitude<-altiDF$altitude
  growtimefilter<-vector(mode="character",length=nrow(dat))
  growtimefilter[which(sunaltitude>=0.5)]<-"keep"
  growtimefilter[which(sunaltitude<0.5)]<-"toss"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  #[5] Growing season
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#Average by DOY - GPP
  dat$GPPmod<-GPPmod
seas<-aggregate(dat$GPPmod[which(!is.na(dat$GPPmod))],
                by=list(dat$DoY[which(!is.na(dat$GPPmod))]),FUN=mean)

#Get a rolling average:
rollingint<-45 #how long for rolling average?
x=0.8 #What values to keep (this is top 50%)

#Savitzky-Golay:
sg<-savitzkyGolay(seas$x,m=0,p=3,w=rollingint) #w = window size; p = polynomial degree
seas$sg<-c(rep(NA,(nrow(seas)-length(sg))/2),
           sg,rep(NA,(nrow(seas)-length(sg))/2))

#Take the top x%
cutoffval<-max(seas$sg,na.rm=T)*x 

#What are the DOY brackets?
growingseas<-seas$Group.1[which(seas$sg>cutoffval)] #This is the growing season. < for NEE

#the filter -- FILTER OUT 1, as with other filters
growseasfilter<-rep("toss",length=nrow(dat))
growseasfilter[which(dat$DoY%in%growingseas)]<-"keep"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#MAKE A FILTERED SPREADSHEET
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#Note that will always want high-quality & only growing season & growing time, but 
#I will not apply the wetness or wind filters for the first analysis
if(site=="US_Ha1"){
  VPDsurfEeqnmod<-rep(NA,length=nrow(dat))
  VPDsurfEe98mod<-rep(NA,length=nrow(dat))
  VPDsurfEe95mod<-rep(NA,length=nrow(dat))
  TsurfEeqnmod<-rep(NA,length=nrow(dat))
  TsurfEe98mod<-rep(NA,length=nrow(dat))
  TsurfEe95mod<-rep(NA,length=nrow(dat))
}
datF<-data.frame(dat$Year,dat$DoY,dat$Hour,
            NEEmod,GPPmod,LEmod,SMAP_rootzoneSM,
            VPDairmod,VPDsurfEeqnmod, VPDsurfEe98mod,VPDsurfEe95mod,
            Tairmod, TsurfEeqnmod, TsurfEe98mod,TsurfEe95mod,rH,
            windfilter, wetfilter, growtimefilter, growseasfilter)
write.csv(datF,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/H2.fin_Filtered_withRH/",site,"_filt.csv"))
}

#write.csv(wetfilterlog,"/Volumes/lss_ech2oLab/AMF_CoreSites/wetfilterlog.csv")

# #Looking over wetfilterlog
# #Compare proportion tossed with and without SMAP
# with<-wetfilterlog$lossSMAP/wetfilterlog$ntot
# wo<-wetfilterlog$lossSWClag2/wetfilterlog$ntot
# 
# with-wo
# summary(with-wo) # median 0.19%; max 1.5% - these are DIFFERENCES, e.g. % of total data omitted because
#   #of SMAP that wouldn't have otherwise been omitted.
#   #I think the take-away here is that I don't need to clip starting at 2015, which is amazzzzeballs.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#NOTES
#[6] Novick filters ############################
# [a - I don't have ET] ET: data were first filtered to remove missing data and 
#ET measurements greater than 2 mm/hour or less than 0.6 mm/hour. 
# [b] Meteorological data were also screened for obvious outliers
# (i.e. air temperature less than -30oC or greater than 50oC,
# VPD < 0, net radiation less than -500 W m-2 or greater than 1500 W m-2).

#dat$NovickQC_filter[which(dat$Rg<0)]<-1 #this is incoming solar, not net!
#dat$NovickQC_filter[which(dat$Rg>2000)]<-1
# [c] Most of the analysis was limited to daytime periods (Rn > 50 W m-2)
# [d] when wind speed exceeded 1 m s-1 and VPD >0.6 kPa (to minimize stability effects).
# ^ 'most' means not actually filtered? But I might want to report on this in 'taking stock'

#dat$NovickQC_filter[which(dat$WS<1)]<-1
#dat$NovickQC_filter[which(dat$VPD_filtered<6)]<-1 #6 rather than .6 because hPa vs kPa
# [e] Also, filtered to periods of stationary LAI: we filtered our data 
# to exclude periods when leaf area was expected to be low or highly dynamic
# (based on temperatures, seems rough. In Mediterranean sites, based on 
# literature or communicating with PI -- I like the GPP way.)
# [f] When parameterizing the model of Eq. 1 in the main text, we
# limited the data to those collected when VPD > 1.0 kPa.
# ^ because of calculating Gs issues; I think that 

#Also Zhang et al. 2019:
#To minimize the contribution of surface evaporation to iWUE,
#and to control for confounding effects from variable leaf area,
#radiation, and windspeed, the data were subjected to the following
#set of filters: 
#(1) we only used original measurements (quality control flag=0) or
#gap filled data of good quality (quality control flag=) as determined by the FLUXNET2015 protocols;
#(2) the analysis was conducted for periods of relatively stationary Normalized Difference Vegetation Index (NDVI)
#(https://modis.ornl.gov/ fixedsite/) with NDVI> 80% NDVImax (the peak of the annual NDVI series).
#The NDVI product has a spatial resolution of 250 m, and we used the pixel where the tower resides;
#(3) the analysis was constrained to 10:00–15:00 local time and only when incoming 
#short-wave radiation (Rs) is higher than 50% of the peak of its annual series (Rsmax),
#i.e. Rs> 50% Rsmax. Note that a discussion of uncertainty introduced by the selection
#of these two thresholds is provided later in this document; 
#(4) data were filtered to exclude records within two days following rain events (>1 mm)
#to minimize the contamination of rain on LE measurements; 
#(5) the wind speed was constrained to exceed 1 m s−1 to select for periods with negligible 
#leaf boundary layer resistance.

#Write out the file
#write.csv(dat,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/H_Filtered/",site,"_filters_v3.csv"))
#write.csv(dat,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/H2_Filtered/",site,"_filters_v4.csv"))
      #^ change to v4 when add thesun angle filter

  #pre-Tuscon filters:
# - growing season, growing time > 0.5 (post-Tuscon: 0.8)
# - precip filter only if rain in last 24 hr (post-Tuscon: also filter if >5mm in last 48 hr)
# - SWC filter if next timestep >0.2 previous (post-Tuscon: >0.02)
# - Did not include WS > 1m/s, VPD > 6hPa, PPFD > 800 µmolPhoton m-2 s-1 (post-Tuscon: includes these)
                     





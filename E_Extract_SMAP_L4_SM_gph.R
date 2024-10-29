# Extract SMAP L4 geophysical soil moisture (surface, rootzone - version 7) 
# Written by M. Johnston 02/2022; Checked, edited for v7, and rerun 05/2023
# Run on Argon: /Shared/lss_ech2oLab/AMF_CoreSites/E2_SMAPOut/Extract_groupX_v2.job -->
# Calls Extract_groupX_v2.R (below - this is the version on the Argon cluster).
# Results are csv files in /Shared/lss_ech2oLab/AMF_CoreSites/E2_SMAPOut/ExtractedData/

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rhdf5")
library(rhdf5)
library(raster)

#Loop through sites/lats/longs
#"Info" is a csv with columns Site (e.g. "US-BRG"), Latitude (decimal degrees), Longitude (decimal degrees).
#I extract 5-6 sites in each job.
info<-read.csv("/Shared/lss_ech2oLab/AMF_CoreSites/E2_SMAPOut/AMFCore_group1_v2.csv") #This changes
site = info$Site
inlat = info$Latitude
inlon = info$Longitude

#List files from which to extract
setwd("/Shared/lss_ech2oLab/SMAP_L4_SM_v7_gph")
fns = list.files(pattern="h5$")

#Initialize empty vectors to hold data
for(j in 1:length(site)){
  print(paste0("Site: ",site[j]))
  yr<-month<-day<-doy<-CompletenessOmission_QC<-DomainConsistency_QC<-
    surfSM<-rzSM<-vector(mode="numeric",length=length(fns))
  time<-vector(mode="character",length=length(fns))
  
  #It seems that all the SMAP grids have the same extent:
  lon <- h5read(fns[1], "cell_lon")
  lat <- h5read(fns[1], "cell_lat")
  #lon2 <- h5read(fns[length(fns)], "cell_lon")
  #lat2 <- h5read(fns[length(fns)], "cell_lat")
  #all.equal(lon,lon2)
  #all.equal(lat,lat2)
  h5closeAll()
  
  #So, just read in the lat/long once (above) to save computation time.
  #Get the index of the lat/lon that I want, so I can extract a single surface & rootzone SM
  lonlat<-cbind(as.vector(lon),as.vector(lat))
  PD<-pointDistance(p1=lonlat,p2=cbind(inlon[j],inlat[j]),lonlat=TRUE)
  nearlon<-lonlat[which(PD==min(PD)),1]
  nearlat<-lonlat[which(PD==min(PD)),2]
  select_index_row<-which(lon==nearlon,arr.ind = TRUE)[1,1]
  select_index_col<-which(lat==nearlat,arr.ind = TRUE)[1,2]
  
  for(i in 1:length(fns)){
    #for(i in c(1,100, 1000, 5000, 10000, 19984)){ #Looks like it's always the same grid.
    #for(i in 1:10){ #for code-testing purposes
    file<-fns[i]
    print(file) #just keeping track.
    
    #Get the data at these lat/longs & put it into a dataframe
    yr[i]<-as.numeric(substr(fns[i],16,19))
    month[i]<-as.numeric(substr(fns[i],20,21))
    day[i]<-as.numeric(substr(fns[i],22,23))
    doy[i]<- strftime(paste(c(yr[i],'-',month[i],'-',day[i]), collapse=''),format='%j')
    time[i]<-paste0(substr(fns[i],25,26),":",substr(fns[i],27,28))
    co<-h5readAttributes(fns[i], name = "/Metadata/DataQuality/SM/CompletenessOmission/")
    dc<-h5readAttributes(fns[i], name = "/Metadata/DataQuality/SM/DomainConsistency/")
    CompletenessOmission_QC[i]<-co$value
    DomainConsistency_QC[i]<-dc$value
    rzSM[i]<-h5read(file,"/Geophysical_Data/sm_rootzone",
                    index=list(select_index_row,select_index_col))
    surfSM[i]<-h5read(file,"/Geophysical_Data/sm_surface",
                      index=list(select_index_row,select_index_col))
    h5closeAll()
  }
  
  datafile<-data.frame(rep(site[j],length(fns)),rep(inlat[j],length(fns)),
                       rep(inlon[j],length(fns)),yr,month,day,doy,time,
                       CompletenessOmission_QC,DomainConsistency_QC,rzSM,surfSM)
  names(datafile)<-c("site","lat","lon","year","month","day","doy","time_intervalcenter_utc",
                     "CompletenessOmission_QC","DomainConsistency_QC","rootzone_SM","surface_SM")
  
  write.csv(datafile,paste0("/Shared/lss_ech2oLab/AMF_CoreSites/E2_SMAPOut/ExtractedData/",site[j],"_","SMAP_L4_SM_v7_gph_withQC.csv"))
  
  rm(datafile)
}
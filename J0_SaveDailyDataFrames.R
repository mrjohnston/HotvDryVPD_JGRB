setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/I2_TadjDF/SVP/")
files<-list.files(pattern="*Tadj.csv*")

for(z in 1:length(files)){ 
  dat<-read.csv(files[z])
  print(paste0("SITE:      ",files[z],"        #",z))
  print(dim(dat))
  names(dat)[which(names(dat)=="Group.2")]<-"dat.Year"
  
  if(outlierfilter==TRUE){
    dat$NEEmod[c(which(dat$NEEmod < (mean(dat$NEEmod,na.rm=T)-(3*sd(dat$NEEmod,na.rm=T)))),
                 which(dat$NEEmod > (mean(dat$NEEmod,na.rm=T)+(3*sd(dat$NEEmod,na.rm=T)))))]<-NA
    dat$GPPmod[c(which(dat$GPPmod < (mean(dat$GPPmod,na.rm=T)-(3*sd(dat$GPPmod,na.rm=T)))),
                 which(dat$GPPmod > (mean(dat$GPPmod,na.rm=T)+(3*sd(dat$GPPmod,na.rm=T)))))]<-NA 
    dat$LEmod[c(which(dat$LEmod < (mean(dat$LEmod,na.rm=T)-(3*sd(dat$LEmod,na.rm=T)))),
                which(dat$LEmod > (mean(dat$LEmod,na.rm=T)+(3*sd(dat$LEmod,na.rm=T)))))]<-NA 
    
    dat$svpAir[c(which(dat$svpAir < (mean(dat$svpAir,na.rm=T)-(3*sd(dat$svpAir,na.rm=T)))),
                 which(dat$svpAir > (mean(dat$svpAir,na.rm=T)+(3*sd(dat$svpAir,na.rm=T)))))]<-NA 
    dat$Tadj_svp_air[c(which(dat$Tadj_svp_air < (mean(dat$Tadj_svp_air,na.rm=T)-(3*sd(dat$Tadj_svp_air,na.rm=T)))),
                       which(dat$Tadj_svp_air > (mean(dat$Tadj_svp_air,na.rm=T)+(3*sd(dat$Tadj_svp_air,na.rm=T)))))]<-NA 
    dat$VPDairmod[c(which(dat$VPDairmod < (mean(dat$VPDairmod,na.rm=T)-(3*sd(dat$VPDairmod,na.rm=T)))),
                    which(dat$VPDairmod > (mean(dat$VPDairmod,na.rm=T)+(3*sd(dat$VPDairmod,na.rm=T)))))]<-NA 
    dat$SMAP_rootzoneSM[c(which(dat$SMAP_rootzoneSM < (mean(dat$SMAP_rootzoneSM,na.rm=T)-(3*sd(dat$SMAP_rootzoneSM,na.rm=T)))),
                          which(dat$SMAP_rootzoneSM > (mean(dat$SMAP_rootzoneSM,na.rm=T)+(3*sd(dat$SMAP_rootzoneSM,na.rm=T)))))]<-NA 
    
    if(z!=2){
      dat$svpLeaf[c(which(dat$svpLeaf < (mean(dat$svpLeaf,na.rm=T)-(3*sd(dat$svpLeaf,na.rm=T)))),
                    which(dat$svpLeaf > (mean(dat$svpLeaf,na.rm=T)+(3*sd(dat$svpLeaf,na.rm=T)))))]<-NA 
      dat$Tadj_svp_leaf[c(which(dat$Tadj_svp_leaf < (mean(dat$Tadj_svp_leaf,na.rm=T)-(3*sd(dat$Tadj_svp_leaf,na.rm=T)))),
                          which(dat$Tadj_svp_leaf > (mean(dat$Tadj_svp_leaf,na.rm=T)+(3*sd(dat$Tadj_svp_leaf,na.rm=T)))))]<-NA 
      dat$VPDsurfEeqnmod[c(which(dat$VPDsurfEeqnmod < (mean(dat$VPDsurfEeqnmod,na.rm=T)-(3*sd(dat$VPDsurfEeqnmod,na.rm=T)))),
                           which(dat$VPDsurfEeqnmod > (mean(dat$VPDsurfEeqnmod,na.rm=T)+(3*sd(dat$VPDsurfEeqnmod,na.rm=T)))))]<-NA 
    }
    dat<-dat[complete.cases(dat),]
  }
  dat$dat.Year<-as.factor(dat$dat.Year)
  write.csv(dat,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/J_DailyFiles/",files[z]))
}

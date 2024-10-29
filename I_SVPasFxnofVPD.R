#SVP as fxn of VPD - this should be linear, easier to fit, less error prone

cc<-function(temp){
  svp<-6.1078*(exp((17.08085*temp)/(temp+234.175)))
  return(svp)
}
ccinv<-function(svp){
  temp<-(234.175*(log(svp/6.1078)))/(17.08085-(log(svp/6.1078)))
  return(temp)
}
outlierfilter=TRUE #DAILY outlier filter for response vars

setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/H2.fin_Filtered_withRH/") #changed to include withRH
files<-list.files()
modelobjectsair<-list()
modelobjectsleaf<-list()

pdf(file="/Users/doris/Documents/UIowa/GAMS_Jun13/notnorm/SVP/plots_noOuts.pdf")
par(mfrow=c(2,3))

for(z in 1:length(files)){
  #for(z in 1:3){
  dat<-read.csv(files[z])
  print(paste0("SITE:      ",files[z],"        #",z))
  
  #Do the filters & complete cases
  datF<-dat[which(dat$windfilter=="keep" & dat$wetfilter=="keep" &
                    dat$growtimefilter=="keep" & dat$growseasfilter=="keep"),]
  if(z!=2){
  datF<-datF[,which(names(datF)%in%c("dat.Year","dat.DoY","VPDairmod","Tairmod",
                                     "TsurfEeqnmod" ,"VPDsurfEeqnmod",
                                     "GPPmod","NEEmod","LEmod","SMAP_rootzoneSM"))]
  datF<-datF[complete.cases(datF),]
  datF$svpAir<-cc(datF$Tairmod)
  datF$svpLeaf<-cc(datF$TsurfEeqnmod)
  datF$ones<-rep(1,length=nrow(datF))
  datD<-aggregate(datF[,3:12],by=list(datF$dat.DoY,datF$dat.Year),FUN=mean)
  datNums<-aggregate(datF$ones,by=list(datF$dat.DoY,datF$dat.Year),FUN=sum)
  datD$nRecord<-datNums$x
  
  
  } else {
  datF<-datF[,which(names(datF)%in%c("dat.Year","dat.DoY","VPDairmod","Tairmod",
                                       "GPPmod","NEEmod","LEmod","SMAP_rootzoneSM"))]
  datF<-datF[complete.cases(datF),]
  datF$svpAir<-cc(datF$Tairmod)
  datF$ones<-rep(1,length=nrow(datF))
  datD<-aggregate(datF[,3:9],by=list(datF$dat.DoY,datF$dat.Year),FUN=mean)
  datNums<-aggregate(datF$ones,by=list(datF$dat.DoY,datF$dat.Year),FUN=sum)
  datD$nRecord<-datNums$x
  }
  
  #datF$tcheck<-ccinv(datF$svp) #looks good
  
  if(outlierfilter==TRUE){
    
    datD$NEEmod[c(which(datD$NEEmod < (mean(datD$NEEmod,na.rm=T)-(3*sd(datD$NEEmod,na.rm=T)))),
                 which(datD$NEEmod > (mean(datD$NEEmod,na.rm=T)+(3*sd(datD$NEEmod,na.rm=T)))))]<-NA
    datD$GPPmod[c(which(datD$GPPmod < (mean(datD$GPPmod,na.rm=T)-(3*sd(datD$GPPmod,na.rm=T)))),
                 which(datD$GPPmod > (mean(datD$GPPmod,na.rm=T)+(3*sd(datD$GPPmod,na.rm=T)))))]<-NA 
    datD$LEmod[c(which(datD$LEmod < (mean(datD$LEmod,na.rm=T)-(3*sd(datD$LEmod,na.rm=T)))),
                which(datD$LEmod > (mean(datD$LEmod,na.rm=T)+(3*sd(datD$LEmod,na.rm=T)))))]<-NA 
    
    datD$svpAir[c(which(datD$svpAir < (mean(datD$svpAir,na.rm=T)-(3*sd(datD$svpAir,na.rm=T)))),
                  which(datD$svpAir > (mean(datD$svpAir,na.rm=T)+(3*sd(datD$svpAir,na.rm=T)))))]<-NA 
    datD$VPDairmod[c(which(datD$VPDairmod < (mean(datD$VPDairmod,na.rm=T)-(3*sd(datD$VPDairmod,na.rm=T)))),
                 which(datD$VPDairmod > (mean(datD$VPDairmod,na.rm=T)+(3*sd(datD$VPDairmod,na.rm=T)))))]<-NA 
    datD$SMAP_rootzoneSM[c(which(datD$SMAP_rootzoneSM < (mean(datD$SMAP_rootzoneSM,na.rm=T)-(3*sd(datD$SMAP_rootzoneSM,na.rm=T)))),
                which(datD$SMAP_rootzoneSM > (mean(datD$SMAP_rootzoneSM,na.rm=T)+(3*sd(datD$SMAP_rootzoneSM,na.rm=T)))))]<-NA 
    
    if(z!=2){
      datD$svpLeaf[c(which(datD$svpLeaf < (mean(datD$svpLeaf,na.rm=T)-(3*sd(datD$svpLeaf,na.rm=T)))),
                    which(datD$svpLeaf > (mean(datD$svpLeaf,na.rm=T)+(3*sd(datD$svpLeaf,na.rm=T)))))]<-NA 
      datD$VPDsurfEeqnmod[c(which(datD$VPDsurfEeqnmod < (mean(datD$VPDsurfEeqnmod,na.rm=T)-(3*sd(datD$VPDsurfEeqnmod,na.rm=T)))),
                    which(datD$VPDsurfEeqnmod > (mean(datD$VPDsurfEeqnmod,na.rm=T)+(3*sd(datD$VPDsurfEeqnmod,na.rm=T)))))]<-NA 
    }
    datD<-datD[complete.cases(datD),]
  }
  
  #MUST use the data argument here, or unexpected predict behavior!
  cutoffN<-3
  datD$color<-"black"
  datD$color[which(datD$nRecord<=cutoffN)]<-"grey"
  mdl1 = lm(svpAir~VPDairmod,data=datD)
  mdl1sub = lm(svpAir~VPDairmod,data=datD[which(datD$nRecord>cutoffN),])
  if(z!=2){
     mdl2 = lm(svpLeaf  ~VPDsurfEeqnmod,data=datD)
     mdl2sub = lm(svpLeaf~VPDsurfEeqnmod,data=datD[which(datD$nRecord>cutoffN),])
  }
  
  #PREDICTION AND PLOT FOR AIR
  VPDAirord<-data.frame(VPDairmod=datD$VPDairmod[order(datD$VPDairmod)])
  p1<-predict(mdl1,newdata=VPDAirord,type="response")
  p1sub<-predict(mdl1sub,newdata=VPDAirord,type="response")
  plot(datD$VPDairmod,datD$svpAir,pch=20,
       col=datD$color,main=substr(files[z],1,6))
    lines(VPDAirord$VPDairmod,p1,col="red",lwd=2)
    lines(VPDAirord$VPDairmod,p1sub,col="blue",lwd=2)
  legend("topleft",legend=c('n<3 days'),pch=20,col="grey",bty='n')
  legend("bottomright",c('all data','n>3 days'),
       lty=1, bty='n',col=c("red","blue"))  
    plot(datD$VPDairmod,residuals(mdl1))
    plot(datD[which(datD$nRecord>cutoffN),]$VPDairmod,residuals(mdl1sub))
    
  #PREDICTION AND PLOT FOR LEAF
    if(z!=2){
  VPDord<-data.frame(VPDsurfEeqnmod=datD$VPDsurfEeqnmod[order(datD$VPDsurfEeqnmod)])
  p2<-predict(mdl2,newdata=VPDord,type="response")
  p2sub<-predict(mdl2sub,newdata=VPDord,type="response")
    plot(datD$VPDsurfEeqnmod,datD$svpLeaf,pch=20,
         col=datD$color,main=substr(files[z],1,6))
    lines(VPDord$VPDsurfEeqnmod,p2,col="red",lwd=2)
    lines(VPDord$VPDsurfEeqnmod,p2sub,col="blue",lwd=2)
  legend("topleft",legend=c('n<3 days'),pch=20,col="grey",bty='n')
  legend("bottomright",c('all data','n>3 days'),
           lty=1, bty='n',col=c("red","blue"))  
    
  plot(datD$VPDsurfEeqnmod,residuals(mdl2))
  plot(datD[which(datD$nRecord>cutoffN),]$VPDsurfEeqnmod,residuals(mdl2sub))
    }
    if(z==2){ #To keep pages neat
      plot(0, xaxt = 'n', yaxt = 'n',  pch = '', ylab = '', xlab = '')
      plot(0, xaxt = 'n', yaxt = 'n', pch = '', ylab = '', xlab = '')
      plot(0, xaxt = 'n', yaxt = 'n', pch = '', ylab = '', xlab = '')  
      }

  #Back out temperature
  VPDAirtopredDF<-data.frame(VPDairmod=datD$VPDairmod)
  predSVPair<-predict(mdl1,newdata=VPDAirtopredDF,type="response")
  predTempair<-ccinv(predSVPair)
  datD$Tadj_svp_air<-datD$Tairmod-predTempair
  
  if(z!=2){
    VPDleaftopredDF<-data.frame(VPDsurfEeqnmod=datD$VPDsurfEeqnmod)
    predSVPleaf<-predict(mdl2,newdata=VPDleaftopredDF,type="response")
    predTempleaf<-ccinv(predSVPleaf)
    datD$Tadj_svp_leaf<-datD$TsurfEeqnmod-predTempleaf
  }
    
  modelobjectsair[[z]]<-mdl1
  if(z!=2){
  modelobjectsleaf[[z]]<-mdl2
  }
  write.csv(datD,paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/I2_TadjDF/SVP/",substr(files[z],1,6),"_D_SVP_Tadj.csv"))

}
dev.off()

saveRDS(modelobjectsair,"/Volumes/lss_ech2oLab/AMF_CoreSites/I2_TadjDF/SVP/SVP_air_Tadjmods_D.rds")
saveRDS(modelobjectsleaf,"/Volumes/lss_ech2oLab/AMF_CoreSites/I2_TadjDF/SVP/SVP_leaf_Tadjmods_D.rds")


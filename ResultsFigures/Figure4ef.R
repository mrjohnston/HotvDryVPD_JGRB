#This code can be sourced; the only thing to change is "choice"
choice<-"air"

setwd("/Users/doris/git/HotvDryVPD/H1aCSVs/")

HHairfiles<-list.files("/Users/doris/git/HotvDryVPD/H1aCSVs/",pattern="*AIR_HH.csv")
HHleaffiles<-list.files("/Users/doris/git/HotvDryVPD/H1aCSVs/",pattern="*LEAF_HH.csv")

nbincutoff<-5 #how many observations necessary to calculate range?

if(choice=="air"){
  HHf<-HHairfiles
  string="/Users/doris/Documents/UIowa/Figures/H1aFig_Feb2024_air.png"
  string2="/Users/doris/Documents/UIowa/Figures/H1aFig_Feb2024_air_vegcol.png"
  Xaxismax<-48
  #Xaxismax<-100
  }
if(choice=="leaf"){
  HHf<-HHleaffiles
  string="/Users/doris/Documents/UIowa/Figures/H1aFig_Feb2024_leaf.png"
  string2="/Users/doris/Documents/UIowa/Figures/H1aFig_Feb2024_leaf_vegcol.png"
  Xaxismax<-110
  }

sites<-substr(HHairfiles,1,6)

#Add colors to siteinfo
library(RColorBrewer)
siteinfo<-read.csv("/Volumes/lss_ech2oLab/AMF_CoreSites/FilteredSiteInfo_v6.csv")
siteinfo<-siteinfo[which(siteinfo$site%in%sites),]
 
 #COLORS BY MAP -- I fixed the scaling, but now IDK how to do colorbar.
 pal = colorRampPalette(c("azure3","cornflowerblue","blue4"))(100)
 siteinfo$MAPcol<-pal[cut(siteinfo$MAP, 100)]
 #Checking the colors -- looks good
 #plot(siteinfo$MAP,pch=20,col=siteinfo$MAPcol)

 #COLORS BY MAT
 pal = colorRampPalette(c("burlywood1","red","darkred"))(100)
 siteinfo$MATcol<-pal[cut(siteinfo$MAT, 100)]
 #Checking the colors -- looks good
 #plot(siteinfo$MAT,pch=20,col=siteinfo$MATcol)
 
 #COLORS BY VEG
 colpal = brewer.pal(length(unique(siteinfo$veg)), "Dark2")
 siteinfo$vegcol<-NA
 for(j in 1:length(unique(siteinfo$veg))){
   siteinfo$vegcol[which(siteinfo$veg==unique(siteinfo$veg)[j])]<-colpal[j]
   }

##########################################

par(xpd=FALSE)
png(filename=string,
    
    width=1200,height=1200,units="px",res=200)
zones <- matrix(c(1,3,4,
                  2,5,6,
                  7,7,0),
                ncol = 3, byrow = TRUE)
layout(zones, widths=c(0.1,1,.15), heights = c(0.8,0.8,0.1))

#The side labels:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Proportion of temperature range",cex=1.6,srt=90)

plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Proportion of RH range",cex=1.6,srt=90)

#par(mar=c(2,2,1.5,2))
par(mar=c(2,2,1.5,1))
for(i in 1:length (HHf)){
  f<-read.csv(HHf[i])
  names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
              "normrangeT_90","normrangerH_90","ninbin",
              "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
  f<-f[which(f$ninbin>nbincutoff),]
 
    lo<-loess(f$normrangeT_90~f$bin,span=0.5)
    pred<-predict(lo)
    if(i==1){
    plot(f$bin,pred, xaxt='n',
         col=siteinfo$MATcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],cex.axis=1.6,
         type="l", lwd=4,ylim=c(0,0.8),xlim=c(0,Xaxismax)) 
    axis(side=1,labels = F,tick = T)
    }
  if(i>1){
    lines(f$bin,pred,col=siteinfo$MATcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],lwd=4) 
  }
}
grid()

#MAT legend
par(mar=c(2,0.1,2,2)) #bottom, left, top, right
MATnums<-siteinfo$MAT
MAThex<-siteinfo$MATcol
MATnums2<-MATnums[order(MATnums)]
MAThex2<-MAThex[order(MATnums)]
legend_image <- rev(as.raster(matrix(MAThex2, ncol=1)))
plot(c(0,1.5),c(min(MATnums2),max(MATnums2)),
     type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.65,y=max(MATnums2)+0.7,expression(paste("MAT (",degree,"C)")),cex=1.2,xpd=NA)
text(x=1.5, y = seq(min(MATnums2),max(MATnums2),by=2), 
     labels = seq(min(MATnums2),max(MATnums2),by=2),cex=1.5,xpd=NA)
rasterImage(legend_image, xleft=-0.2, ybottom=min(MATnums2),
            xright=0.8, ytop=max(MATnums2),interpolate=T)

#MAP plot
#par(mar=c(2,2,1.5,2))
par(mar=c(2,2,1.5,1))
for(i in 1:length (HHf)){
  f<-read.csv(HHf[i])
  names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
              "normrangeT_90","normrangerH_90","ninbin",
              "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
  f<-f[which(f$ninbin>nbincutoff),]
  
  lo<-loess(f$normrangerH_90~f$bin,span=0.5)
  pred<-predict(lo)
  if(i==1){
    plot(f$bin,pred, xaxt='n',
         col=siteinfo$MAPcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],cex.axis=1.6,
         type="l", lwd=4,ylim=c(0,0.8),xlim=c(0,Xaxismax)) 
    axis(side=1,labels = T,tick = T,cex.axis=1.6)
  }
  if(i>1){
    lines(f$bin,pred,col=siteinfo$MAPcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],lwd=4) 
  }
}
grid()

#MAP legend
par(mar=c(2,0.1,2,2)) #bottom, left, top, right
MAPnums<-siteinfo$MAP
MAPhex<-siteinfo$MAPcol
MAPnums2<-MAPnums[order(MAPnums)]
MAPhex2<-MAPhex[order(MAPnums)]
legend_image <- rev(as.raster(matrix(MAPhex2, ncol=1)))
plot(c(0,1.5),c(min(MAPnums2),max(MAPnums2)),
     type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.65,y=max(MAPnums2)+40,expression(paste("MAP (mm)")),cex=1.2,xpd=NA)
text(x=1.5, y = seq(round(min(MAPnums2),0),max(MAPnums2),by=100), 
     labels = seq(round(min(MAPnums2),0),max(MAPnums2),by=100),cex=1.5,xpd=NA)
rasterImage(legend_image, xleft=-0.2, ybottom=min(MAPnums2),
            xright=0.8, ytop=max(MAPnums2),interpolate=T)

#The bottom label:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,-.2,"VPD (hPa)",cex=1.6)
dev.off()

######################################################
png(filename=string2,
    width=1200,height=1200,units="px",res=200)
zones <- matrix(c(1,3,5,
                  2,4,5,
                  6,6,0),
                ncol = 3, byrow = TRUE)
layout(zones, widths=c(0.1,1,.2), heights = c(0.8,0.8,0.1))

#The side labels:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Proportion of temperature range",cex=1.2,srt=90)

plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Proportion of RH range",cex=1.2,srt=90)

par(mar=c(2,2,1.5,2))
for(i in 1:length (HHf)){
  f<-read.csv(HHf[i])
  names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
              "normrangeT_90","normrangerH_90","ninbin",
              "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
  f<-f[which(f$ninbin>nbincutoff),]
  
  lo<-loess(f$normrangeT_90~f$bin,span=0.5)
  pred<-predict(lo)
  if(i==1){
    plot(f$bin,pred, xaxt='n',
         col=siteinfo$vegcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],cex.axis=1.6,
         type="l", lwd=4,ylim=c(0,0.8),xlim=c(0,Xaxismax)) 
    axis(side=1,labels = F,tick = T, cex.axis=1.6)
  }
  if(i>1){
    lines(f$bin,pred,col=siteinfo$vegcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],lwd=4) 
  }
}
grid()

#MAP plot
par(mar=c(2,2,1.5,2))
for(i in 1:length (HHf)){
  f<-read.csv(HHf[i])
  names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
              "normrangeT_90","normrangerH_90","ninbin",
              "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
  f<-f[which(f$ninbin>3),]
  
  lo<-loess(f$normrangerH_90~f$bin,span=0.5)
  pred<-predict(lo)
  if(i==1){
    plot(f$bin,pred, xaxt='n',
         col=siteinfo$vegcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],cex.axis=1.6,
         type="l", lwd=4,ylim=c(0,0.8),xlim=c(0,Xaxismax)) 
    axis(side=1,labels = T,tick = T,cex.axis=1.6)
  }
  if(i>1){
    lines(f$bin,pred,col=siteinfo$vegcol[which(siteinfo$site==substr(HHairfiles[i],1,6))],lwd=4) 
  }
}
grid()
#VEG legend
par(mar=c(2,0.5,2,2)) #bottom, left, top, right

veglabels<-unique(siteinfo$veg)
hex<-unique(siteinfo$vegcol)
#Double chekcing 
#unique(sub$vegcol[which(sub$veg=="WSA")])
legend_image <- rev(as.raster(matrix(hex, ncol=1)))
plot(c(0,1.5),c(1,7),
     type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = seq(1.1,7,by=0.9), 
     labels = veglabels,cex=0.8,xpd=NA)
rasterImage(legend_image, xleft=0, ybottom=0.75,
            xright=1, ytop=7,interpolate=F)

#The bottom label:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,-.2,"VPD (hPa)",cex=1.2)
dev.off()

#####################

##PLOTS WITH POINTS
# par(mfrow=c(4,4))
# for(i in 1:length (HHf)){
#   f<-read.csv(HHf[i])
#   names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
#               "normrangeT_90","normrangerH_90","ninbin",
#               "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
#   f<-f[which(f$ninbin>nbincutoff),]
#   
#   lo<-loess(f$normrangerH_90~f$bin,span=.5)
#   pred<-predict(lo)
#   
#   plot(f$bin,f$normrangerH_90,main=f$site[1])
#   lines(f$bin,pred,col="green",lwd=3)
# }

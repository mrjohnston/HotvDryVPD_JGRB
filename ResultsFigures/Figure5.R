#Check which models have significant interactions
library(mgcv)
library(scales)

airresults<-read.csv("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/airresults_ALL_modsummary_with3way.csv")
leafresults<-read.csv("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/leafresults_ALL_modsummary_with3way.csv")
info<-read.csv("/Volumes/lss_ech2oLab/AMF_CoreSites/FilteredSiteInfo_v5.csv")

infosub<-info[which(substr(info$site,4,6)%in%airresults$site),]

substr(infosub$site,4,6)==airresults$site #ok good
#substr(infosub$site,4,6)==leafresults$site #nope - bacause missing Ha

#Add growing season T and growing season VPD to infosub
infosub$GSVPD<-NA
infosub$GSairtemp<-NA
for (q in 1:nrow(infosub)){
  site<-substr(infosub$site,4,6)[q]
  dat<-read.csv(paste0("/Volumes/lss_ech2oLab/AMF_CoreSites/J_DailyFiles/US_",site,"_D_SVP_Tadj.csv"))
  infosub$GSVPD[q]<-mean(dat$VPDairmod)
  infosub$GSairtemp[q]<-mean(dat$Tairmod)
}
#Just checking
plot(infosub$MAT,infosub$GSairtemp); abline(0,1)

infosub$GppMod<-rep("black",nrow(infosub))
infosub$GppMod[which(airresults$GPPchoice==1)]<-"green"
infosub$GppMod[which(airresults$GPPchoice==3)]<-"green"

infosub$GppModDE<-airresults$GPPdevexpl
infosub$GppMod_pch<-NA
infosub$GppMod_pch[which(airresults$GPPchoice==1)]<-19
infosub$GppMod_pch[which(airresults$GPPchoice==3)]<-19
infosub$GppMod_pch[which(airresults$GPPchoice==2)]<-1

infosub$NeeMod<-rep("black",nrow(infosub))
infosub$NeeMod[which(airresults$NEEchoice==1)]<-"brown"
infosub$NeeMod[which(airresults$NEEchoice==3)]<-"brown"
infosub$NeeModDE<-airresults$NEEdevexpl
infosub$NeeMod_pch<-NA
infosub$NeeMod_pch[which(airresults$NEEchoice==1)]<-19
infosub$NeeMod_pch[which(airresults$NEEchoice==3)]<-19
infosub$NeeMod_pch[which(airresults$NEEchoice==2)]<-1

infosub$LeMod<-rep("black",nrow(infosub))
infosub$LeMod[which(airresults$LEchoice==1)]<-"blue"
infosub$LeMod[which(airresults$LEchoice==3)]<-"blue"
infosub$LeModDE<-airresults$LEdevexpl
infosub$LeMod_pch<-NA
infosub$LeMod_pch[which(airresults$LEchoice==1)]<-19
infosub$LeMod_pch[which(airresults$LEchoice==3)]<-19
infosub$LeMod_pch[which(airresults$LEchoice==2)]<-1

#size = consistently scaled
DE<-c(infosub$GppModDE,infosub$NeeModDE,infosub$LeModDE)
DEscale<-rescale(DE,to=c(0.5,5))
infosub$GppModScale<-DEscale[1:26]
infosub$NeeModScale<-DEscale[27:52]
infosub$LeModScale<-DEscale[53:length(DEscale)]

# png(filename="/Users/doris/Documents/UIowa/Figures/ModelDEInt_Summary2.png",
#     width=2100,height=750,units="px",res=220)
# zones <- matrix(c(1,3,4,5,0,
#                   0,2,2,2,0),
#                 ncol = 5, byrow = TRUE)
# layout(zones, widths=c(0.1,1,1,1,.8),
#        heights = c(1,.05))
# 
# #The side labels:
# par(mar = c(0,0,0,0))
# plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
# text(0,0,"Mean annual precipitation (mm)",cex=1.2,srt=90)
# 
# plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
# text(0,0,expression(paste("Mean annual temperature (",degree,"C)")),cex=1.2)
# 
# par(mar=c(2,2,2,1)) #bottom, left, top, right
# plot(infosub$MAT,infosub$MAP,col=infosub$GppMod,pch=infosub$GppMod_pch,cex=infosub$GppModScale,
#      xlab="",  ylab="",cex.lab=1.5,main="GPP",xlim=c(0,20),lwd=2)
# 
# plot(infosub$MAT,infosub$MAP,col=infosub$NeeMod,pch=infosub$NeeMod_pch,cex=infosub$NeeModScale,
#      xlab="",xlim=c(0,20),lwd=2,
#      ylab="",cex.lab=1.5,main="NEP",yaxt='n')
# axis(2,tick=T,labels = FALSE)
# plot(infosub$MAT,infosub$MAP,col=infosub$LeMod,pch=infosub$LeMod_pch,cex=infosub$LeModScale,
#      xlab="",yaxt='n',xlim=c(0,20),lwd=2,
#      ylab="",main="LE")
# axis(2,tick=T,labels = FALSE)
# 
# #This plot has to have the same scale? - easier to use xpd
# par(xpd=NA)
# points(c(24,  23,24,25, 24,24,24),
#        c(600, 700,700,700, 800,900,1000),
#        pch=c(1,16,16,16,16,16,16),lwd=2,
#        cex=c(median(DEscale),median(DEscale),median(DEscale),median(DEscale),
#              min(DEscale),median(DEscale),max(DEscale)),
#        col=c("black","green","brown","blue","black","black","black"))
# text(x = 26,y=1000,paste0("Deviance Expl. = ",round(max(DE)*100,0),"%"),cex=1.2,pos=4)
# text(x = 26,y=900, paste0("Deviance Expl. = ",round(median(DE)*100,0),"%"),cex=1.2,pos=4)
# text(x = 26,y=800, paste0("Deviance Expl. = ",round(min(DE)*100,2),"%"),cex=1.2,pos=4)
# text(x = 26,y=700, "Interaction",cex=1.2,pos=4)
# text(x = 26,y=600, "No Interaction",cex=1.2,pos=4)
# 
# dev.off()



############# IGBP

library(RColorBrewer)
colpal = brewer.pal(length(unique(infosub$veg)), "Dark2")
infosub$vegcol<-NA
infosub$vegcolR<-NA
infosub$vegcolG<-NA
infosub$vegcolB<-NA
for(j in 1:length(unique(infosub$veg))){
  infosub$vegcol[which(infosub$veg==unique(infosub$veg)[j])]<-colpal[j]
  }

png(filename="/Users/doris/Documents/UIowa/Figures/ModelDEInt_Summary_vegcol_AIR.png",
    width=2100,height=750,units="px",res=220)
zones <- matrix(c(1,3,4,5,0,
                  0,2,2,2,0),
                ncol = 5, byrow = TRUE)
layout(zones, widths=c(0.1,1,1,1,.8),
       heights = c(1,.05))

#The side labels:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Mean annual precipitation (mm)",cex=1.2,srt=90)

plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
#text(0,0,expression(paste("Mean annual temperature (",degree,"C)")),cex=1.2)
text(0,0,expression(paste("Mean growing season temperature (",degree,"C)")),cex=1.2)

par(mar=c(2,2,2,1)) #bottom, left, top, right
#plot(infosub$MAT,infosub$MAP,col=infosub$vegcol,pch=infosub$GppMod_pch,cex=infosub$GppModScale,
#     xlab="",  ylab="",cex.lab=1.5,main="GPP",xlim=c(0,20),lwd=2)
plot(infosub$GSairtemp,infosub$MAP,col=infosub$vegcol,pch=infosub$GppMod_pch,cex=infosub$GppModScale,
     xlab="",  ylab="",cex.lab=1.5,main="GPP",xlim=c(12,33),lwd=2)

#plot(infosub$MAT,infosub$MAP,col=infosub$vegcol,pch=infosub$NeeMod_pch,cex=infosub$NeeModScale,
#     xlab="",xlim=c(0,20),lwd=2,
#     ylab="",cex.lab=1.5,main="NEP",yaxt='n')
plot(infosub$GSairtemp,infosub$MAP,col=infosub$vegcol,pch=infosub$NeeMod_pch,cex=infosub$NeeModScale,
     xlab="",xlim=c(12,33),lwd=2,
     ylab="",cex.lab=1.5,main="NEP",yaxt='n')
axis(2,tick=T,labels = FALSE)

#plot(infosub$MAT,infosub$MAP,col=infosub$vegcol,pch=infosub$LeMod_pch,cex=infosub$LeModScale,
#     xlab="",yaxt='n',xlim=c(0,20),lwd=2,
#     ylab="",main="LE")
plot(infosub$GSairtemp,infosub$MAP,col=infosub$vegcol,pch=infosub$LeMod_pch,cex=infosub$LeModScale,
     xlab="",yaxt='n',xlim=c(12,33),lwd=2,
     ylab="",main="LE")
axis(2,tick=T,labels = FALSE)

#This plot has to have the same scale? - easier to use xpd
par(xpd=NA)
points(c(36,  35,36,37, 36,36,36),
       c(600, 700,700,700, 800,900,1000),
       pch=c(1,16,16,16,16,16,16),lwd=2,
       cex=c(median(DEscale),median(DEscale),median(DEscale),median(DEscale),
             min(DEscale),median(DEscale),max(DEscale)),
       col=c("black",unique(infosub$vegcol)[1],
             unique(infosub$vegcol)[2],unique(infosub$vegcol)[3],"black","black","black"))
# text(x = 26,y=1000,paste0("Deviance Expl. = ",round(max(DE)*100,0),"%"),cex=1.2,pos=4)
# text(x = 26,y=900, paste0("Deviance Expl. = ",round(median(DE)*100,0),"%"),cex=1.2,pos=4)
# text(x = 26,y=800, paste0("Deviance Expl. = ",round(min(DE)*100,2),"%"),cex=1.2,pos=4)
# text(x = 26,y=700, "Interaction",cex=1.2,pos=4)
# text(x = 26,y=600, "No Interaction",cex=1.2,pos=4)
text(x = 38,y=1000,paste0("Deviance Expl. = ",round(max(DE)*100,0),"%"),cex=1.2,pos=4)
text(x = 38,y=900, paste0("Deviance Expl. = ",round(median(DE)*100,0),"%"),cex=1.2,pos=4)
text(x = 38,y=800, paste0("Deviance Expl. = ",round(min(DE)*100,2),"%"),cex=1.2,pos=4)
text(x = 38,y=700, "Interaction",cex=1.2,pos=4)
text(x = 38,y=600, "No Interaction",cex=1.2,pos=4)

dev.off()

#VEG legend
png(filename="/Users/doris/Documents/UIowa/Figures/Veglegend.png",
    width=1000,height=140,units="px",res=300)
par(mar=c(0.5,0.5,0.5,0.5)) #bottom, left, top, right
veglabels<-unique(infosub$veg)
hex<-unique(infosub$vegcol)
legend_image <- (as.raster(matrix(hex, ncol=1))) #this used to be within "rev", why?
  #I think it's right now...
plot(#c(0,1.5),c(1,7),
  c(1,8),c(0,1.5),
     type = 'n', axes = F,xlab = '', ylab = '')
rasterImage(t(legend_image), xleft=1, ybottom=0,
            xright=8, ytop=1.7,interpolate=F)
text(#x=1.5, y = seq(1.1,7,by=0.9), 
  y=.75,x=seq(1.5,7.5,by=1),
  labels = veglabels,cex=0.8,xpd=NA)
dev.off()

# ##############
# png(filename="/Users/doris/Documents/UIowa/Figures/ModelDEInt_Summary.png",
#     width=2100,height=800,units="px",res=220)
# zones <- matrix(c(1,3,4,5,0,
#                   0,2,2,2,0),
#                 ncol = 5, byrow = TRUE)
# layout(zones, widths=c(0.1,1,1,1,.8),
#        heights = c(1,.05))
# 
# #The side labels:
# par(mar = c(0,0,0,0))
# plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
# text(0,0,"Mean annual precipitation (mm)",cex=1.2,srt=90)
# 
# plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
# text(0,0,expression(paste("Mean annual temperature (",degree,"C)")),cex=1.2)
# 
# par(mar=c(2,2,2,1)) #bottom, left, top, right
# plot(infosub$MAT,infosub$MAP,col=infosub$GppMod,pch=16,cex=infosub$GppModScale,
#      xlab="",  ylab="",cex.lab=1.5,main="GPP")
# 
# plot(infosub$MAT,infosub$MAP,col=infosub$NeeMod,pch=16,cex=infosub$NeeModScale,
#      xlab="",
#      ylab="",cex.lab=1.5,main="NEE",yaxt='n')
# axis(2,tick=T,labels = FALSE)
# plot(infosub$MAT,infosub$MAP,col=infosub$LeMod,pch=16,cex=infosub$LeModScale,
#      xlab="",yaxt='n',
#      ylab="",main="LE")
# axis(2,tick=T,labels = FALSE)
# 
# #This plot has to have the same scale? - easier to use xpd
# par(xpd=NA)
# points(c(22,  21,22,23, 22,22,22),
#        c(600, 700,700,700, 800,900,1000),
#        pch=16,
#        cex=c(median(DEscale),median(DEscale),median(DEscale),median(DEscale),
#              min(DEscale),median(DEscale),max(DEscale)),
#        col=c("black","green","brown","blue","black","black","black"))
# text(x = 24,y=1000,paste0("Deviance Expl. = ",round(max(DE)*100,0),"%"),cex=1.2,pos=4)
# text(x = 24,y=900, paste0("Deviance Expl. = ",round(median(DE)*100,0),"%"),cex=1.2,pos=4)
# text(x = 24,y=800, paste0("Deviance Expl. = ",round(min(DE)*100,2),"%"),cex=1.2,pos=4)
# text(x = 24,y=700, "Interaction",cex=1.2,pos=4)
# text(x = 24,y=600, "No Interaction",cex=1.2,pos=4)
# 
# dev.off()
# 
# par(mar = c(0,0,0,0))
# plot(c(0,-0.3,0,0.3,0,0,0),c(-3,-2,-2,-2,-1,0,1),pch=16,col=c("black","green","brown","blue","black","black","black"),
#      ylim=c(-6,3),xlim=c(0,4),cex=c(3.2,3.2,3.2,3.2,.5,3.2,5),axes=FALSE,
#      xlab="",ylab="")
# text(x = 0.5,y=-1, "Dev. Expl. = 8.6%",cex=1.5,pos=4)
# text(x = 0.5,y=0, "Dev. Expl. = 55%",cex=1.5,pos=4)
# text(x = 0.5,y=1, "Dev. Expl. = 85%",cex=1.5,pos=4)
# text(x = 0.5,y=-2, "Interaction",cex=1.5,pos=4)
# text(x = 0.5,y=-3, "No Interaction",cex=1.5,pos=4)

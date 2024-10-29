#This code can be sourced; the only thing to change is "choice"
choice<-"air"

library(viridis)
library(plyr)

setwd("/Users/doris/git/HotvDryVPD/H1aCSVs/")

HHairfiles<-list.files("/Users/doris/git/HotvDryVPD/H1aCSVs/",pattern="*AIR_HH.csv")
HHleaffiles<-list.files("/Users/doris/git/HotvDryVPD/H1aCSVs/",pattern="*LEAF_HH.csv")

nbincutoff<-5 #how many observations necessary to calculate range?

if(choice=="air"){
  HHf<-HHairfiles
  filename1<-"/Users/doris/Documents/UIowa/Figures/H1aSiteRangeAll_Feb2024_air.png"
  filename2<-"/Users/doris/Documents/UIowa/Figures/H1aSiteRangeSubset_Feb2024_air.png"
}
if(choice=="leaf"){
  HHf<-HHleaffiles
  filename1<-"/Users/doris/Documents/UIowa/Figures/H1aSiteRangeAll_Feb2024_leaf.png"
  filename2<-"/Users/doris/Documents/UIowa/Figures/H1aSiteRangeSubset_Feb2024_leaf.png"
}

sites<-substr(HHairfiles,1,6)

####################
#Colors for VPD
interval=0.1 #VPD interval, in hPa

#Read in all the files and stick them together to get the correct, consistent coloring
for(i in 1:length(HHf)){
  f<-read.csv(HHf[i])
  names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
              "normrangeT_90","normrangerH_90","ninbin",
              "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
  f<-f[which(f$ninbin>nbincutoff),]
 
  if(i==1){
    all<-f
  }
  if(i>1){
    all<-rbind(all,f)
  }
}

infcolors<-inferno(length(unique(all$bin)))
#Trying out some other color schemes for Joel
#infcolors<-viridis(length(unique(all$bin)))
#infcolors<-magma(length(unique(all$bin)))

#https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal
  #^Get around these very annoying tolerance issues by mapping factors instead:
colormatchdf<-data.frame(unique(all$bin)[order(unique(all$bin))],infcolors)
names(colormatchdf)<-c("bin","inf")
#plot(colormatchdf[,1],col=colormatchdf[,2]) #good
colormatchdf$bin<-as.factor(colormatchdf$bin)
toreplace<-as.factor(all$bin)
all$inf<-mapvalues(toreplace, from=colormatchdf$bin, to=colormatchdf$inf)

#########################################################
#Plot with all sites
png(filename=filename1,
     width=1600,height=1000,units="px",res=200)
 zones <- matrix(c(1,2,3,4,5,6,7,8,
                   1,9,10,11,12,13,14,15,
                   1,16,17,18,19,20,21,22,
                   1,23,24,25,26,27,28,0,
                   0,29,29,29,29,29,29,29),
                 ncol = 8, byrow = TRUE)
 layout(zones, widths=c(.2,1,1,1,1,1,1,1), heights = c(1,1,1,1,.2))
#The side label:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Relative humidity (%)",cex=1.2,srt=90)

par(mar = c(1,1.5,1.8,.3),mgp=c(3,0.5,0)) 
  #mar: Bottom, left, top, right
  #mgp: axis title, axis labels, axis line [default: (3,1,0)]

#Produces a plot for each site:
for(i in 1:length(sites)){
  
  #site & data subset gets its own panel:
  d<-all[which(all$site==sites[i]),]
  d$inf<-as.character(d$inf)
  
  #empty plot:
  if(i==2 & choice=="leaf"){
    plot(1, type="n", xlab="",ylab="",main="Ha1",xaxt='n',yaxt='n')
    grid()
    next
  } 
  
  xlimit=c(min(d$sq05T)-2,max(d$sq95T)+2)
  ylimit=c(min(d$sq05rH)-5,max(d$sq95rH)+5)
  plot(1, type="n", xlab="",ylab="",main=substr(unique(d$site),4,6),#xaxt='n',yaxt='n',
       xlim=xlimit,ylim=ylimit)

  #Add the data - median points & ranges
  points(d$sq50T,d$sq50rH,col=d$inf,cex=0.2)
    segments(x0=d$sq50T,y0=d$sq05rH,x1=d$sq50T,y1=d$sq95rH,col=d$inf)
    segments(x0=d$sq05T,y0=d$sq50rH,x1=d$sq95T,y1=d$sq50rH,col=d$inf)

  grid()
}

#VPD legend
VPDmeds<-unique(all$bin)[order(unique(all$bin))]
par(mar = c(1,2,2,2)) #Bottom, left, top, right
legend_image <- rev(as.raster(matrix(infcolors, ncol=1)))
plot(c(0,1.5),c(min(VPDmeds),max(VPDmeds)),
     type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.75,y=max(VPDmeds)+4,"VPD (hPa)",cex=1,xpd=NA)
text(x=1.5, y = seq(min(VPDmeds),max(VPDmeds),by=interval*100), 
     labels = seq(min(VPDmeds),max(VPDmeds),by=interval*100),cex=0.8,xpd=NA)
rasterImage(legend_image, xleft=0, ybottom=min(VPDmeds),
            xright=1, ytop=max(VPDmeds))

#The bottom label:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
if(choice=="air"){
  text(0,-.2,expression(paste("Air temperature (",degree,"C)")),cex=1.2)
}
if(choice=="leaf"){
  text(0,-.2,expression(paste("Surface temperature (",degree,"C)")),cex=1.2)
}
dev.off()

#############################################################
#############################################################
#Only Wkg, Ton, UMB, Ho1 -- in this order

#Colors for VPD
#Read in all the files and stick them together to get the correct, consistent coloring
for(i in c(grep("Wkg",HHf),grep("Ton",HHf),grep("UMB",HHf),grep("Ho1",HHf))){
  f<-read.csv(HHf[i])
  names(f)<-c("site","VPDmed","bin","srangeT_90","srangerH_90",
              "normrangeT_90","normrangerH_90","ninbin",
              "sq05T","sq50T","sq95T","sq05rH","sq50rH","sq95rH")
  f<-f[which(f$ninbin>nbincutoff),]
  
  if(i==grep("Wkg",HHf)){
    all2<-f
  } else {
    all2<-rbind(all2,f)
  }
}

infcolors<-inferno(length(unique(all2$bin)))
#infcolors<-viridis(length(unique(all2$bin)))
#infcolors<-magma(length(unique(all2$bin)))
#https://stackoverflow.com/questions/9508518/why-are-these-numbers-not-equal
#^Get around these very annoying tolerance issues by mapping factors instead:
colormatchdf<-data.frame(unique(all2$bin)[order(unique(all2$bin))],infcolors)
names(colormatchdf)<-c("bin","inf")
#plot(colormatchdf[,1],col=colormatchdf[,2]) #good
colormatchdf$bin<-as.factor(colormatchdf$bin)
toreplace<-as.factor(all2$bin)
all2$inf<-mapvalues(toreplace, from=colormatchdf$bin, to=colormatchdf$inf)

png(filename=filename2,
    width=1200,height=1200,units="px",res=200)
zones <- matrix(c(1,2,3,6,
                  1,4,5,6,
                  0,7,7,0),
                ncol = 4, byrow = TRUE)
layout(zones, widths=c(.2,1,1,.3), heights = c(1,1,.2))

#The side label:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')
text(0,0,"Relative humidity (%)",cex=1.6,srt=90)

par(mar = c(1,1.5,1.8,.3),mgp=c(3,0.5,0)) 
#mar: Bottom, left, top, right
#mgp: axis title, axis labels, axis line [default: (3,1,0)]


#Produces a plot for sites of interest:
sites<-c("US_Wkg","US_Ton","US_UMB","US_Ho1")
for(i in 1:4){
  
  #site & data subset gets its own panel:
  d<-all2[which(all2$site==sites[i]),]
  d$inf<-as.character(d$inf)
  
  #empty plot:
  if(i==1){
    plot(1, type="n", xlab="",ylab="",main=substr(unique(d$site),4,6),xaxt='n',
         xlim=c(min(all2$sq05T),max(all2$sq95T)),
         ylim=c(min(all2$sq05rH),max(all2$sq95rH)),cex.axis=1.6)
    axis(side=1,labels = F,tick = T)
  }
  if(i==2){
    plot(1, type="n", xlab="",ylab="",main=substr(unique(d$site),4,6),xaxt='n',yaxt='n',
         xlim=c(min(all2$sq05T),max(all2$sq95T)),
         ylim=c(min(all2$sq05rH),max(all2$sq95rH)),cex.axis=1.6)
    axis(side=1,labels = F,tick = T)
    axis(side=2,labels = F,tick = T)
  }
  if(i==3){
    plot(1, type="n", xlab="",ylab="",main=substr(unique(d$site),4,6),
         xlim=c(min(all2$sq05T),max(all2$sq95T)),
         ylim=c(min(all2$sq05rH),max(all2$sq95rH)),cex.axis=1.6)
  }
  if(i==4){
    plot(1, type="n", xlab="",ylab="",main=substr(unique(d$site),4,6),yaxt='n',
         xlim=c(min(all2$sq05T),max(all2$sq95T)),
         ylim=c(min(all2$sq05rH),max(all2$sq95rH)),cex.axis=1.6)
    axis(side=2,labels = F,tick = T)
  }
  
  #Add the data - median points & ranges
  points(d$sq50T,d$sq50rH,col=d$inf,cex=0.2)
    segments(x0=d$sq50T,y0=d$sq05rH,x1=d$sq50T,y1=d$sq95rH,col=d$inf)
    segments(x0=d$sq05T,y0=d$sq50rH,x1=d$sq95T,y1=d$sq50rH,col=d$inf)

  grid()
}

#VPD legend
VPDmeds<-unique(all2$bin)[order(unique(all2$bin))]
par(mar = c(1,.5,2,1)) #Bottom, left, top, right
legend_image <- rev(as.raster(matrix(infcolors, ncol=1)))
plot(c(0,1.5),c(0,max(VPDmeds)),
     type = 'n', axes = F,xlab = '', ylab = '')
text(x=0.75,y=max(VPDmeds)+2,"VPD (hPa)",cex=1.2,xpd=NA)
text(x=1.5, y = seq(0,max(VPDmeds),by=interval*100), 
     labels = seq(0,max(VPDmeds),by=interval*100),cex=1.5,xpd=NA)
rasterImage(legend_image, xleft=0, ybottom=0,
            xright=1, ytop=max(VPDmeds))

#The bottom label:
par(mar = c(0,0,0,0))
plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1),yaxt="n",xaxt="n",bty='n')

if(choice=="air"){
  text(0,-.2,expression(paste("Air temperature (",degree,"C)")),cex=1.6)
}
if(choice=="leaf"){
  text(0,-.2,expression(paste("Surface temperature (",degree,"C)")),cex=1.6)
}

dev.off()

####################################################


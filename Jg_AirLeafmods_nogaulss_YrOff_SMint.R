#SAME AS Jb, (no gaulss), but year is offset ONLY
#AND: include SM interaction

library(mgcv)
gam.control(maxit=300) #default maxit=200 - this didn't work
#300 didn't converge
#400 yielded Error in t(coef1) %*% St : 436 arguments passed to '%*%' which requires 2
options(warn=2) #break loop for warning (default is options(warn=0))
#I specified k = ... but NOT in the gaulss term - does that matter? 
#check the models and see how they look.
#why specify k to deal with lack of convergence? 
#https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.convergence.html
kchoice=12
kchoice2=5 #for the interaction, will effectively be squared
kchoice3 = 3 # for the three-way interaction

setwd("/Volumes/lss_ech2oLab/AMF_CoreSites/J_DailyFiles/")
files<-list.files(pattern="*Tadj.csv*")

for(z in 1:length(files)){ 
  dat<-read.csv(files[z])
  dat$dat.Year<-as.factor(dat$dat.Year) #UGH, why didn't this work when I tried to save as part of J0??
  dat$NEEmod<-dat$NEEmod*-1
  print(paste0("SITE:      ",files[z],"        #",z))
  print(dim(dat))

  # ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("running GPP mods")
  #no smooths
  AirGPP0<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + VPDairmod + Tadj_svp_air,data=dat,method="REML")
  #one smooth
  AirGPP1<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDairmod,k=kchoice) + Tadj_svp_air,data=dat,method="REML")
  AirGPP2<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + VPDairmod + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  AirGPPa<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDairmod + Tadj_svp_air,data=dat,method="REML")
  #two smooths
  AirGPPb<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + Tadj_svp_air, data=dat,method="REML")
  AirGPPc<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDairmod + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  AirGPP3<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  #all smooths
  AirGPPd<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  #2-way interactions: VPD & Tadj  / SM & VPD  / SM & Tadj  
  AirGPPe<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(VPDairmod,Tadj_svp_air,k=kchoice2),data=dat,method="REML")
  AirGPP4<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(VPDairmod,SMAP_rootzoneSM,k=kchoice2),data=dat,method="REML")
  AirGPP5<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(SMAP_rootzoneSM,Tadj_svp_air,k=kchoice2),data=dat,method="REML")
  #combinations of 2-way interactions
  AirGPP6<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(VPDairmod,Tadj_svp_air,k=kchoice2)+
                 ti(VPDairmod,SMAP_rootzoneSM,k=kchoice2),
                 data=dat,method="REML")
  AirGPP7<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(SMAP_rootzoneSM,Tadj_svp_air,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2),
                 data=dat,method="REML")  
  AirGPP8<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                ti(Tadj_svp_air,VPDairmod,k=kchoice2),
                data=dat,method="REML")
  AirGPP9<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2),
               data=dat,method="REML")
   #3-way interaction
  AirGPPf<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2)+
                 ti(VPDairmod,Tadj_svp_air,SMAP_rootzoneSM,k=kchoice3),data=dat,method="REML")
  
  AirGPPmods<-list(AirGPPa,AirGPPb,AirGPPc,AirGPPd,AirGPPe,AirGPPf,AirGPP0,AirGPP1,AirGPP2,AirGPP3,AirGPP4,AirGPP5,AirGPP6,AirGPP7,AirGPP8,AirGPP9)
  saveRDS(AirGPPmods,paste0("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/GPP_D_Air_",substr(files[z],4,6),".rds"))
  
  if(z!=2){
    LeafGPP0<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + VPDsurfEeqnmod + Tadj_svp_leaf,data=dat,method="REML")
    #one smooth
    LeafGPP1<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + Tadj_svp_leaf,data=dat,method="REML")
    LeafGPP2<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + VPDsurfEeqnmod + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    LeafGPPa<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDsurfEeqnmod + Tadj_svp_leaf,data=dat,method="REML")
    #two smooths
    LeafGPPb<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + Tadj_svp_leaf, data=dat,method="REML")
    LeafGPPc<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDsurfEeqnmod + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    LeafGPP3<-gam(GPPmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    #all smooths
    LeafGPPd<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    #2-way interactions: VPD & Tadj  / SM & VPD  / SM & Tadj 
    LeafGPPe<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(VPDsurfEeqnmod,Tadj_svp_leaf,k=kchoice2),data=dat,method="REML")
    LeafGPP4<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(VPDsurfEeqnmod,SMAP_rootzoneSM,k=kchoice2),data=dat,method="REML")
    LeafGPP5<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(SMAP_rootzoneSM,Tadj_svp_leaf,k=kchoice2),data=dat,method="REML")
    
    #combinations of 2-way interactions
    LeafGPP6<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                   ti(VPDsurfEeqnmod,Tadj_svp_leaf,k=kchoice2)+
                   ti(VPDsurfEeqnmod,SMAP_rootzoneSM,k=kchoice2),
                 data=dat,method="REML")
    LeafGPP7<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                   ti(SMAP_rootzoneSM,Tadj_svp_leaf,k=kchoice2)+
                   ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2),
                 data=dat,method="REML")  
    LeafGPP8<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                   ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                   ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2),
                 data=dat,method="REML")
    LeafGPP9<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                   ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                   ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2)+
                   ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2),
                 data=dat,method="REML")
    #3-way interaction
    LeafGPPf<-gam(GPPmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                    ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2)+
                    ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2)+
                    ti(VPDsurfEeqnmod,Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice3),data=dat,method="REML")
    
    LeafGPPmods<-list(LeafGPPa,LeafGPPb,LeafGPPc,LeafGPPd,LeafGPPe,LeafGPPf,LeafGPP0,LeafGPP1,LeafGPP2,LeafGPP3,LeafGPP4,LeafGPP5,LeafGPP6,LeafGPP7,LeafGPP8,LeafGPP9)
    saveRDS(LeafGPPmods,paste0("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/GPP_D_Leaf_",substr(files[z],4,6),".rds"))
  }
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("running NEE mods")
  #no smooths
  AirNEE0<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + VPDairmod + Tadj_svp_air,data=dat,method="REML")
  #one smooth
  AirNEE1<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDairmod,k=kchoice) + Tadj_svp_air,data=dat,method="REML")
  AirNEE2<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + VPDairmod + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  AirNEEa<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDairmod + Tadj_svp_air,data=dat,method="REML")
  #two smooths
  AirNEEb<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + Tadj_svp_air, data=dat,method="REML")
  AirNEEc<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDairmod + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  AirNEE3<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  #all smooths
  AirNEEd<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  #2-way interactions: VPD & Tadj  / SM & VPD  / SM & Tadj 
  AirNEEe<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(VPDairmod,Tadj_svp_air,k=kchoice2),data=dat,method="REML")
  AirNEE4<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(VPDairmod,SMAP_rootzoneSM,k=kchoice2),data=dat,method="REML")
  AirNEE5<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(SMAP_rootzoneSM,Tadj_svp_air,k=kchoice2),data=dat,method="REML")
  #combinations of 2-way interactions
  AirNEE6<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(VPDairmod,Tadj_svp_air,k=kchoice2)+
                 ti(VPDairmod,SMAP_rootzoneSM,k=kchoice2),
               data=dat,method="REML")
  AirNEE7<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(SMAP_rootzoneSM,Tadj_svp_air,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2),
               data=dat,method="REML")  
  AirNEE8<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2),
               data=dat,method="REML")
  AirNEE9<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2),
               data=dat,method="REML")
  #3-way interaction
  AirNEEf<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2)+
                 ti(VPDairmod,Tadj_svp_air,SMAP_rootzoneSM,k=kchoice3),data=dat,method="REML")
  
  
   AirNEEmods<-list(AirNEEa,AirNEEb,AirNEEc,AirNEEd,AirNEEe,AirNEEf,AirNEE0,AirNEE1,AirNEE2,AirNEE3,AirNEE4,AirNEE5,AirNEE6,AirNEE7,AirNEE8,AirNEE9)
  saveRDS(AirNEEmods,paste0("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/NEE_D_Air_",substr(files[z],4,6),".rds"))

  if(z!=2){
    LeafNEE0<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + VPDsurfEeqnmod + Tadj_svp_leaf,data=dat,method="REML")
    #one smooth
    LeafNEE1<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + Tadj_svp_leaf,data=dat,method="REML")
    LeafNEE2<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + VPDsurfEeqnmod + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    LeafNEEa<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDsurfEeqnmod + Tadj_svp_leaf,data=dat,method="REML")
    #two smooths
    LeafNEEb<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + Tadj_svp_leaf, data=dat,method="REML")
    LeafNEEc<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDsurfEeqnmod + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    LeafNEE3<-gam(NEEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    #all smooths
    LeafNEEd<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    #2-way interactions: VPD & Tadj  / SM & VPD  / SM & Tadj  & ther combos
    LeafNEEe<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(VPDsurfEeqnmod,Tadj_svp_leaf,k=kchoice2),data=dat,method="REML")
    LeafNEE4<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(VPDsurfEeqnmod,SMAP_rootzoneSM,k=kchoice2),data=dat,method="REML")
    LeafNEE5<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(SMAP_rootzoneSM,Tadj_svp_leaf,k=kchoice2),data=dat,method="REML")
    #combinations of 2-way interactions
    LeafNEE6<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(VPDsurfEeqnmod,Tadj_svp_leaf,k=kchoice2)+
                    ti(VPDsurfEeqnmod,SMAP_rootzoneSM,k=kchoice2),
                  data=dat,method="REML")
    LeafNEE7<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(SMAP_rootzoneSM,Tadj_svp_leaf,k=kchoice2)+
                    ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2),
                  data=dat,method="REML")  
    LeafNEE8<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                    ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2),
                  data=dat,method="REML")
    LeafNEE9<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                    ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2)+
                    ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2),
                  data=dat,method="REML")
    
    #3-way interaction
    LeafNEEf<-gam(NEEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                    ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2)+
                    ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2)+
                    ti(VPDsurfEeqnmod,Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice3),data=dat,method="REML")
    
    LeafNEEmods<-list(LeafNEEa,LeafNEEb,LeafNEEc,LeafNEEd,LeafNEEe,LeafNEEf,LeafNEE0,LeafNEE1,LeafNEE2,LeafNEE3,LeafNEE4,LeafNEE5,LeafNEE6,LeafNEE7,LeafNEE8,LeafNEE9)
    saveRDS(LeafNEEmods,paste0("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/NEE_D_Leaf_",substr(files[z],4,6),".rds"))
    } 

  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  print("running LE mods")
  #no smooths
  AirLE0<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + VPDairmod + Tadj_svp_air,data=dat,method="REML")
  #one smooth
  AirLE1<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDairmod,k=kchoice) + Tadj_svp_air,data=dat,method="REML")
  AirLE2<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + VPDairmod + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  AirLEa<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDairmod + Tadj_svp_air,data=dat,method="REML")
  #two smooths
  AirLEb<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + Tadj_svp_air, data=dat,method="REML")
  AirLEc<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDairmod + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  AirLE3<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  #all smooths
  AirLEd<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice),data=dat,method="REML")
  #2-way interactions: VPD & Tadj  / SM & VPD  / SM & Tadj 
  AirLEe<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(VPDairmod,Tadj_svp_air,k=kchoice2),data=dat,method="REML")
  AirLE4<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(VPDairmod,SMAP_rootzoneSM,k=kchoice2),data=dat,method="REML")
  AirLE5<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) + ti(SMAP_rootzoneSM,Tadj_svp_air,k=kchoice2),data=dat,method="REML")
  
  #combinations of 2-way interactions
  AirLE6<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(VPDairmod,Tadj_svp_air,k=kchoice2)+
                 ti(VPDairmod,SMAP_rootzoneSM,k=kchoice2),
               data=dat,method="REML")
  AirLE7<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(SMAP_rootzoneSM,Tadj_svp_air,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2),
               data=dat,method="REML")  
  AirLE8<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2),
               data=dat,method="REML")
  AirLE9<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                 ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                 ti(Tadj_svp_air,VPDairmod,k=kchoice2)+
                 ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2),
               data=dat,method="REML")
  #3-way interaction
  AirLEf<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDairmod,k=kchoice) + ti(Tadj_svp_air,k=kchoice) +
                ti(Tadj_svp_air,SMAP_rootzoneSM,k=kchoice2)+
                ti(Tadj_svp_air,VPDairmod,k=kchoice2)+
                ti(SMAP_rootzoneSM,VPDairmod,k=kchoice2)+
                ti(VPDairmod,Tadj_svp_air,SMAP_rootzoneSM,k=kchoice3),data=dat,method="REML")
  
  
  AirLEmods<-list(AirLEa,AirLEb,AirLEc,AirLEd,AirLEe,AirLEf,AirLE0,AirLE1,AirLE2,AirLE3,AirLE4,AirLE5,AirLE6,AirLE7,AirLE8,AirLE9)
    saveRDS(AirLEmods,paste0("/Users/doris/Documents/UIowa/GAMs_Jul7/Air_moreMods/LE_D_Air_",substr(files[z],4,6),".rds"))

  if(z!=2){
    LeafLE0<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + VPDsurfEeqnmod + Tadj_svp_leaf,data=dat,method="REML")
    #one smooth
    LeafLE1<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + Tadj_svp_leaf,data=dat,method="REML")
    LeafLE2<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + VPDsurfEeqnmod + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    LeafLEa<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDsurfEeqnmod + Tadj_svp_leaf,data=dat,method="REML")
    #two smooths
    LeafLEb<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + Tadj_svp_leaf, data=dat,method="REML")
    LeafLEc<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + VPDsurfEeqnmod + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    LeafLE3<-gam(LEmod ~ SMAP_rootzoneSM + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    #all smooths
    LeafLEd<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice),data=dat,method="REML")
    #2-way interactions: VPD & Tadj  / SM & VPD  / SM & Tadj 
    LeafLEe<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(VPDsurfEeqnmod,Tadj_svp_leaf,k=kchoice2),data=dat,method="REML")
    LeafLE4<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(VPDsurfEeqnmod,SMAP_rootzoneSM,k=kchoice2),data=dat,method="REML")
    LeafLE5<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) + ti(SMAP_rootzoneSM,Tadj_svp_leaf,k=kchoice2),data=dat,method="REML")
    #combinations of 2-way interactions
    LeafLE6<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(VPDsurfEeqnmod,Tadj_svp_leaf,k=kchoice2)+
                    ti(VPDsurfEeqnmod,SMAP_rootzoneSM,k=kchoice2),
                  data=dat,method="REML")
    LeafLE7<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(SMAP_rootzoneSM,Tadj_svp_leaf,k=kchoice2)+
                    ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2),
                  data=dat,method="REML")  
    LeafLE8<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                    ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2),
                  data=dat,method="REML")
    LeafLE9<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                    ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                    ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2)+
                    ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2),
                  data=dat,method="REML")
    #3-way interaction
    LeafLEf<-gam(LEmod ~ ti(SMAP_rootzoneSM,k=kchoice) + dat.Year + ti(VPDsurfEeqnmod,k=kchoice) + ti(Tadj_svp_leaf,k=kchoice) +
                   ti(Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice2)+
                   ti(Tadj_svp_leaf,VPDsurfEeqnmod,k=kchoice2)+
                   ti(SMAP_rootzoneSM,VPDsurfEeqnmod,k=kchoice2)+
                   ti(VPDsurfEeqnmod,Tadj_svp_leaf,SMAP_rootzoneSM,k=kchoice3),data=dat,method="REML")
    
    LeafLEmods<-list(LeafLEa,LeafLEb,LeafLEc,LeafLEd,LeafLEe,LeafLEf,LeafLE0,LeafLE1,LeafLE2,LeafLE3,LeafLE4,LeafLE5,LeafLE6,LeafLE7,LeafLE8,LeafLE9)

    saveRDS(LeafLEmods,paste0("/Users/doris/Documents/UIowa/GAMs_Jul7/Leaf_moreMods/LE_D_Leaf_",substr(files[z],4,6),".rds"))
    }
    
gc()
}






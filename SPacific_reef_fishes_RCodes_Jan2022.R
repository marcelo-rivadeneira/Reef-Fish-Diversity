rm(list = ls())
library(vegan)
library(iNEXT)
library(nlme)
library(ape)
library(FD)
library(MuMIn)
library(plotrix)
library(beanplot)
library(plyr)
library(picante)
library(vegan)
library(tidyverse)
library(RVAideMemoire)
setwd("~/Dropbox/Proyectos/Current Projects/Reef Fishes/Script Jan2022")
census=read.csv("SPacific_reef_fishes_census.csv",T,sep=";")
traits=read.csv("SPacific_reef_fishes_traits.csv",T,sep=";",row.names = 1,
                colClasses=c("factor","factor","factor","factor","numeric",
                             "numeric","factor","factor","factor"))
coord=read.csv("SPacific_reef_fishes_coordsites.csv",T,sep=";")
names(census)[1:10]
names(traits)
str(traits)

agg.Location=aggregate(census[,9:141],by=list(census$Location),sum);agg.Location$Scale=rep("Location",nrow(agg.Location));agg.Location$Site.name=agg.Location$Group.1; agg.Location$Realm=census$Realm[match(agg.Location$Site.name,census$Location)]; agg.Location$Province=census$Province[match(agg.Location$Site.name,census$Location)]
agg.Site=aggregate(census[,9:141],by=list(census$Site),sum);agg.Site$Scale=rep("Site",nrow(agg.Site));agg.Site$Site.name=agg.Site$Group.1; agg.Site$Realm=census$Realm[match(agg.Site$Site.name,census$Site)];agg.Site$Province=census$Province[match(agg.Site$Site.name,census$Site)]
agg.Transect=aggregate(census[,9:141],by=list(census$Transect),sum);agg.Transect$Scale=rep("Transect",nrow(agg.Transect));agg.Transect$Site.name=census$Site[match(agg.Transect$Group.1,census$Transect)];agg.Transect$Realm=census$Realm[match(agg.Transect$Site.name,census$Site)];agg.Transect$Province=census$Province[match(agg.Transect$Site.name,census$Site)]
agg.scales=rbind(agg.Site,agg.Transect,agg.Location)
names(agg.scales)


## Calculating Chao 
table.diversity=data.frame(t(estimateR(agg.scales[,2:134])))

## Coverage-based rarefaction (using iNEXT)
## https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
agg2.Location=t(agg.Location[,2:134]); colnames(agg2.Location)=agg.Location$Site.name
agg2.Site=t(agg.Site[,2:134]); colnames(agg2.Site)=agg.Site$Site.name
agg2.Transect=t(agg.Transect[,2:134]); colnames(agg2.Transect)=agg.Transect$Site.name

cov.Location=iNEXT(agg2.Location, q=0, datatype="abundance", endpoint=10000)  
cov.Site=iNEXT(agg2.Site, q=0, datatype="abundance", endpoint=2500)
cov.Transect=iNEXT(agg2.Transect, q=0, datatype="abundance", endpoint=1000)

# coverage
mean(cov.Location$DataInfo$SC); sd(cov.Location$DataInfo$SC); quantile(cov.Location$DataInfo$SC,0.05)
mean(cov.Site$DataInfo$SC); sd(cov.Site$DataInfo$SC);quantile(cov.Site$DataInfo$SC,0.05)
mean(cov.Transect$DataInfo$SC); sd(cov.Transect$DataInfo$SC);quantile(cov.Transect$DataInfo$SC,0.05)

chao.Location=ChaoSpecies(agg2.Location, datatype="abundance")
chao.Site=ChaoSpecies(agg2.Site, datatype="abundance")
chao.Transect=ChaoSpecies(agg2.Transect, datatype="abundance")

chao.all=rbind(chao.Site,chao.Transect,chao.Location)

cor.test(log(chao.Location$Observed),log(chao.Location$Estimator))
cor.test(log(chao.Site$Observed),log(chao.Site$Estimator))
cor.test(log(chao.Transect$Observed),log(chao.Transect$Estimator))

estimateD.Location=estimateD(agg2.Location,datatype = "abundance", base = "coverage", level = 0.95)
estimateD.Location=subset(estimateD.Location,estimateD.Location$order==0)
estimateD.Site=estimateD(agg2.Site,datatype = "abundance", base = "coverage", level = 0.95)
estimateD.Site=subset(estimateD.Site,estimateD.Site$order==0)
estimateD.Transect=estimateD(agg2.Transect,datatype = "abundance", base = "coverage", level = 0.95)
estimateD.Transect=subset(estimateD.Transect,estimateD.Transect$order==0)
estimateD.all=rbind(estimateD.Site,estimateD.Transect,estimateD.Location)



#######################################################################################
######
###### Code for Table 2
######
#######################################################################################

AllFD=dbFD(traits, agg.scales[,2:134],corr="cailliez")
table.diversity=cbind.data.frame(chao.all, FRic=AllFD$FRic,FEve=AllFD$FEve,FDiv=AllFD$FDiv,FDis=AllFD$FDis,FRed=diversity(agg.scales[,2:134],"simpson")-AllFD$RaoQ)
table.diversity=cbind.data.frame(table.diversity,agg.scales[,135:138])
table.diversity$Lat=census$Lat[match(table.diversity$Site.name,census$Site)]
table.diversity$Lon=census$Lon[match(table.diversity$Site.name,census$Site)]
table.diversity$Lon2=census$Lon2[match(table.diversity$Site.name,census$Site)]
table.diversity$abundance=apply(agg.scales[,2:134],1,sum)
table.diversity$sampling.unit=agg.scales$Group.1

## LMM, only for Site and Transect scale

indices.diversidad=c("log10(Estimator)","FRic","FEve","FDiv","FDis","FRed")
escalas=c("Site","Transect")
combos=expand.grid(indices.diversidad,escalas)
tablon.summary=matrix(nrow=12,ncol=7)


for (w in 1:6)
{
  print(w)
  fmula=formula(paste(combos[w,1],"~1"))
  fmula2=formula(paste(combos[w,1],"~ Realm"))
  lme.site.random=lme(fmula, random=~ 1|Realm/Province,data=table.diversity[table.diversity$Scale==combos[w,2],],na.action = na.omit)
  tablon.summary[w,1:2]=varcomp(lme.site.random, TRUE, F)[1:2]
  tablon.summary[w,4]=varcomp(lme.site.random, TRUE, F)[3]
  lme.site.mixed=lme(fmula2, random=~ 1|Province/Site.name,data=table.diversity[table.diversity$Scale==combos[w,2],],na.action = na.omit)
  atemp=summary(lme.site.mixed)
  tablon.summary[w,5]=atemp$tTable[2]
  tablon.summary[w,6]=atemp$tTable[10]
  tablon.summary[w,7]=r.squaredGLMM(lme.site.mixed)[1]
}

for (w in 7:12)
{
  fmula=formula(paste(combos[w,1],"~1"))
  fmula2=formula(paste(combos[w,1],"~ Realm"))
  lme.site.random=lme(fmula, random=~ 1|Realm/Province/Site.name,data=table.diversity[table.diversity$Scale==combos[w,2],],na.action = na.omit)
  tablon.summary[w,1:4]=varcomp(lme.site.random, TRUE, F)
  lme.site.mixed=lme(fmula2, random=~ 1|Province/Site.name/sampling.unit,data=table.diversity[table.diversity$Scale==combos[w,2],],na.action = na.omit)
  atemp=summary(lme.site.mixed)
  tablon.summary[w,5]=atemp$tTable[2]
  tablon.summary[w,6]=atemp$tTable[10]
  tablon.summary[w,7]=r.squaredGLMM(lme.site.mixed)[1]
}

table2=data.frame(combos, round(tablon.summary,3))
names(table2)=c("scale","index","%var.Realm","%var.Prov","%var.transect","%var.residual","Realm.coeff","Realm.pvalue","PseudoR2_marginal")
write.csv(table2,"MS_Table2.csv")


#######################################################################################
######
###### Code for Fig. 2
######
#######################################################################################

names.provinces.location=factor(table.diversity$Province[table.diversity$Scale=="Location"],levels=c("Southeast Australian Shelf", "Northern New Zealand", 
                                                                                                     "Southern New Zealand","Juan Fernandez and Desventuradas","Warm Temperate Southeastern Pacific"))
names.provinces.site=factor(table.diversity$Province[table.diversity$Scale=="Site"],levels=c("Southeast Australian Shelf", "Northern New Zealand", 
                                                                                             "Southern New Zealand","Juan Fernandez and Desventuradas","Warm Temperate Southeastern Pacific"))
names.provinces.transect=factor(table.diversity$Province[table.diversity$Scale=="Transect"],levels=c("Southeast Australian Shelf", "Northern New Zealand", 
                                                                                                     "Southern New Zealand","Juan Fernandez and Desventuradas","Warm Temperate Southeastern Pacific"))
abb.provinces=c("SAS", "NNZ","SNZ","JFD","WTSP")

color.provinces=data.frame(provinces=levels(as.factor(agg.Location$Province)),
                           color=c("cadetblue1", "coral1","brown4","chartreuse3","deepskyblue3"))
color.provinces=color.provinces[c(3,2,4,1,5),]


tablediv.location=subset(table.diversity,table.diversity$Scale=="Location")
tablediv.location=tablediv.location[levels(names.provinces.location),]

pdf("Jan2022_Fig3.pdf",width=11,height=8.5)
par(mfrow=c(3,3),mar=c(4,4.5,2,1),oma=c(1,1,1,2))
plot(NULL, xlab="Individuals", ylab="Species richness", xlim=c(1, 10000), ylim=c(1, 150),log="xy",cex.axis=1.5,cex.lab=1.5)
for(i in 1:ncol(agg2.Location))
{
  cov.target=cov.Location$iNextEst[[i]]
  cov.target$province=agg.Location$Province[i]
  cov.target$linetype=ifelse(cov.target$method=="interpolated" | cov.target$method=="observed",1,2)
  cov.target$pch=ifelse(cov.target$method=="observed",21,NA)
  cov.target$color=color.provinces$color[match(cov.target$province,color.provinces$provinces)]
  cov.target.int=subset(cov.target,cov.target$linetype==1)
  cov.target.ext=subset(cov.target,cov.target$linetype==2)
  lines(cov.target.int$m,cov.target.int$qD,lty=1,col=cov.target.int$color,lwd=1)
  lines(cov.target.ext$m,cov.target.ext$qD,lty=3,col=cov.target.ext$color,lwd=1)
  points(cov.target$m,cov.target$qD,pch=cov.target$pch,col=cov.target$color,bg=cov.target$color)
}
mtext("(A) Locality",side=3,adj=0,line=0,cex=1.5)
legend("topleft",c("SAS","NNZ","SNZ","JFD","WTSP"),lty=1,col=c("brown4","coral1","chartreuse3",
                                                                "cadetblue1","deepskyblue3"),cex=1.2)
legend("bottomright",c("interpolated","extrapolated","observed"),lty=c(1,2,NA),
       pch=c(NA,NA,21),pt.bg = c(NA,NA,"black"), cex=1.2)

plot(NULL, xlab="Individuals", ylab="", xlim=c(1, 2500), ylim=c(1, 50),log="xy",cex.axis=1.5,cex.lab=1.5)
for(i in 1:ncol(agg2.Site))
{
  cov.target=cov.Site$iNextEst[[i]]
  cov.target$province=agg.Site$Province[i]
  cov.target$linetype=ifelse(cov.target$method=="interpolated" | cov.target$method=="observed",1,2)
  cov.target$pch=ifelse(cov.target$method=="observed",21,NA)
  cov.target$color=color.provinces$color[match(cov.target$province,color.provinces$provinces)]
  cov.target.int=subset(cov.target,cov.target$linetype==1)
  cov.target.ext=subset(cov.target,cov.target$linetype==2)
  lines(cov.target.int$m,cov.target.int$qD,lty=1,col=cov.target.int$color,lwd=1)
  lines(cov.target.ext$m,cov.target.ext$qD,lty=3,col=cov.target.ext$color,lwd=1)
  points(cov.target$m,cov.target$qD,pch=cov.target$pch,col=cov.target$color,bg=cov.target$color)
}
mtext("(B) Site",side=3,adj=0,line=0,cex=1.5)

plot(NULL, xlab="Individuals", ylab="", xlim=c(1, 1000), ylim=c(1, 50),log="xy",cex.axis=1.5,cex.lab=1.5)
for(i in 1:ncol(agg2.Transect))
{
  cov.target=cov.Transect$iNextEst[[i]]
  cov.target$province=agg.Transect$Province[i]
  cov.target$linetype=ifelse(cov.target$method=="interpolated" | cov.target$method=="observed",1,2)
  cov.target$pch=ifelse(cov.target$method=="observed",21,NA)
  cov.target$color=color.provinces$color[match(cov.target$province,color.provinces$provinces)]
  cov.target.int=subset(cov.target,cov.target$linetype==1)
  cov.target.ext=subset(cov.target,cov.target$linetype==2)
  lines(cov.target.int$m,cov.target.int$qD,lty=1,col=cov.target.int$color,lwd=1)
  lines(cov.target.ext$m,cov.target.ext$qD,lty=3,col=cov.target.ext$color,lwd=1)
  points(cov.target$m,cov.target$qD,pch=cov.target$pch,col=cov.target$color,bg=cov.target$color)
}
mtext("(C) Transect",side=3,adj=0,line=0,cex=1.5)
mtext("Rarefaction/Extrapolation",side=4,cex=1.2,line=1)

plot(NULL, xlab="Sample coverage", ylab="Species richness", xlim=c(0.1, 1), ylim=c(1, 150),log="xy",cex.axis=1.5,cex.lab=1.5)
for(i in 1:ncol(agg2.Location))
{
  cov.target=cov.Location$iNextEst[[i]]
  cov.target$province=agg.Location$Province[i]
  cov.target$linetype=ifelse(cov.target$method=="interpolated" | cov.target$method=="observed",1,2)
  cov.target$pch=ifelse(cov.target$method=="observed",21,NA)
  cov.target$color=color.provinces$color[match(cov.target$province,color.provinces$provinces)]
  cov.target.int=subset(cov.target,cov.target$linetype==1)
  cov.target.ext=subset(cov.target,cov.target$linetype==2)
  lines(cov.target.int$SC,cov.target.int$qD,lty=1,col=cov.target.int$color,lwd=1)
  lines(cov.target.ext$SC,cov.target.ext$qD,lty=3,col=cov.target.ext$color,lwd=1)
  points(cov.target$SC,cov.target$qD,pch=cov.target$pch,col=cov.target$color,bg=cov.target$color)
}

plot(NULL, xlab="Sample coverage", ylab="", xlim=c(0.1, 1), ylim=c(1, 50),log="xy",cex.axis=1.5,cex.lab=1.5)
for(i in 1:ncol(agg2.Site))
{
  cov.target=cov.Site$iNextEst[[i]]
  cov.target$province=agg.Site$Province[i]
  cov.target$linetype=ifelse(cov.target$method=="interpolated" | cov.target$method=="observed",1,2)
  cov.target$pch=ifelse(cov.target$method=="observed",21,NA)
  cov.target$color=color.provinces$color[match(cov.target$province,color.provinces$provinces)]
  cov.target.int=subset(cov.target,cov.target$linetype==1)
  cov.target.ext=subset(cov.target,cov.target$linetype==2)
  lines(cov.target.int$SC,cov.target.int$qD,lty=1,col=cov.target.int$color,lwd=1)
  lines(cov.target.ext$SC,cov.target.ext$qD,lty=3,col=cov.target.ext$color,lwd=1)
  points(cov.target$SC,cov.target$qD,pch=cov.target$pch,col=cov.target$color,bg=cov.target$color)
}

plot(NULL, xlab="Sample coverage", ylab="", xlim=c(0.1, 1), ylim=c(1, 50),log="xy",cex.axis=1.5,cex.lab=1.5)
for(i in 1:ncol(agg2.Transect))
{
  cov.target=cov.Transect$iNextEst[[i]]
  cov.target$province=agg.Transect$Province[i]
  cov.target$linetype=ifelse(cov.target$method=="interpolated" | cov.target$method=="observed",1,2)
  cov.target$pch=ifelse(cov.target$method=="observed",21,NA)
  cov.target$color=color.provinces$color[match(cov.target$province,color.provinces$provinces)]
  cov.target.int=subset(cov.target,cov.target$linetype==1)
  cov.target.ext=subset(cov.target,cov.target$linetype==2)
  lines(cov.target.int$SC,cov.target.int$qD,lty=1,col=cov.target.int$color,lwd=1)
  lines(cov.target.ext$SC,cov.target.ext$qD,lty=3,col=cov.target.ext$color,lwd=1)
  points(cov.target$SC,cov.target$qD,pch=cov.target$pch,col=cov.target$color,bg=cov.target$color)
}
mtext("Coverage",side=4,cex=1.2,line=1)

plot(Estimator~names.provinces.location,data=table.diversity[table.diversity$Scale=="Location",],
     col=color.provinces$color,xaxt="n",xlab="",ylab="Species richness (Chao1)",
     cex.axis=1.5,cex.lab=1.5,ylim=c(1,150),log="y")
lines(x=c(3.5,3.5),y=c(1,150),lty=2,col="gray20")
axis(1,las=1,cex.axis=1.2,labels=abb.provinces,at=1:5)
text("West",x=1,y=1.5,cex=1.5,col="gray40")
text("East",x=5,y=1.5,cex=1.5,col="gray40")

plot(Estimator~names.provinces.site,data=table.diversity[table.diversity$Scale=="Site",],
     col=color.provinces$color,xaxt="n",xlab="",ylab="",
     cex.axis=1.5,cex.lab=1.5,ylim=c(1,150),log="y")
lines(x=c(3.5,3.5),y=c(1,150),lty=2,col="gray20")
axis(1,las=1,cex.axis=1.2,labels=abb.provinces,at=1:5)
text("West",x=1,y=1.5,cex=1.5,col="gray40")
text("East",x=5,y=1.5,cex=1.5,col="gray40")

plot(Estimator~names.provinces.transect,data=table.diversity[table.diversity$Scale=="Transect",],
     ylim=c(1,150),col=color.provinces$color,xaxt="n",xlab="",ylab="",cex.axis=1.5,cex.lab=1.5,log="y")
axis(1,las=1,cex.axis=1.2,labels=abb.provinces,at=1:5)
lines(x=c(3.5,3.5),y=c(1,150),lty=2,col="gray20")
text("West",x=1,y=1.5,cex=1.5,col="gray40")
text("East",x=5,y=1.5,cex=1.5,col="gray40")
mtext("Asymptotic spp. richness",side=4,cex=1.2,line=1)
dev.off()



#######################################################################################
######
###### Figure S1
######
#######################################################################################

obs.FDs=ddply(table.diversity,~Scale+Province,summarize,mFRic=mean(FRic,na.rm=T),
            mFEve=mean(FEve,na.rm=T),mFDiv=mean(FDiv,na.rm=T),mFDis=mean(FDis,na.rm=T),
            mFRed=mean(FRed,na.rm=T))

obsCI.FDs=ddply(table.diversity,~Scale+Province,summarize,mFRic=1.96*sd(FRic,na.rm=T)/sqrt(length(FRic)),
                mFEve=1.96*sd(FEve,na.rm=T)/sqrt(length(FEve)),mFDiv=1.96*sd(FDiv,na.rm=T)/sqrt(length(FDiv)),
                mFDis=1.96*sd(FDis,na.rm=T)/sqrt(length(FDis)),mFRed=1.96*sd(FRed,na.rm=T)/sqrt(length(FRed)))

obs.FDs
obs.FDs$Lon=table.diversity$Lon[match(obs.FDs$Province,table.diversity$Province)]
obs.FDs=obs.FDs[order(obs.FDs$Scale,obs.FDs$Lon),]

obsCI.FDs$Lon=table.diversity$Lon[match(obsCI.FDs$Province,table.diversity$Province)]
obsCI.FDs=obsCI.FDs[order(obsCI.FDs$Scale,obsCI.FDs$Lon),]


plot.indices=c("FRic","FEve","FDiv","FDis","FRed")
plot.scales=c("Site","Transect")
plot.indscales=expand.grid(plot.scales,plot.indices)
names(plot.indscales)=c("Scale","FD")
scale.title1=c("(A) Site","","","","")
scale.title2=c("(B) Transect","","","","")
fd.titles=c("functional richness","functional eveness","functional divergence","functional dispersion",
            "functional redundacy")
abb.provinces.titles=matrix(ncol=5,nrow=10)
abb.provinces.titles[6:10,]=rep((abb.provinces),each=5)
abb.provinces.titles[is.na(abb.provinces.titles)]=""
ylab.titles=c("FRic","FEve","FDiv","FDis","FRed")
realm.titles.W=c("","","",""," West")
realm.titles.E=c("","","","","East ")

pdf("Jan2022_FigS1.pdf",width=11,height=6.8)
par(mfrow=c(2,5),mar=c(4,4,1,1),oma=c(1,1,2,1),new=F)
for(t in 1:5){
  plotCI(6:10,obs.FDs[1:5,t+2],uiw=obsCI.FDs[6:10,t+2],xaxt="n",xlab="",ylab=ylab.titles[t],scol=color.provinces$color)
  points(6:10,obs.FDs[1:5,t+2],pch=21,bg=color.provinces$color,cex=2)
  par(new = TRUE)  
  plot(6:10,1:5,xaxt="n",cex=0,axes=F,xlab="",ylab="")
  axis(1,las=1,at=6:10,labels=abb.provinces,cex.axis=0.9)
  mtext(fd.titles[t],cex=0.8,adj=0)
  mtext(scale.title1[t],cex=1,adj=0,line=1)
  mtext(realm.titles.W[t],cex=0.7,adj=0,line=-1,col="gray20")
  mtext(realm.titles.E[t],cex=0.7,adj=1,line=-1,col="gray20")
  lines(x=c(8.5,8.5),y=c(-25,25),lty=2,col="red")
  }
for(t in 1:5){
  plotCI(11:15,obs.FDs[11:15,t+2],uiw=obsCI.FDs[11:15,t+2],xaxt="n",xlab="",ylab=ylab.titles[t],scol=color.provinces$color)
  points(11:15,obs.FDs[11:15,t+2],pch=21,bg=color.provinces$color,cex=2)
  par(new = TRUE)  
  plot(1:5,1:5,xaxt="n",cex=0,axes=F,xlab="",ylab="")
  axis(1,las=1,at=1:5,labels=abb.provinces,cex.axis=0.9)
  mtext(scale.title2[t],cex=1,adj=0,line=1)
  mtext(realm.titles.W[t],cex=0.7,adj=0,line=-1,col="gray20")
  mtext(realm.titles.E[t],cex=0.7,adj=1,line=-1,col="gray20")
  lines(x=c(3.5,3.5),y=c(-25,25),lty=2,col="red")
  }
dev.off()


#######################################################################################
######
###### Figure 4, species traits weighted by abundance
######
#######################################################################################
agg.Province=aggregate(census[,9:141],by=list(census$Province),sum)
rownames(agg.Province)=agg.Province[,1];agg.Province=agg.Province[,-1]

cwmDOM=functcomp(traits, as.matrix(agg.Province),CWM.type="dom")
cwmDOM$Lon=census$Lon[match(rownames(cwmDOM),census$Province)]  
cwmDOM=cwmDOM[order(cwmDOM$Lon),]
  
cwm=functcomp(traits, as.matrix(agg.Province),CWM.type="all")
cwm$Lon=census$Lon[match(rownames(cwm),census$Province)]  
cwm=cwm[order(cwm$Lon),]

cwm1=cwm[,1:4]
cwm2=cwm[,5:10]
cwm3=cwm[,11:21]
cwm4=cwm[,24:27]
cwm5=cwm[,28:30]
cwm6=cwm[,31:33]
cwm7=cwm[,22]
cwm8=cwm[,23]

#### isolate, maybe get rid of it
spp.Province=aggregate(ifelse(census[,9:141]>0,1,0),by=list(census$Province),sum)
rownames(spp.Province)=spp.Province$Group.1;spp.Province=spp.Province[,-1]
spp.Province=ifelse(spp.Province>0,1,0)
traits.spp.Province=data.frame(t(spp.Province))
traits.spp.Province2=as.data.frame(as.table(as.matrix(traits.spp.Province[,1:5])))
traits.spp.Province2=cbind(traits.spp.Province2,traits[match(traits.spp.Province2$Var1,rownames(traits)),])
traits.spp.Province2=subset(traits.spp.Province2,traits.spp.Province2$Freq>0)
names(traits.spp.Province2)[1:2]=c("Species","Province")
#####

spp.Province.ind=as.data.frame(t(agg.Province))
spp.Province.ind=as.data.frame(as.table(as.matrix(spp.Province.ind)))
spp.Province.ind$TrophicLevel=rep(traits$TrophicLevel,5)
spp.Province.ind$MaxLenght=rep(traits$MaxLenght,5)
spp.Province.ind[spp.Province.ind==0]=NA; spp.Province.ind=na.exclude(spp.Province.ind)

a=NULL
b=NULL
c=NULL
for(q in 1:nrow(spp.Province.ind))
{
  print(q)
  a=c(a,rep(spp.Province.ind$TrophicLevel[q],spp.Province.ind$Freq[q]))
  b=c(b,rep(spp.Province.ind$MaxLenght[q],spp.Province.ind$Freq[q]))
  c=c(c,rep(as.character(spp.Province.ind$Var2[q]),spp.Province.ind$Freq[q]))
}
ind.traits=data.frame(Province=c,TrophicLevel=a,MaxLenght=b)

beanplot(ind.traits$TrophicLevel~as.factor(ind.traits$Province),ylim=c(1,6),col=c("firebrick2"),maxstripline=0,
         xaxt="n", ylab="Level",cex.lab=1.2,cex.axis=1.2,what=c(0,1,0,0))




col.TrophicGroup=c("darkolivegreen","darkolivegreen4","darkolivegreen3","darkolivegreen2")
col.Trophic.category=c("brown4","brown1","coral", "cornsilk4","cornsilk3","cornsilk")
col.FineTrophicGroup=c("black","gray30","gray60","gray90","white",
                       "orange","yellow4","yellow3","yellow2","yellow1","lightyellow")
col.Water.Column.Position=c("turquoise4","turquoise3","turquoise2","turquoise")
col.Active.Period=c("slategray4","slategray3","slategray2","slategray1")
col.Gregariousness=c("navajowhite4","navajowhite3","navajowhite")

pdf("Jan2022_Fig4.pdf",width=11,height=8.5)
par(mfrow=c(2,4),mar=c(5,5,2,1),oma=c(1,1,1,1))
a1=barplot(t(cwm1),names.arg=abb.provinces,col=col.TrophicGroup,ylab="Proportion",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,las=1,at=a1,labels=abb.provinces,cex.axis=0.8)
box()
legend(x=.26,y=0.65, legend=levels(traits$TrophicGroup),fill=col.TrophicGroup,bg="white",cex=0.6)
mtext("(A) Trophic group",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.7,3.7),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.6,y=0.1,cex=1,col="gray80")
text("East",x=5.5,y=0.1,cex=1,col="gray80")

a2=barplot(t(cwm2),names.arg=abb.provinces,col=col.Trophic.category, ylab="Proportion",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,las=1,at=a2,labels=abb.provinces,cex.axis=0.8)
box()
legend(x=0.26,y=0.85, legend=levels(traits$Trophic.category),fill=col.Trophic.category,bg="white",cex=0.6)
mtext("(B) Trophic category",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.7,3.7),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.6,y=0.1,cex=1,col="gray80")
text("East",x=5.5,y=0.1,cex=1,col="gray80")

a3=barplot(t(cwm3),names.arg=abb.provinces,col=col.FineTrophicGroup, ylab="Proportion",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,las=1,at=a3,labels=abb.provinces,cex.axis=0.8)
box()
lines(x=c(3.7,3.7),y=c(-10,50),lty=2,col="gray20")
legend(x=0.5,y=0.95, legend=levels(traits$FineTrophicGroup),fill=col.FineTrophicGroup,bg="white",cex=0.6,ncol=2)
mtext("(C) Fine Trophic Group",side=3, line=0.5,adj=0,cex=1)
text("West",x=0.6,y=0.1,cex=1,col="black")
text("East",x=5.5,y=0.1,cex=1,col="black")

a4=barplot(t(cwm4),names.arg=abb.provinces,col=col.Water.Column.Position, ylab="Proportion",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,las=1,at=a4,labels=abb.provinces,cex.axis=0.8)
box()
legend(x=0.26,y=0.65, legend=levels(traits$Water.Column.Position),fill=col.Water.Column.Position,bg="white",cex=0.6)
mtext("(D) Water column position",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.7,3.7),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.6,y=0.1,cex=1,col="black")
text("East",x=5.5,y=0.1,cex=1,col="black")

a5=barplot(t(cwm5),names.arg=abb.provinces,col=col.Active.Period, ylab="Proportion",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,las=1,at=a5,labels=abb.provinces,cex.axis=0.8)
box()
legend(x=0.26,y=0.45, legend=levels(traits$Active.Period),fill=col.Active.Period,bg="white",cex=0.6)
mtext("(E) Active period",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.7,3.7),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.6,y=0.1,cex=1,col="black")
text("East",x=5.5,y=0.1,cex=1,col="black")

a6=barplot(t(cwm6),names.arg=abb.provinces,col=col.Gregariousness, ylab="Proportion",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,las=1,at=a6,labels=abb.provinces,cex.axis=0.8)
box()
legend(x=0.26,y=0.45, legend=levels(traits$Gregariousness),fill=col.Gregariousness,bg="white",cex=0.6)
mtext("(F) Gregariousness",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.7,3.7),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.6,y=0.1,cex=1,col="gray80")
text("East",x=5.5,y=0.1,cex=1,col="gray80")

beanplot(spp.Province.ind$TrophicLevel~as.factor(spp.Province.ind$Var2),ylim=c(1,6),col=c("firebrick2"),maxstripline=0,
         xaxt="n", ylab="Level",cex.lab=1.2,cex.axis=1.2,what=c(0,1,0,0))
points(x=1:5,y=cwm7,pch=3,cex=4,lwd=3)
axis(1,las=1,at=1:5,labels=abb.provinces,cex.axis=0.8)
mtext("(G) Trophic level",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.5,3.5),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.8,y=1.4,cex=1,col="gray20")
text("East",x=5,y=1.4,cex=1,col="gray20")

beanplot(spp.Province.ind$MaxLenght~as.factor(spp.Province.ind$Var2),log="y",col=c("aquamarine4"),maxstripline=0,
         xaxt="n", ylab="body length (cm)",cex.lab=1.2,cex.axis=1.2,what=c(0,1,0,0))
points(x=1:5,y=cwm8,pch=3,cex=4,lwd=3)
axis(1,las=1,at=1:5,labels=abb.provinces,cex.axis=0.8)
mtext("(H) Body size",side=3, line=0.5,adj=0,cex=1)
lines(x=c(3.5,3.5),y=c(-10,50),lty=2,col="gray20")
text("West",x=0.8,y=2,cex=1,col="gray20")
text("East",x=5,y=2,cex=1,col="gray20")
dev.off()

#######################################################################################
######
###### Figure 5, CAPSCALE and PERMANOVA
######
#######################################################################################

functional.abb=read.csv("functionalgroups_abb.csv",T,sep=";")

cwm.Site=functcomp(traits, as.matrix(agg.Site[,2:134]),CWM.type="all")
cwm.Site$TrophicLevel=(cwm.Site$TrophicLevel-min(cwm.Site$TrophicLevel))/(max(cwm.Site$TrophicLevel)-min(cwm.Site$TrophicLevel))
cwm.Site$MaxLenght=(cwm.Site$MaxLenght-min(cwm.Site$MaxLenght))/(max(cwm.Site$MaxLenght)-min(cwm.Site$MaxLenght))

abb.provinces=c("SAS", "NNZ","SNZ","JFD","WTSP")
color.provinces=data.frame(provinces=levels(as.factor(agg.Location$Province)),
                           color=c("cadetblue1", "coral1","brown4","chartreuse3","deepskyblue3"))
color.provinces$realm=c("east","west","west","west","east")
fact.province=data.frame(Province=agg.Site$Province)
fact.province$color=color.provinces$color[match(fact.province$Province,color.provinces$provinces)]
fact.province$realm=color.provinces$realm[match(fact.province$Province,color.provinces$provinces)]

## CAPSCALE
cap.Site=capscale(cwm.Site ~factor(Province),data = fact.province,distance = "bray")
eigen.Site=cap.Site %>%eigenvals(model="constrained") %>% 
  decostand(., "total", MARGIN = 2) %>% round(2)
scores.stand=as.data.frame(scores(cap.Site)$sites) %>%decostand("max")
traits.loads=summary(cap.Site)$species[,1:2] %>% data.frame() %>%
  filter(abs(CAP1) > 0.1 | abs(CAP2) > 0.1)
top.threshold=0.2
traits.loads.all=summary(cap.Site)$species[,1:2] %>% data.frame()
traits.loads.all$top=ifelse(abs(traits.loads.all$CAP1)>top.threshold | abs(traits.loads.all$CAP2)>top.threshold,
"top","non")
traits.loads.all$abb=functional.abb$GroupCode[match(rownames(traits.loads.all),functional.abb$Name)]
traits.loads.all[traits.loads.all$top=="top",]

pdf("Jan2022_Fig5.pdf",height=5,width = 7)
plot(scores.stand,col="white",xlab="CAP1 (74 %)",ylab="CAP2 (12 %)")
points(scores.stand, bg=fact.province$color, col=fact.province$color,pch = 21,cex=0.8)
ordihull(scores.stand,groups=fact.province$Province,draw="polygon",col=color.provinces$color,
         border=color.provinces$color,label=F,lwd=0.5)
legend("topright",abb.provinces,pch=21,pt.bg=color.provinces$color[c(3,2,4,1,5)],col=color.provinces$color[c(3,2,4,1,5)])
arrows(0,0,traits.loads.all[,1], traits.loads.all[,2], length=0.1, lty=1, 
       col=ifelse(traits.loads.all$top=="top","gray20","gray80"),
       lwd=ifelse(traits.loads.all$top=="top",1,0.5),angle=15, code=2)
text(jitter(traits.loads.all[,1][traits.loads.all$top=="top"],amount=0.05), jitter(traits.loads.all[,2][traits.loads.all$top=="top"],amount=0.05),
     traits.loads.all$abb[traits.loads.all$top=="top"],cex=0.7,pos=c(1,2,3,4))
dev.off()

## PERMANOVA
perm=how(nperm = 10000)
setBlocks(perm)=with(fact.province, realm)
model.perma=adonis2(cwm.Site~Province,data=fact.province,permutations = perm)
model.perma

dist.trait=vegdist(cwm.Site)
pairwise.perm.manova(dist.trait,fact.province$Province,nperm=10000,p.method="BH")




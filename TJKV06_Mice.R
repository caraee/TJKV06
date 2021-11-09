library("here")
library("tidyverse")
library("extrafont")
library("scales")
library("grid")
library("gtable")
library("data.table")
library("Amelia")
library("cmdstanr")
library("brms")
library("future")
library("gt")

weekly_monitoring<-read.table("weekly_monitoring.txt",sep="\t",stringsAsFactors = F, 
                              header=T, check.names=F,as.is=T)
GSIS<-read.table("GSIS.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)
ArgTT<-read.table("ArgTT.txt",sep="\t",stringsAsFactors = F, 
                  header=T, check.names=F,as.is=T)
animalTable<-read.table("animalTable.txt",sep="\t",stringsAsFactors = F, 
                        header=T, check.names=F,as.is=T)

GroupPalette<-c(  "#d03d00",
                  "#7a59ff",
                  "#ffbd47",
                  "#e4007b",
                  "#003d25",
                  "#908595") 
GroupShapes<-c(8,23,22,25,24,21)

IndivPalette<-c("#5c008d","#afe500","#3f3ae2","#fff32d","#c858ff","#beff72",
                "#f200ab","#50ffd1","#ff253b","#01bb7a","#8f0090","#dac801",
                "#012c8e","#ffef6f","#120040","#baff97","#ee005e","#00c5b9",
                "#d60139","#008240","#fe8fff","#81a500","#bc9aff","#ff941c",
                "#75a2ff","#d46c00","#00b3ee","#ff6732","#004d86","#ffb573",
                "#07000a","#b6ffdc","#860047","#004f09","#b2bfff","#9d4b00",
                "#aff0ff","#721400","#00919f","#ff8271","#002e17","#fff3f9",
                "#324100","#ffbedc","#006150","#ffeab2","#004858","#5e4c00")

##### Weekly Monitoring #####
#No missing values
animalTable$Group<-gsub("PEC:","",animalTable$Group)
# animalTable$Sex<-gsub("Male","♂",animalTable$Sex)
# animalTable$Sex<-gsub("Female","♀", animalTable$Sex)

weekly_monitoring$Sex<-sapply(weekly_monitoring$AnimalID,function(a){
  animalTable$Sex[animalTable$AnimalID==a]})
weekly_monitoring$Group<-sapply(weekly_monitoring$AnimalID,function(a){
  animalTable$Group[animalTable$AnimalID==a]})
weekly_monitoring$Date<-as.POSIXct(strptime(weekly_monitoring$Date,"%Y/%m/%d",tz="America/Vancouver"))
weekly_monitoring$Days <- sapply(1:nrow(weekly_monitoring),function(i){
  surgDate<-as.POSIXct(strptime("2018/02/20","%Y/%m/%d",tz="America/Vancouver"))
  d<-weekly_monitoring$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
weekly_monitoring$Weeks<-round(weekly_monitoring$Days/7)
weekly_monitoring<-weekly_monitoring[order(weekly_monitoring$Weeks,
                                           weekly_monitoring$Sex,weekly_monitoring$Group,weekly_monitoring$AnimalID),]
weekly_monitoring$AnimalID <- factor(weekly_monitoring$AnimalID, levels = as.character(unique(weekly_monitoring$AnimalID)))
weekly_monitoring$Group<-factor(weekly_monitoring$Group,levels = unique(weekly_monitoring$Group))
weekly_monitoring$Sex<-factor(weekly_monitoring$Sex,levels = unique(weekly_monitoring$Sex))
weekly_monitoring$FullGroup<-paste(weekly_monitoring$Sex,"_",weekly_monitoring$Group,sep="")
weekly_monitoring$FullGroup<-factor(weekly_monitoring$FullGroup,levels = c("Female_FP",
                                                                           "Female_KC",
                                                                           "Female_TC",
                                                                           "Male_FP",
                                                                           "Male_KC",
                                                                           "Male_TC"))
names(GroupPalette)<-levels(weekly_monitoring$FullGroup)
names(GroupShapes)<-levels(weekly_monitoring$FullGroup)
names(IndivPalette)<-c(levels(weekly_monitoring$AnimalID))

#animal died during surgery
weekly_monitoring<-weekly_monitoring[weekly_monitoring$AnimalID!="TJKV06M.09",]

##### Weekly Monitoring: Statistics #####
###### Blood Glucose ######
weekly_monitoring <- transform(weekly_monitoring, DaysNom = factor(Days))

fit_BG <- brm(BG~FullGroup*DaysNom+(DaysNom|AnimalID), 
              data = weekly_monitoring, 
              control = list(adapt_delta = 0.99),
              chains = 12, 
              iter = 6000,
              cores=4)
save(fit_BG,file="fit_BG.RData")

load("./fit_BG.RData")
newdata = data.frame(Group = rep(levels(weekly_monitoring$Group),
                                 38)|>str_sort(),
                     Sex=rep(levels(weekly_monitoring$Sex),19)|>rep(3),
                     DaysNom = rep(levels(weekly_monitoring$DaysNom),each=2)|>rep(3),
                     Days = rep(unique(weekly_monitoring$Days),each=2)|>rep(3))
newdata$FullGroup<-paste(newdata$Sex, newdata$Group,sep="_")

fit_weeklyBG <- fitted(fit_BG,
                       newdata = newdata, 
                       re_formula = NA) # extract the full MCMC

ff_weeklyBG <- fit_weeklyBG |>
  as_tibble() |>
  bind_cols(newdata)

ff_weeklyBG<-ff_weeklyBG[ff_weeklyBG$Days>=-2,]

ff_weeklyBG$FullGroup<-factor(ff_weeklyBG$FullGroup,levels=levels(weekly_monitoring$FullGroup))
ff_weeklyBG$Group<-factor(ff_weeklyBG$Group,levels=levels(weekly_monitoring$Group))
ff_weeklyBG$Sex<-factor(ff_weeklyBG$Sex,levels=levels(weekly_monitoring$Sex))
ff_weeklyBG$DaysNom<-factor(ff_weeklyBG$DaysNom,levels=levels(weekly_monitoring$DaysNom))
ff_weeklyBG$Weeks<-round(ff_weeklyBG$Days/7)

write_csv(ff_weeklyBG,file="./weeklyBG.csv")

p_weeklyBG<-ggplot(weekly_monitoring[weekly_monitoring$Days>=-2,],
                   aes(x=Weeks,y=BG,group=Group))+
  geom_line(aes(colour=FullGroup,group=AnimalID),linetype=2,size=0.3,alpha=0.8)+
  facet_grid(Sex~.,scales="free_x",space="free_x")+
  labs(x="Weeks post implant",y="Blood Glucose (mM)")+
  scale_y_continuous(limits=c(0,15))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  scale_x_continuous(limits=c(0,28),
                     breaks=pretty_breaks(n=6))+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank())+
  # ,
  #       legend.position = c(0.5,0.50),
  #       legend.key.width = unit(0.5,"cm"))+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = ff_weeklyBG,
             aes(y = Estimate,shape=FullGroup,
                 color=FullGroup,fill=FullGroup),
             na.rm=T,size=2,alpha=0.7)+
  geom_smooth(data = ff_weeklyBG,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup, colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1)

p_weeklyBG
ggsave("BG_Bayes.png",path="./Figures/weekly",width = 20, height = 12, units = "cm",dpi=300)


###### Comparisons #####
BGtab1<-rbind(
  hypothesis(fit_BG,"DaysNom10+FullGroupFemale_TC:DaysNom10-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom15+FullGroupFemale_TC:DaysNom15-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom22+FullGroupFemale_TC:DaysNom22-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom29+FullGroupFemale_TC:DaysNom29-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom36+FullGroupFemale_TC:DaysNom36-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom43+FullGroupFemale_TC:DaysNom43-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom50+FullGroupFemale_TC:DaysNom50-DaysNomM2-FullGroupFemale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom64+FullGroupFemale_TC:DaysNom64-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom78+FullGroupFemale_TC:DaysNom78-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom92+FullGroupFemale_TC:DaysNom92-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom99+FullGroupFemale_TC:DaysNom99-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom130+FullGroupFemale_TC:DaysNom130-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom143+FullGroupFemale_TC:DaysNom143-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom164+FullGroupFemale_TC:DaysNom164-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom183+FullGroupFemale_TC:DaysNom183-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom191+FullGroupFemale_TC:DaysNom191-DaysNomM2-FullGroupFemale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom10+FullGroupFemale_KC:DaysNom10-DaysNomM2-FullGroupFemale_KC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom15+FullGroupFemale_KC:DaysNom15-DaysNomM2-FullGroupFemale_KC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom22+FullGroupFemale_KC:DaysNom22-DaysNomM2-FullGroupFemale_KC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom29+FullGroupFemale_KC:DaysNom29-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom36+FullGroupFemale_KC:DaysNom36-DaysNomM2-FullGroupFemale_KC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom50+FullGroupFemale_KC:DaysNom50-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom64+FullGroupFemale_KC:DaysNom64-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom78+FullGroupFemale_KC:DaysNom78-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom92+FullGroupFemale_KC:DaysNom92-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom99+FullGroupFemale_KC:DaysNom99-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom130+FullGroupFemale_KC:DaysNom130-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom143+FullGroupFemale_KC:DaysNom143-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom164+FullGroupFemale_KC:DaysNom164-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom183+FullGroupFemale_KC:DaysNom183-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom191+FullGroupFemale_KC:DaysNom191-DaysNomM2-FullGroupFemale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom191-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom183-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom164-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom143-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom130-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom99-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom92-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom78-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom64-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom50-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom43-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom36-DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom10+FullGroupMale_TC:DaysNom10-DaysNomM2-FullGroupMale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom15+FullGroupMale_TC:DaysNom15-DaysNomM2-FullGroupMale_TC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom22+FullGroupMale_TC:DaysNom22-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom29+FullGroupMale_TC:DaysNom29-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom36+FullGroupMale_TC:DaysNom36-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom43+FullGroupMale_TC:DaysNom43-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom50+FullGroupMale_TC:DaysNom50-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom64+FullGroupMale_TC:DaysNom64-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom78+FullGroupMale_TC:DaysNom78-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom92+FullGroupMale_TC:DaysNom92-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom99+FullGroupMale_TC:DaysNom99-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom130+FullGroupMale_TC:DaysNom130-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom143+FullGroupMale_TC:DaysNom143-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom164+FullGroupMale_TC:DaysNom164-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom183+FullGroupMale_TC:DaysNom183-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom191+FullGroupMale_TC:DaysNom191-DaysNomM2-FullGroupMale_TC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom10+FullGroupMale_KC:DaysNom10-DaysNomM2-FullGroupMale_KC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom15+FullGroupMale_KC:DaysNom15-DaysNomM2-FullGroupMale_KC:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom22+FullGroupMale_KC:DaysNom22-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom29+FullGroupMale_KC:DaysNom29-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom36+FullGroupMale_KC:DaysNom36-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom43+FullGroupMale_KC:DaysNom43-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom50+FullGroupMale_KC:DaysNom50-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom64+FullGroupMale_KC:DaysNom64-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom78+FullGroupMale_KC:DaysNom78-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom92+FullGroupMale_KC:DaysNom92-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom99+FullGroupMale_KC:DaysNom99-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom130+FullGroupMale_KC:DaysNom130-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom143+FullGroupMale_KC:DaysNom143-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom164+FullGroupMale_KC:DaysNom164-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom183+FullGroupMale_KC:DaysNom183-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom191+FullGroupMale_KC:DaysNom191-DaysNomM2-FullGroupMale_KC:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom10+FullGroupMale_FP:DaysNom10-DaysNomM2-FullGroupMale_FP:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom15+FullGroupMale_FP:DaysNom15-DaysNomM2-FullGroupMale_FP:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom15+FullGroupMale_FP:DaysNom22-DaysNomM2-FullGroupMale_FP:DaysNomM2=0")$hypothesis,
  hypothesis(fit_BG,"DaysNom29+FullGroupMale_FP:DaysNom29-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom36+FullGroupMale_FP:DaysNom36-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom43+FullGroupMale_FP:DaysNom43-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom50+FullGroupMale_FP:DaysNom50-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom64+FullGroupMale_FP:DaysNom64-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom78+FullGroupMale_FP:DaysNom78-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom92+FullGroupMale_FP:DaysNom92-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom99+FullGroupMale_FP:DaysNom99-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom130+FullGroupMale_FP:DaysNom130-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom143+FullGroupMale_FP:DaysNom143-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom164+FullGroupMale_FP:DaysNom164-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom183+FullGroupMale_FP:DaysNom183-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis,
  hypothesis(fit_BG,"DaysNom191+FullGroupMale_FP:DaysNom191-DaysNomM2-FullGroupMale_FP:DaysNomM2<0")$hypothesis
)

BGtab2<-rbind(
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom36-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom36>0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom50-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom50>0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom64-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom64=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom78-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom78<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom92-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom92<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom99-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom99=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom130-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom130=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_KC+FullGroupFemale_KC:DaysNom22)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_KC+FullGroupFemale_KC:DaysNom29)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_KC+FullGroupFemale_KC:DaysNom36)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_KC+FullGroupFemale_KC:DaysNom50)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_KC+FullGroupFemale_KC:DaysNom64)=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom78=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom92=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom99<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom130<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom143<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom164<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom183<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_KC+FullGroupFemale_KC:DaysNom191<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_TC+FullGroupFemale_TC:DaysNom22)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_TC+FullGroupFemale_TC:DaysNom29)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_TC+FullGroupFemale_TC:DaysNom36)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_TC+FullGroupFemale_TC:DaysNom50)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupFemale_TC+FullGroupFemale_TC:DaysNom64)<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom78<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom92<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom99<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom130<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom143<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom164<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom183<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom191<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom10-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom10>0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom15-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom15>0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom22-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom22=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom29-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom29=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom36-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom36=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom50-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom50=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom64-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom64=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom78-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom78=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom92-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom92=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom99-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom99=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom130-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom130=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom143-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom143=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom164-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom164=0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom183-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom183<0")$hypothesis,
  hypothesis(fit_BG,"FullGroupMale_TC+FullGroupMale_TC:DaysNom191-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom191<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom10)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom10)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom15)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom15)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom22)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom22)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom29)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom29)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom36)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom36)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom50)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom50)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom64)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom64)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom78)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom78)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom92)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom92)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom99)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom99)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom130)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom130)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom143)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom143)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom164)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom164)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom183)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom183)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_TC+FullGroupMale_TC:DaysNom191)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom191)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom10)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom10)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom15)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom15)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom22)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom22)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom29)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom29)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom36)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom36)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom43)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom43)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom50)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom50)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom64)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom64)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom78)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom78)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom92)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom92)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom99)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom99)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom130)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom130)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom143)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom143)<0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom164)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom164)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom183)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom183)=0")$hypothesis,
  hypothesis(fit_BG,"(FullGroupMale_KC+FullGroupMale_KC:DaysNom191)-
           (FullGroupMale_FP+FullGroupMale_FP:DaysNom191)=0")$hypothesis)

###### Body Weight #######
fit_BW <- brm(BW~FullGroup*DaysNom+(DaysNom|AnimalID), 
              data = weekly_monitoring, 
              chains = 12, 
              iter = 6000,
              cores=4)
#re-run later
save(fit_BW,file="fit_BW.RData")

load("fit_BW.RData")

fit_weeklyBW <- fitted(fit_BW,
                       newdata = newdata, 
                       re_formula = NA) # extract the full MCMC

ff_weeklyBW <- fit_weeklyBW |>
  as_tibble() |>
  bind_cols(newdata)

ff_weeklyBW$FullGroup<-factor(ff_weeklyBW$FullGroup,levels=levels(weekly_monitoring$FullGroup))
ff_weeklyBW$Group<-factor(ff_weeklyBW$Group,levels=levels(weekly_monitoring$Group))
ff_weeklyBW$Sex<-factor(ff_weeklyBW$Sex,levels=levels(weekly_monitoring$Sex))
ff_weeklyBW$DaysNom<-factor(ff_weeklyBW$DaysNom,levels=levels(weekly_monitoring$DaysNom))
ff_weeklyBW$Weeks<-round(ff_weeklyBW$Days/7)

write_csv(ff_weeklyBW,file="./weeklyBW.csv")

p_weeklyBW<-ggplot(weekly_monitoring[!is.na(weekly_monitoring$BW),],
          aes(x=Weeks,y=BW,group=FullGroup))+
  facet_grid(Sex~.,scales="free")+
  geom_line(aes(colour=FullGroup,group=AnimalID),linetype=2,size=0.3,alpha=0.8,
            show.legend = F)+
  labs(x="Weeks post implant",y="Body Weight (g)")+
  scale_y_continuous(limits=c(0,35))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  scale_x_continuous(breaks=pretty_breaks(n=6))+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Male: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Male: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Male: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Male: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Male: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Male: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: TheraCyte"))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank(),
        legend.background = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = ff_weeklyBW,
             aes(y = Estimate,shape=FullGroup,
                 color=FullGroup,fill=FullGroup),
             na.rm=T,size=2,alpha=0.7,
             show.legend = F)+
  geom_smooth(data = ff_weeklyBW,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup, colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1,
              show.legend = F)

p_weeklyBW
ggsave("BW_Bayes.png",path="./Figures/weekly",width = 20, height = 12, units = "cm",dpi=300)

BWtab1<-rbind(
  hypothesis(fit_BW,"FullGroupFemale_TC+FullGroupFemale_TC:DaysNom15-
           FullGroupFemale_KC-FullGroupFemale_KC:DaysNom15=0")$hypothesis,
  hypothesis(fit_BW,"(FullGroupFemale_TC+FullGroupFemale_TC:DaysNom15)=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom183-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom183<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom191-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom191<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom10-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom10>0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom15-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom15>0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom22-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom22>0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom29-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom29>0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom36-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom36>0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom43-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom43>0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom50-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom50=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom64-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom64=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom78-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom78=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom92-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom92=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom99-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom99=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom130-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom130=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom143-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom143=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom164-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom164=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom183-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom183=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom191-
           FullGroupMale_KC-FullGroupMale_KC:DaysNom191=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom10-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom10=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom15-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom15=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom22-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom22=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom29-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom29<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom36-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom36<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom43-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom43<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom50-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom50<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom64-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom64<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom78-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom78<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom92-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom92<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom99-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom99<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom130-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom130<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom143-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom143<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom164-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom164<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom183-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom183<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_KC+FullGroupMale_KC:DaysNom191-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom191<0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom10-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom10=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom15-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom15=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom22-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom22=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom29-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom29=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom36-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom36=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom43-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom43=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom50-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom50=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom64-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom64=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom78-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom78=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom92-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom92=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom99-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom99=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom130-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom130=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom143-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom143=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom164-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom164=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom183-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom183=0")$hypothesis,
  hypothesis(fit_BW,"FullGroupMale_TC+FullGroupMale_TC:DaysNom191-
           FullGroupMale_FP-FullGroupMale_FP:DaysNom191=0")$hypothesis)

#### GSIS: Imputations and Figures ####
GSIS$Sex<-sapply(GSIS$AnimalID,function(a){
  animalTable$Sex[animalTable$AnimalID==a]})
GSIS$Group<-sapply(GSIS$AnimalID,function(a){
  animalTable$Group[animalTable$AnimalID==a]})
GSIS$Sex<-factor(GSIS$Sex)
GSIS$Date<-as.POSIXct(strptime(GSIS$Date,"%Y/%m/%d",tz="America/Vancouver"))
GSIS$Group<-factor(GSIS$Group,levels = levels(weekly_monitoring$Group))
GSIS$FullGroup<-paste(GSIS$Sex,"_",GSIS$Group,sep="")
GSIS$FullGroup<-factor(GSIS$FullGroup,levels = levels(weekly_monitoring$FullGroup))
GSIS$Days <- sapply(1:nrow(GSIS),function(i){
  surgDate<-as.POSIXct(strptime("2018/02/20","%Y/%m/%d",tz="America/Vancouver"))
  d<-GSIS$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
GSIS$Weeks<-round(GSIS$Days/7)
GSIS<-GSIS[order(GSIS$Date,
                 GSIS$Sex,GSIS$Group,GSIS$AnimalID),]
GSIS$AnimalID <- factor(GSIS$AnimalID, levels = as.character(unique(GSIS$AnimalID)))
GSIS$CP.End[GSIS$AnimalID=="TJKV06M.24"&GSIS$Time==30&GSIS$Weeks==21]<-0.0083330
GSIS$CP<-ifelse(GSIS$CP.End==GSIS$CP.Start,GSIS$CP.End,NA)
GSIS$BG<-ifelse(GSIS$BG.End==GSIS$BG.Start,GSIS$BG.End,NA)

ITT<-GSIS[GSIS$Delivery=="ITT",]
GSIS<-GSIS[GSIS$Delivery!="ITT",]

GSIS$Delivery<-NULL
ITT$Delivery<-NULL


save.image("./preGSISAmelia.RData")
pr <- matrix(c(0,13,-10,log(0.0083330*2),0.9,
               0,14,33.3,147,0.995),byrow=T, nrow=2,ncol = 5)

# pr <- matrix(c(0,14,33.3,147,0.999),byrow=T, nrow=1,ncol = 5)
# bounds <- matrix(c(13,-10000,log(0.0083330*2)),ncol=3)

set.seed(12345)
GSIS.imp <- amelia(GSIS, m = 38, p2s=1, idvars = c("Date","Days",
                                                   "CP.End",
                                                   "CP.Start",
                                                   "BG.End",
                                                   "BG.Start",
                                                   "Sex","Group"),
                   ts="Time",cs="AnimalID",intercs = F,
                   ords="Weeks",polytime=1,
                   logs="CP",
                   noms = c("FullGroup"),
                   priors=pr,
                   #bounds=bounds,
                   multicore=8,
                   empri=0.001*nrow(GSIS)) 

plot(GSIS.imp,which.vars = 13:14)

head(GSIS.imp$imputations[[1]])
lapply(GSIS.imp$imputations, function(im){
  max(im$CP)
})

#### GSIS: Statistical Analyses ####
summary(GSIS.imp)
GSIS.imp <- transform(GSIS.imp, WeeksNom = factor(Weeks))
GSIS.imp <- transform(GSIS.imp, TimeNom = factor(Time))

save(GSIS,GSIS.imp,imps_IP,imps_IP_CP,file = "./prebrms.RData")

###### GSIS: Blood Glucose Stats ######
imps_IP<-lapply(GSIS.imp$imputations, function(im){
  im
}) #list of dataframes

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)

save.image("./prebrms.RData")


fitA <- brm_multiple(log(BG)~FullGroup*(WeeksNom/TimeNom)+(1|AnimalID),
                     iter=1000,
                     warmup=200,
                     family="gaussian",
                     prior=c(set_prior("normal(0,5)", class = "b")),
                     control = list(adapt_delta=0.5),
                     data = imps_IP[1:2],
                     chains=1,
                     future=T,
                     silent=0,
                     save_pars = save_pars(all = TRUE))

fitB <- brm_multiple(log(BG)~FullGroup*(WeeksNom/TimeNom)+(TimeNom|AnimalID),
                     iter=1000,
                     warmup=200,
                     family="gaussian",
                     prior=c(set_prior("normal(0,5)", class = "b")),
                     control = list(adapt_delta=0.5),
                     data = imps_IP[1:2],
                     chains=1,
                     future=T,
                     silent=0,
                     save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo") #much better
fitB <- add_criterion(fitB, "loo")

loo_compare(fitA,fitB) #gamma better than skew normal, with random effects is better

pp_check(fitA)
pp_check(fitB)

fit_GSIS_BG <- brm_multiple(BG~FullGroup*(WeeksNom/TimeNom)+
                              (WeeksNom|AnimalID),
                            iter=8000,
                            warmup=1000,
                            thin=4,
                            chains=8,
                            family="gamma",
                            prior=set_prior("normal(-0.5,5)", class = "b"),
                            control = list(max_treedepth = 12,
                                           adapt_delta=0.9),
                            data = imps_IP,
                            future=T,
                            silent=0)

# looking for ~similar plots of observed and predicted values
load("./fit_GSIS_BG.RData")
pp_check(fit_GSIS_BG)
GSIS_BG_summary
max(fit_GSIS_BG$rhats)

newdata2 = data.frame(Group = factor(rep(levels(GSIS$Group),44)|>str_sort()),
                      Sex=factor(c(rep(levels(weekly_monitoring$Sex),each=3)|>rep(4),
                                   rep(levels(weekly_monitoring$Sex),each=5)|>rep(2))|>rep(3)),
                      WeeksNom = c(rep(c(4,8,10,12),each=6),rep(c(16,21),each=10))|>rep(3),
                      TimeNom = rep(c(rep(c(0,30,60),8),rep(c(0,30,60,90,120),4)),3),
                      Weeks=c(rep(c(4,8,10,12),each=6),rep(c(16,21),each=10))|>rep(3))
newdata2$FullGroup<-paste(newdata2$Sex, newdata2$Group,sep="_")

fGSIS_BG <- fitted(fit_GSIS_BG,
               newdata = newdata2, 
               re_formula = NA) # extract the full MCMC

ff_GSIS_BG <- fGSIS_BG |>
  as_tibble() |>
  bind_cols(newdata2)

ff_GSIS_BG$Weeks<-c(rep(c(4,8,10,12),each=6),rep(c(16,21),each=10))|>rep(3)
ff_GSIS_BG$Time<-rep(c(rep(c(0,30,60),8),rep(c(0,30,60,90,120),4)),3)
ff_GSIS_BG$Group<-factor(ff_GSIS_BG$Group,levels=levels(GSIS$Group))

write_csv(ff_GSIS_BG,file="./GSIS_BG.csv")

gridtab<-matrix(c(2,3,4,5,6,7,8,9,10,11,12,13,14,7,9,11,13,15,17,19,21,23,25,27,29,31),ncol=2)
limits<-c(-5,125)
pd <- position_dodge(width=2)
breaks<-c(0,30,60,90,120)

p_GSIS_BG<-ggplot(GSIS,aes(x=Time,y=BG.Start,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~Weeks,scales = "fixed")+
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),size=0.3,linetype=2,alpha=0.8)+
  geom_hline(data=data.frame(yint = 33.3),aes(yintercept=yint),colour="#666666",linetype=3)+
  labs(x="Time (minutes)",y="Blood Glucose (mM)")+
  #geom_hline(data=data.frame(yint = 33.3,Species="Mouse"),aes(yintercept=yint),colour="#666666",linetype=3)+
  coord_cartesian(ylim=c(0, 35))+
  scale_x_continuous(breaks = breaks)+
  ggtitle("Weeks post implant")+
  theme(plot.title = element_text(family = "Arial",color="black", size=8, hjust=0.5))+
  #facet_grid(~Weeks,scales="free_x",space="free_x")+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  # guides(fill=guide_legend(ncol=3),
  #        colour=guide_legend(ncol=3),
  #        shape=guide_legend(ncol=3))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank())+ 
# panel.spacing.y = unit(2,"cm"),
# legend.position = c(0.5,0.5)
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_smooth(data = ff_GSIS_BG,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1) +
  geom_point(data=ff_GSIS_BG,
             aes(y=Estimate,shape=FullGroup),position=pd,size=2,alpha=0.7)
p_GSIS_BG

n<-length(unique(GSIS$Weeks))
nn<-gridtab[n-1,2]
z <- ggplotGrob(p_GSIS_BG)
z <- gtable_add_padding(z, unit(0.5, "cm"))
z <- gtable_add_rows(z, unit(z$heights[[3]], 'cm'), 2)
z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, fill = gray(0.85)),height=unit(2.5,"npc")),
                          textGrob("Weeks post implant", gp = gpar(col = gray(0),
                                                                   fontfamily="Arial",fontsize=8))),
                     2, 6, 1,nn+1, name = paste(runif(2))) #4n+1
z <- gtable_add_rows(z, unit(2/8, "line"), 3)
grid.newpage()
grid.draw(z)

ggsave("BG_IP_Bayes.png",plot=z,
       path="./Figures/GSIS",width = 30, height = 14, units = "cm")

#### Comparisons ####
GSISBGtab1<-rbind(
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC)-
             exp(Intercept)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)-
             exp(Intercept+WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)-
             exp(Intercept+WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)-
             exp(Intercept+WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)-
             exp(Intercept+WeeksNom8+WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)-
             exp(Intercept+WeeksNom8+WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)-
             exp(Intercept+WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)-
             exp(Intercept+WeeksNom10+WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)-
             exp(Intercept+WeeksNom10+WeeksNom10:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)-
             exp(Intercept+WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)-
             exp(Intercept+WeeksNom12+WeeksNom12:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)-
             exp(Intercept+WeeksNom12+WeeksNom12:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)-
             exp(Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)-
             exp(Intercept+WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC)-
             exp(Intercept)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             exp(Intercept+WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             exp(Intercept+WeeksNom4:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)-
             exp(Intercept+WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             exp(Intercept+WeeksNom8+WeeksNom8:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             exp(Intercept+WeeksNom8+WeeksNom8:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)-
             exp(Intercept+WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             exp(Intercept+WeeksNom10+WeeksNom10:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             exp(Intercept+WeeksNom10+WeeksNom10:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)-
             exp(Intercept+WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             exp(Intercept+WeeksNom12+WeeksNom12:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             exp(Intercept+WeeksNom12+WeeksNom12:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)-
             exp(Intercept+WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             exp(Intercept+WeeksNom16+WeeksNom16:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)-
             exp(Intercept+WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             exp(Intercept+WeeksNom21+WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC)-
             exp(Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC)-
             exp(Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC)-
             exp(Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC)-
             exp(Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)-
             exp(Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)-
             exp(Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)-
             exp(Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)<0")$hypothesis)

#Comparisons to time 0 within the week
GSISBGtab2<-rbind(
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom4:TimeNom30)-
             exp(Intercept)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom4:TimeNom60)-
             exp(Intercept)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+WeeksNom8:TimeNom30)-
             exp(Intercept+WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+WeeksNom8:TimeNom60)-
             exp(Intercept+WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+WeeksNom10:TimeNom30)-
             exp(Intercept+WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+WeeksNom10:TimeNom60)-
             exp(Intercept+WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+WeeksNom12:TimeNom30)-
             exp(Intercept+WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+WeeksNom12:TimeNom60)-
             exp(Intercept+WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+WeeksNom16:TimeNom30)-
             exp(Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+WeeksNom16:TimeNom60)-
             exp(Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+WeeksNom16:TimeNom90)-
             exp(Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+WeeksNom16:TimeNom120)-
             exp(Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+WeeksNom21:TimeNom30)-
             exp(Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+WeeksNom21:TimeNom60)-
             exp(Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+WeeksNom21:TimeNom90)-
             exp(Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+WeeksNom21:TimeNom120)-
             exp(Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupFemale_KC+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupFemale_KC+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupFemale_KC+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupFemale_KC+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupFemale_KC+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupFemale_KC+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupFemale_TC+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupFemale_TC+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupFemale_TC+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupFemale_TC+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupFemale_TC+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupFemale_TC+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupMale_FP+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupMale_FP+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupMale_FP+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupMale_FP+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupMale_FP+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupMale_FP+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupMale_KC+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupMale_KC+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupMale_KC+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupMale_KC+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupMale_KC+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupMale_KC+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)-
             exp(Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)-
             exp(Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupMale_TC+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)-
             exp(Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom8+FullGroupMale_TC+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)-
             exp(Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupMale_TC+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)-
             exp(Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom10+FullGroupMale_TC+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)-
             exp(Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupMale_TC+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)-
             exp(Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom12+FullGroupMale_TC+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)-
             exp(Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)-
             exp(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)-
             exp(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)-
             exp(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)-
             exp(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)-
             exp(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)-
             exp(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)-
             exp(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_BG,"exp(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)-
             exp(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis)

###### GSIS: C-Peptide Stats ######
options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)

imps_IP_CP<-lapply(imps_IP, function(im){
  im<-im[!is.na(im$CP.Start),]
  im$CPnM<-im$CP*1000/3020.29
  im
})

fitA <- brm(bf(CPnM~FullGroup*(WeeksNom/TimeNom)+(1|AnimalID),
               sigma~FullGroup),
            iter=4000,
            warmup=1000,
            family="skew_normal",
            # control = list(adapt_delta=0.5),
            prior=c(set_prior("normal(0,3)", class = "b"),
                    set_prior("normal(0,3)", class = "Intercept")),
            data = imps_IP_CP[[1]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitB <-brm(bf(CPnM~FullGroup*(WeeksNom/TimeNom)+(WeeksNom|AnimalID),
              sigma~FullGroup),
           iter=4000,
           warmup=1000,
           family="skew_normal",
           prior=c(set_prior("normal(0,5)", class = "b"),
                    set_prior("normal(0,5)", class = "Intercept")),
           control = list(adapt_delta=0.9,
                          max_treedepth = 15),
           data = imps_IP_CP[[1]],
           silent=0,
           backend = "cmdstanr",
           threads = threading(parallel::detectCores()),
           save_pars = save_pars(all = TRUE))

bayestestR::check_prior(fitB)

fitA <- add_criterion(fitA, "loo") 
fitB <- add_criterion(fitB, "loo")

loo_compare(fitA,fitB,criterion = "loo")

pp_check(fitB)+scale_y_log10()

ff_GSIS_CP<-bind_rows(ff_GSIS_CP_early,ff_GSIS_CP_late)|>arrange(Group,Sex,WeeksNom,TimeNom)

ff_GSIS_CP$Weeks<-c(rep(c(4,8,10,12),each=3),rep(c(16,21),each=5))|>rep(3*2)
ff_GSIS_CP$Time<-rep(c(rep(c(0,30,60),4),rep(c(0,30,60,90,120),2)),3*2)
ff_GSIS_CP$Group<-factor(ff_GSIS_CP$Group,levels=levels(GSIS$Group))

GSIS_IP_CP_sum<-summary(fit_GSIS_IP_CP)

fit_GSIS_IP_CP <- brm_multiple(bf(CPnM ~ FullGroup * (WeeksNom/TimeNom) + (1 | AnimalID), 
                                  sigma ~ FullGroup),
                               iter=8000,
                               warmup=1000,
                               thin=5,
                               family = "skew_normal",
                               chains=8,
                               prior=c(set_prior("normal(0,5)", class = "b"),
                                       set_prior("normal(0,5)", class = "Intercept")),
                               control = list(max_treedepth = 15,
                                              adapt_delta=0.9),
                               data = imps_IP_CP,
                               silent=0,
                               backend = "cmdstanr",
                               threads = threading(parallel::detectCores()),
                               save_pars = save_pars(all = TRUE))


save(GSIS_IP_CP_sum,fit_GSIS_IP_CP,file="./fit_GSIS_IP_CP.RData")

load("./fit_GSIS_IP_CP2.RData")

GSIS_IP_CP_sum<-summary(fit_GSIS_IP_CP)
summary(fit_GSIS_IP_CP)
pp_check(fit_GSIS_IP_CP)  # ~similar plots of observed and predicted values

max(fit_GSIS_IP_CP$rhats)

# already built newdata2 in Blood Glucose section
newdata2 = data.frame(Group = factor(rep(levels(GSIS$Group),44)|>str_sort()),
                      Sex=factor(c(rep(levels(weekly_monitoring$Sex),each=3)|>rep(4),
                                   rep(levels(weekly_monitoring$Sex),each=5)|>rep(2))|>rep(3)),
                      WeeksNom = factor(c(rep(c(4,8,10,12),each=6),rep(c(16,21),each=10))|>rep(3)),
                      TimeNom = factor(rep(c(rep(c(0,30,60),8),rep(c(0,30,60,90,120),4)),3)))
newdata2$FullGroup<-paste(newdata2$Sex, newdata2$Group,sep="_")

IP_CP_mod <- fitted(fit_GSIS_IP_CP,
               newdata = newdata2, 
               re_formula = NA) # extract the full MCMC

ff_GSIS_CP <- IP_CP_mod |>
  as_tibble() |>
  bind_cols(newdata2)

ff_GSIS_CP$Weeks<-c(rep(c(4,8,10,12),each=6),rep(c(16,21),each=10))|>rep(3)
ff_GSIS_CP$Time<-rep(c(rep(c(0,30,60),8),rep(c(0,30,60,90,120),4)),3)
ff_GSIS_CP$Group<-factor(ff_GSIS_CP$Group,levels=levels(GSIS$Group))

write_csv(ff_GSIS_CP,file="./GSIS_CP.csv")

gridtab<-matrix(c(2,3,4,5,6,7,8,9,10,11,12,13,14,7,9,11,13,15,17,19,21,23,25,27,29,31),ncol=2)
pd <- position_dodge(width=2)
limits<-c(-5,125)
breaks<-c(0,30,60,90,120)

p_GSIS_CP<-ggplot(GSIS,aes(x=Time,y=CP.End*1000/3020.29,colour=FullGroup,
                   group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~Weeks,scales = "fixed")+
  geom_linerange(aes(ymin=CP.Start*1000/3020.29,ymax=CP.End*1000/3020.29,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),size=0.3,linetype=2,alpha=0.8)+
  labs(x="Time (minutes)",y="Human C-peptide (nM)")+
  geom_hline(yintercept = 0.00833*1000/3020.29,linetype=3,alpha=0.7)+
  scale_x_continuous(limits=limits,breaks = breaks)+ #
  scale_y_continuous(breaks=pretty_breaks(6))+
  ggtitle("Weeks post implant")+
  theme(plot.title = element_text(family = "Arial",color="black", size=8, hjust=0.5))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_smooth(data = ff_GSIS_CP, 
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup,colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1)+
  geom_point(data=ff_GSIS_CP,
             aes(y=Estimate,shape=FullGroup),position=pd,size=2,alpha=0.8)
p_GSIS_CP

# n<-length(unique(GSIS$Weeks))
# nn<-gridtab[n-1,2]
# p<-p+facet_grid(Sex~Weeks,scales = "fixed")
# z <- ggplotGrob(p)
# z <- gtable_add_padding(z, unit(0.5, "cm"))
# z <- gtable_add_rows(z, unit(z$heights[[3]], 'cm'), 2)
# z <- gtable_add_grob(z, 
#                      list(rectGrob(gp = gpar(col = NA, fill = gray(0.85)),height=unit(2.5,"npc")),
#                           textGrob("Weeks post implant", gp = gpar(col = gray(0),
#                                                                    fontfamily="Arial",fontsize=8))),
#                      2, 6, 1,nn+1, name = paste(runif(2))) #4n+1 #t, l, b, r
# z <- gtable_add_rows(z, unit(2/8, "line"), 3)
# grid.newpage()
# grid.draw(z)

ggsave("CP_IP_Bayes.png", plot=z,
       path="./Figures/GSIS",width = 30, height = 14, units = "cm")

#### Comparisons ####
GSISCPtab1<-rbind(
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)-
             (Intercept+WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)-
             (Intercept+WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)-
             (Intercept+WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)-
             (Intercept+WeeksNom8+WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)-
             (Intercept+WeeksNom8+WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)-
             (Intercept+WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)-
             (Intercept+WeeksNom10+WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)-
             (Intercept+WeeksNom10+WeeksNom10:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)-
             (Intercept+WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)-
             (Intercept+WeeksNom12+WeeksNom12:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)-
             (Intercept+WeeksNom12+WeeksNom12:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)-
             (Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)-
             (Intercept+WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             (Intercept+WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             (Intercept+WeeksNom4:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)-
             (Intercept+WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             (Intercept+WeeksNom8+WeeksNom8:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             (Intercept+WeeksNom8+WeeksNom8:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)-
             (Intercept+WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             (Intercept+WeeksNom10+WeeksNom10:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             (Intercept+WeeksNom10+WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)-
             (Intercept+WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             (Intercept+WeeksNom12+WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             (Intercept+WeeksNom12+WeeksNom12:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)-
             (Intercept+WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             (Intercept+WeeksNom16+WeeksNom16:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)-
             (Intercept+WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             (Intercept+WeeksNom21+WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)-
             (Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)-
             (Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)-
             (Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)<0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)-
             (Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)-
             (Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)-
             (Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)-
             (Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)-
             (Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)-
             (Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)-
             (Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)-
             (Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)-
             (Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)-
             (Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)-
             (Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)-
             (Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)-
             (Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)-
             (Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)-
             (Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)=0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)=0")$hypothesis)

#Comparisons to time 0 within the week
GSISCPtab2<-rbind(
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom4:TimeNom30)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom4:TimeNom60)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+WeeksNom8:TimeNom30)-
             (Intercept+WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+WeeksNom8:TimeNom60)-
             (Intercept+WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+WeeksNom10:TimeNom30)-
             (Intercept+WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+WeeksNom10:TimeNom60)-
             (Intercept+WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+WeeksNom12:TimeNom30)-
             (Intercept+WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+WeeksNom12:TimeNom60)-
             (Intercept+WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom30)-
             (Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom60)-
             (Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom90)-
             (Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom120)-
             (Intercept+WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom30)-
             (Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom60)-
             (Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom90)-
             (Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom120)-
             (Intercept+WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupFemale_KC+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupFemale_KC+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupFemale_KC+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupFemale_KC+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupFemale_KC+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupFemale_KC+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_KC+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_KC+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupFemale_TC+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupFemale_TC+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupFemale_TC+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupFemale_TC+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupFemale_TC+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupFemale_TC+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupFemale_TC+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupFemale_TC+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupMale_FP+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupMale_FP+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom8+FullGroupMale_FP:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupMale_FP+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupMale_FP+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom10+FullGroupMale_FP:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupMale_FP+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupMale_FP+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom12+FullGroupMale_FP:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_FP+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom16+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_FP+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom21+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupMale_KC+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupMale_KC+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupMale_KC+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupMale_KC+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupMale_KC+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupMale_KC+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_KC+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_KC+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupMale_TC+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+FullGroupMale_TC+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupMale_TC+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+FullGroupMale_TC+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupMale_TC+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+FullGroupMale_TC+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+FullGroupMale_TC+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+FullGroupMale_TC+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis)

#Sex differences (absolute differences)
GSISCPtab3<-rbind(
  hypothesis(fit_GSIS_IP_CP,"(Intercept)-(Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom4:TimeNom30+FullGroupMale_FP:WeeksNom4:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom4:TimeNom60+FullGroupMale_FP:WeeksNom4:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8)-(Intercept+WeeksNom8+FullGroupMale_FP+FullGroupMale_FP:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom8+WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_FP:WeeksNom8+FullGroupMale_FP:WeeksNom8:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10)-(Intercept+WeeksNom10+FullGroupMale_FP+FullGroupMale_FP:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom10+WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_FP:WeeksNom10+FullGroupMale_FP:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12)-(Intercept+WeeksNom12+FullGroupMale_FP+FullGroupMale_FP:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom12+WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_FP:WeeksNom12+FullGroupMale_FP:WeeksNom12:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16)-(Intercept+WeeksNom16+FullGroupMale_FP+FullGroupMale_FP:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom16+WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_FP:WeeksNom16+FullGroupMale_FP:WeeksNom16:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21)-(Intercept+WeeksNom21+FullGroupMale_FP+FullGroupMale_FP:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+WeeksNom21+WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_FP+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_FP:WeeksNom21+FullGroupMale_FP:WeeksNom21:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC)-(Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom30+FullGroupFemale_KC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom4:TimeNom30+FullGroupMale_KC:WeeksNom4:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom4:TimeNom60+FullGroupFemale_KC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom4:TimeNom60+FullGroupMale_KC:WeeksNom4:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom8+FullGroupFemale_KC:WeeksNom8)-
             (Intercept+FullGroupMale_KC+WeeksNom8+FullGroupMale_KC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_KC:WeeksNom8+FullGroupFemale_KC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_KC:WeeksNom8+FullGroupMale_KC:WeeksNom8:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom10+FullGroupFemale_KC:WeeksNom10)-
             (Intercept+FullGroupMale_KC+WeeksNom10+FullGroupMale_KC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_KC:WeeksNom10+FullGroupFemale_KC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_KC:WeeksNom10+FullGroupMale_KC:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+FullGroupFemale_KC:WeeksNom12)-
             (Intercept+FullGroupMale_KC+WeeksNom12+FullGroupMale_KC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_KC:WeeksNom12+FullGroupFemale_KC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_KC:WeeksNom12+FullGroupMale_KC:WeeksNom12:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+FullGroupFemale_KC:WeeksNom16)-
             (Intercept+FullGroupMale_KC+WeeksNom16+FullGroupMale_KC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_KC:WeeksNom16+FullGroupFemale_KC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_KC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_KC:WeeksNom16+FullGroupMale_KC:WeeksNom16:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+FullGroupFemale_KC:WeeksNom21)-
             (Intercept+FullGroupMale_KC+WeeksNom21+FullGroupMale_KC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom90)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_KC:WeeksNom21+FullGroupFemale_KC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_KC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_KC:WeeksNom21+FullGroupMale_KC:WeeksNom21:TimeNom120)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC)-(Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom30+FullGroupFemale_TC:WeeksNom4:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom4:TimeNom30+FullGroupMale_TC:WeeksNom4:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom4:TimeNom60+FullGroupFemale_TC:WeeksNom4:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom4:TimeNom60+FullGroupMale_TC:WeeksNom4:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+FullGroupFemale_TC:WeeksNom8)-
             (Intercept+FullGroupMale_TC+WeeksNom8+FullGroupMale_TC:WeeksNom8)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom30+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupFemale_TC:WeeksNom8+FullGroupFemale_TC:WeeksNom8:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom8+WeeksNom8:TimeNom60+FullGroupMale_TC:WeeksNom8+FullGroupMale_TC:WeeksNom8:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+FullGroupFemale_TC:WeeksNom10)-
             (Intercept+FullGroupMale_TC+WeeksNom10+FullGroupMale_TC:WeeksNom10)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom30+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupFemale_TC:WeeksNom10+FullGroupFemale_TC:WeeksNom10:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom10+WeeksNom10:TimeNom60+FullGroupMale_TC:WeeksNom10+FullGroupMale_TC:WeeksNom10:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+FullGroupFemale_TC:WeeksNom12)-
             (Intercept+FullGroupMale_TC+WeeksNom12+FullGroupMale_TC:WeeksNom12)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom30+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupFemale_TC:WeeksNom12+FullGroupFemale_TC:WeeksNom12:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom12+WeeksNom12:TimeNom60+FullGroupMale_TC:WeeksNom12+FullGroupMale_TC:WeeksNom12:TimeNom60)<0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+FullGroupFemale_TC:WeeksNom16)-
             (Intercept+FullGroupMale_TC+WeeksNom16+FullGroupMale_TC:WeeksNom16)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom30+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom60+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom90)-
             (Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom90+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupFemale_TC:WeeksNom16+FullGroupFemale_TC:WeeksNom16:TimeNom120)-
             (Intercept+FullGroupMale_TC+WeeksNom16+WeeksNom16:TimeNom120+FullGroupMale_TC:WeeksNom16+FullGroupMale_TC:WeeksNom16:TimeNom120)<0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+FullGroupFemale_TC:WeeksNom21)-
             (Intercept+FullGroupMale_TC+WeeksNom21+FullGroupMale_TC:WeeksNom21)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom30)-
             (Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom30+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom30)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom60)-
             (Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom60+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom60)>0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom90)-
             (Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom90+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom90)<0")$hypothesis,
  hypothesis(fit_GSIS_IP_CP,"(Intercept+FullGroupFemale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupFemale_TC:WeeksNom21+FullGroupFemale_TC:WeeksNom21:TimeNom120)-
             (Intercept+FullGroupMale_TC+WeeksNom21+WeeksNom21:TimeNom120+FullGroupMale_TC:WeeksNom21+FullGroupMale_TC:WeeksNom21:TimeNom120)<0")$hypothesis)

##### ITT: Imputations #####
#ITT defined in GSIS section
# ITT$FullGroup<-paste(ITT$Sex,"_",ITT$Group,sep="")
# ITT$FullGroup<-factor(ITT$FullGroup,levels = levels(weekly_monitoring$FullGroup))
ITT<-ITT[!is.na(ITT$BG.Start),]
ITT$Sex<-sapply(ITT$AnimalID,function(a){
  animalTable$Sex[animalTable$AnimalID==a]})

#drop rescued animal
ITT<-filter(ITT,!(AnimalID=="TJKV06M.03"&Time>60))

CProws<-which(ITT$Time%in%c(0,30,60))
pr<-matrix(c(CProws,rep(c(13,-10,log(0.008333),0.9),each=length(CProws))),ncol=5)
pr<-rbind(pr,c(0,14,0,1.1,0.95))

#bounds<-matrix(c(14,0,1.1*1.5),ncol=3,byrow = T) #Blood glucose with a little wiggle room (mM)

save.image("./preAmelia.RData")

set.seed(12345)
ITT.imp <- amelia(ITT, m = 60, p2s=1, idvars = c("Date","Days","Weeks",
                                                    "CP.End",
                                                    "CP.Start",
                                                    "BG.End","BG.Start",
                                                    "Sex","Group"),
                     ts="Time",
                     cs="AnimalID",intercs = F,polytime=1,
                     logs=c("CP"),
                     noms = c("FullGroup"),
                     priors=pr,
                     multicore=8)

plot(ITT.imp,which.vars = 13:14)

lapply(ITT.imp$imputations, function(im){
  im<-filter(im,Time%in%c(0,30,60))
  im[which.max(im$CP),]
})

#### ITT: Statistical Analyses ####
summary(ITT.imp)
ITT.imp <- transform(ITT.imp, TimeNom = factor(Time))

imps_ITT<-lapply(ITT.imp$imputations, function(im){
  im$CPnM<-im$CP*1000/3020.29
  im
}) #list of dataframes

save(imps_ITT,imps_ITT_CP,ITT,ITT.imp,file="./prebrms_ITT.RData")

###### ITT: Blood Glucose Stats ######
options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)

fitA <- brm(BG~FullGroup*TimeNom+(1|AnimalID), 
                     #iter=1000,
                     #warmup=200,
                     thin=3,
                     family="skew_normal",
                     #prior=c(set_prior("normal(0,5)", class = "b")),
                     #control = list(adapt_delta=0.5),
                     data = imps_ITT[[1]],
                     #chains=1,
                     #future=T,
                     silent=0,
                     backend = "cmdstanr",
                     threads = threading(parallel::detectCores()),
                     save_pars = save_pars(all = TRUE))

fitB <- brm(BG~FullGroup*TimeNom+(TimeNom|AnimalID),
                     #iter=1000,
                     #warmup=200,
                     thin=3,
                     family="skew_normal",
                     #prior=c(set_prior("normal(0,5)", class = "b")),
                     #control = list(adapt_delta=0.5),
                     data = imps_ITT[[1]],
                     #chains=1,
                     silent=0,
                     backend = "cmdstanr",
                     threads = threading(parallel::detectCores()),
                     save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo")
fitB <- add_criterion(fitB, "loo")

loo_compare(fitA,fitB,criterion = "loo")

pp_check(fitA)
pp_check(fitB)

fit_ITT_BG <- brm_multiple(BG~FullGroup*TimeNom+(1|AnimalID),
                           data = imps_ITT,
                           iter=8000,
                           warmup=1000,
                           thin=5,
                           family="skew_normal",
                           prior=c(set_prior("normal(10,5)", class = "Intercept"),
                                   set_prior("normal(0,5)", class = "b")),
                           control = list(max_treedepth = 12,
                                          adapt_delta=0.99),
                           chains=4,
                           future=T)

ITT_BG_summary<-summary(fit_ITT_BG)

save(ITT_BG_summary,fit_ITT_BG,file="fit_ITT_BG.RData")

load("./fit_ITT_BG.RData")

pp_check(fit_ITT_BG) # ~similar plots of observed and predicted values
summary(fit_ITT_BG)

newdata3 = data.frame(Group = factor(rep(levels(ITT$Group),2*7)|>str_sort()),
                      Sex=factor(c(rep(levels(weekly_monitoring$Sex),7*3))),
                      TimeNom = factor(rep(c(0,10,20,30,60,90,120),each=2)|>rep(3)))
newdata3$FullGroup<-paste(newdata3$Sex, newdata3$Group,sep="_")

fitITT1 <- fitted(fit_ITT_BG,
                  newdata = newdata3, 
                  re_formula = NA) # extract the full MCMC

ff_ITT_BG <- fitITT1 |>
  as_tibble() |>
  bind_cols(newdata3)

ff_ITT_BG$Time<-rep(c(0,10,20,30,60,90,120),each=2)|>rep(3)
ff_ITT_BG$FullGroup<-factor(ff_ITT_BG$FullGroup,levels=levels(ITT$FullGroup))
ff_ITT_BG$Group<-factor(ff_ITT_BG$Group,levels=levels(ITT$Group))
ff_ITT_BG$Sex<-factor(ff_ITT_BG$Sex,levels=levels(ITT$Sex))

write_csv(ff_ITT_BG,file="./ITT_BG.csv")

limits<-c(-5,125)
breaks<-c(0,10,20,30,60,90,120)

p_ITT_BG<-ggplot(ITT,aes(x=Time,y=BG.Start,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "fixed")+
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,group=AnimalID),
                 linetype=3,size=0.3,alpha=0.8,show.legend = F)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),
            position=pd,linetype=2,size=0.3,alpha=0.8,show.legend = F)+
  geom_hline(data=data.frame(yint = 1.1),aes(yintercept=yint),colour="#666666",linetype=3)+
  geom_point(data = ff_ITT_BG,
             aes(y = Estimate, shape=FullGroup),size=2,alpha=0.7,show.legend = F)+
  geom_smooth(data = ff_ITT_BG,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1,show.legend = F) +
  labs(x="Time (minutes)",y="Blood Glucose (mM)")+
  coord_cartesian(ylim=c(0, 15))+
  scale_x_continuous(breaks = breaks)+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank(),
        legend.background = element_blank())+
        # ,
        # legend.position = c(0.5,0.85))+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p_ITT_BG

ggsave("BG_ITT_Bayes.png",path="./Figures/GSIS",width = 20, height = 12, units = "cm")

##### Comparisons #####
#Within sex, between sites
#Reference: FullGroupFemale_FP TimeNom0
ITTBGtab1<-rbind(
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC)-(Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom10+FullGroupFemale_KC:TimeNom10)-
             (Intercept+TimeNom10)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom20+FullGroupFemale_KC:TimeNom20)-
             (Intercept+TimeNom20)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
             (Intercept+TimeNom30)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+TimeNom60)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom90+FullGroupFemale_KC:TimeNom90)-
             (Intercept+TimeNom90)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom120+FullGroupFemale_KC:TimeNom120)-
             (Intercept+TimeNom120)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC)-(Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom10+FullGroupFemale_TC:TimeNom10)-
             (Intercept+TimeNom10)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom20+FullGroupFemale_TC:TimeNom20)-
             (Intercept+TimeNom20)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+TimeNom30)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+TimeNom60)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom90+FullGroupFemale_TC:TimeNom90)-
             (Intercept+TimeNom90)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom120+FullGroupFemale_TC:TimeNom120)-
             (Intercept+TimeNom120)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC)-(Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom10+FullGroupFemale_TC:TimeNom10)-
             (Intercept+FullGroupFemale_KC+TimeNom10+FullGroupFemale_KC:TimeNom10)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom20+FullGroupFemale_TC:TimeNom20)-
             (Intercept+FullGroupFemale_KC+TimeNom20+FullGroupFemale_KC:TimeNom20)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom90+FullGroupFemale_TC:TimeNom90)-
             (Intercept+FullGroupFemale_KC+TimeNom90+FullGroupFemale_KC:TimeNom90)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom120+FullGroupFemale_TC:TimeNom120)-
             (Intercept+FullGroupFemale_KC+TimeNom120+FullGroupFemale_KC:TimeNom120)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC)-(Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom10+FullGroupMale_KC:TimeNom10)-
             (Intercept+FullGroupMale_FP+TimeNom10+FullGroupMale_FP:TimeNom10)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom20+FullGroupMale_KC:TimeNom20)-
             (Intercept+FullGroupMale_FP+TimeNom20+FullGroupMale_FP:TimeNom20)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom90+FullGroupMale_KC:TimeNom90)-
             (Intercept+FullGroupMale_FP+TimeNom90+FullGroupMale_FP:TimeNom90)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom120+FullGroupMale_KC:TimeNom120)-
             (Intercept+FullGroupMale_FP+TimeNom120+FullGroupMale_FP:TimeNom120)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom10+FullGroupMale_TC:TimeNom10)-
             (Intercept+FullGroupMale_FP+TimeNom10+FullGroupMale_FP:TimeNom10)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom20+FullGroupMale_TC:TimeNom20)-
             (Intercept+FullGroupMale_FP+TimeNom20+FullGroupMale_FP:TimeNom20)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom90+FullGroupMale_TC:TimeNom90)-
             (Intercept+FullGroupMale_FP+TimeNom90+FullGroupMale_FP:TimeNom90)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom120+FullGroupMale_TC:TimeNom120)-
             (Intercept+FullGroupMale_FP+TimeNom120+FullGroupMale_FP:TimeNom120)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom10+FullGroupMale_TC:TimeNom10)-
             (Intercept+FullGroupMale_KC+TimeNom10+FullGroupMale_KC:TimeNom10)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom20+FullGroupMale_TC:TimeNom20)-
             (Intercept+FullGroupMale_KC+TimeNom20+FullGroupMale_KC:TimeNom20)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom90+FullGroupMale_TC:TimeNom90)-
             (Intercept+FullGroupMale_KC+TimeNom90+FullGroupMale_KC:TimeNom90)=0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom120+FullGroupMale_TC:TimeNom120)-
             (Intercept+FullGroupMale_KC+TimeNom120+FullGroupMale_KC:TimeNom120)=0")$hypothesis)

#compared to 0 within group
ITTBGtab2<-rbind(
  hypothesis(fit_ITT_BG,"(Intercept+TimeNom10)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+TimeNom20)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+TimeNom30)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+TimeNom60)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+TimeNom90)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+TimeNom120)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom10+FullGroupFemale_KC:TimeNom10)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom20+FullGroupFemale_KC:TimeNom20)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom90+FullGroupFemale_KC:TimeNom90)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_KC+TimeNom120+FullGroupFemale_KC:TimeNom120)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom10+FullGroupFemale_TC:TimeNom10)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom20+FullGroupFemale_TC:TimeNom20)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom90+FullGroupFemale_TC:TimeNom90)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupFemale_TC+TimeNom120+FullGroupFemale_TC:TimeNom120)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom10+FullGroupMale_FP:TimeNom10)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom20+FullGroupMale_FP:TimeNom20)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom90+FullGroupMale_FP:TimeNom90)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom120+FullGroupMale_FP:TimeNom120)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom10+FullGroupMale_KC:TimeNom10)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom20+FullGroupMale_KC:TimeNom20)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom90+FullGroupMale_KC:TimeNom90)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom120+FullGroupMale_KC:TimeNom120)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom10+FullGroupMale_TC:TimeNom10)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom20+FullGroupMale_TC:TimeNom20)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom90+FullGroupMale_TC:TimeNom90)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom120+FullGroupMale_TC:TimeNom120)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis)

#Sex differences (absolute values) - don't actually bother with this
ITTBGtab3<-rbind(
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)-
             (Intercept+TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_BG,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)>0")$hypothesis)

###### ITT: C-Peptide Stats ######
library("future")
options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)

imps_ITT_CP<-lapply(imps_ITT, function(im){
  im<-filter(im,Time%in%c(0,30,60))
  im
}) #list of dataframes


fitA <- brm(bf(CPnM~FullGroup*TimeNom+(TimeNom|AnimalID)), 
            iter=4000,
            warmup=1000,
            thin=1,
            family="shifted_lognormal",
            prior=c(set_prior("normal(1,5)", class = "Intercept"),
                    set_prior("normal(0,5)", class = "b")),
            control = list(adapt_delta=0.9999,
                           max_treedepth=15),
            data = imps_ITT_CP[[1]],
            #chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

mix <- mixture(gaussian, gaussian)
prior <- c(
  prior(normal(0, 0.5), Intercept, dpar = mu1),
  prior(normal(1, 2), Intercept, dpar = mu2),
  prior(normal(0,5), b, dpar=mu1),
  prior(normal(0,5), b, dpar=mu2)
  # prior(normal(0, log(5)), class = Intercept, dpar = sigma1),
  # prior(normal(0, 5),
  #       class = sd, group = Device,
  #       dpar = sigma1),
  # prior(normal(0, log(5)), class = Intercept, dpar = sigma2),
  # prior(normal(0, 5),
  #       class = sd, group = Device,
  #       dpar = sigma2)
)

fitB <- brm(bf(CPnM~FullGroup*TimeNom+(TimeNom|AnimalID)), 
            iter=4000,
            warmup=1000,
            thin=1,
            family="skew_normal",
            prior=c(set_prior("normal(1,5)", class = "Intercept"),
                    set_prior("normal(0,5)", class = "b")),
            control = list(adapt_delta=0.9999,
                           max_treedepth=15),
            data = imps_ITT_CP[[1]],
            #chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo")
fitB <- add_criterion(fitB, "loo")

loo_compare(fitA,fitB,criterion = "loo")

pp_check(fitA)
pp_check(fitB)

bayestestR::check_prior(fit_ITT_CP)

fit_ITT_CP <- brm_multiple(bf(CPnM~FullGroup*TimeNom+(TimeNom|AnimalID)), 
                           iter=4000,
                           warmup=1000,
                           thin=1,
                           family="gaussian",
                           prior=c(set_prior("normal(1,5)", class = "Intercept"),
                                   set_prior("normal(0,5)", class = "b")),
                           control = list(adapt_delta=0.95,
                                          max_treedepth=15),
                           data = imps_ITT_CP,
                           #chains=1,
                           #future=T,
                           silent=0,
                           backend = "cmdstanr",
                           threads = threading(parallel::detectCores()),
                           save_pars = save_pars(all = TRUE))


#save(fit_ITT_CP,file="fit_ITT_CP_log.RData")


#save(fit_ITT_CP,file="./fit_ITT_CP.RData")
load("./fit_ITT_CP.RData")
pp_check(fit_ITT_CP) # ~similar plots of observed and predicted values

ITT_CP_summary

fit_ITT_CP

max(fit_ITT_CP$rhats)

newdata3 = data.frame(Group = factor(rep(levels(ITT$Group),6)|>str_sort()),
                      Sex=factor(c(rep(levels(weekly_monitoring$Sex),3))),
                      TimeNom = factor(rep(c(0,30,60),each=2)|>rep(3)))
newdata3$FullGroup<-paste(newdata3$Sex, newdata3$Group,sep="_")

ITT_CP_mod <- fitted(fit_ITT_CP,
                  newdata = newdata3, 
                  re_formula = NA) # extract the full MCMC

ff_ITT_CP <- ITT_CP_mod |>
  as_tibble() |>
  bind_cols(newdata3)


ff_ITT_CP$Time<-rep(c(0,30,60),each=2)|>rep(3)

ff_ITT_CP$FullGroup<-factor(ff_ITT_CP$FullGroup,levels=levels(ITT$FullGroup))
ff_ITT_CP$Group<-factor(ff_ITT_CP$Group,levels=levels(ITT$Group))
ff_ITT_CP$Sex<-factor(ff_ITT_CP$Sex,levels=levels(ITT$Sex))

write_csv(ff_ITT_CP,file="./ITT_CP.csv")

breaks<-c(0,30,60);limits<-c(-5,65)

p_ITT_CP<-ggplot(filter(ITT,Time%in%breaks),
          aes(x=Time,y=CP.End*1000/3020.29,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "fixed")+
  geom_linerange(aes(ymin=CP.Start*1000/3020.29,ymax=CP.End*1000/3020.29,group=AnimalID),
                 linetype=3,size=0.3,alpha=0.8)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),size=0.3,linetype=2,alpha=0.8)+
  labs(x="Time (minutes)",y="Human C-peptide (nM)")+
  geom_hline(yintercept = 0.00833*1000/3020.29,linetype=3)+
  scale_x_continuous(breaks = breaks,limits=limits)+ #
  scale_y_continuous(breaks=pretty_breaks(6))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank(),
        legend.background = element_blank())+
  geom_smooth(data = ff_ITT_CP,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup,colour=FullGroup),
              stat = "identity",
              alpha = 1/4, size = 1)+
  geom_point(data = ff_ITT_CP,
             aes(y = Estimate,shape=FullGroup),
             position=pd,size=2,alpha=0.7)
p_ITT_CP
ggsave("CP_ITT_Bayes.png", path="./Figures/GSIS",width = 22, height = 14, units = "cm")


##### Comparisons #####
ITTCPtab1<-rbind(
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_KC)-(Intercept)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
             (Intercept+TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC)-(Intercept)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC)-(Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC)-(Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)>0")$hypothesis)

ITTCPtab2<-rbind(
  hypothesis(fit_ITT_CP,"(Intercept+TimeNom30)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+TimeNom60)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis)

#Sex differences (absolute values)
ITTCPtab3<-rbind(
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_FP)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)-
             (Intercept+TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)>0")$hypothesis,
  hypothesis(fit_ITT_CP,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)>0")$hypothesis)

#Percent change: comparing same timepoint between sexes
hypothesis(fit_ITT_CP,"(TimeNom30)/Intercept-
           (TimeNom30+FullGroupMale_FP:TimeNom30)/
           (Intercept+FullGroupMale_FP)<0")
hypothesis(fit_ITT_CP,"(TimeNom60)/Intercept-
           (TimeNom60+FullGroupMale_FP:TimeNom60)/
           (Intercept+FullGroupMale_FP)<0")

hypothesis(fit_ITT_CP,"(Intercept+TimeNom30+FullGroupFemale_TC+FullGroupFemale_TC:TimeNom30)/
           (Intercept+FullGroupFemale_TC)-
           (Intercept+TimeNom30+FullGroupMale_TC+FullGroupMale_TC:TimeNom30)/
           (Intercept+FullGroupMale_TC)<0")
hypothesis(fit_ITT_CP,"(Intercept+TimeNom60+FullGroupFemale_TC+FullGroupFemale_TC:TimeNom60)/
           (Intercept+FullGroupFemale_TC)-
           (Intercept+TimeNom60+FullGroupMale_TC+FullGroupMale_TC:TimeNom60)/
           (Intercept+FullGroupMale_TC)<0")


hypothesis(fit_ITT_CP,"(Intercept+TimeNom30+FullGroupFemale_KC+FullGroupFemale_KC:TimeNom30)/
           (Intercept+FullGroupFemale_KC)-
           (Intercept+TimeNom30+FullGroupMale_KC+FullGroupMale_KC:TimeNom30)/
           (Intercept+FullGroupMale_KC)<0")
hypothesis(fit_ITT_CP,"(Intercept+TimeNom60+FullGroupFemale_KC+FullGroupFemale_KC:TimeNom60)/
           (Intercept+FullGroupFemale_KC)-
           (Intercept+TimeNom60+FullGroupMale_KC+FullGroupMale_KC:TimeNom60)/
           (Intercept+FullGroupMale_KC)<0")

########## ArgTT: Imputations and Figures ##########
ArgTT$Sex<-sapply(ArgTT$AnimalID,function(a){
  animalTable$Sex[animalTable$AnimalID==a]})
ArgTT$Group<-sapply(ArgTT$AnimalID,function(a){
  animalTable$Group[animalTable$AnimalID==a]})
ArgTT$Sex<-factor(ArgTT$Sex)
ArgTT$Date<-as.POSIXct(strptime(ArgTT$Date,"%Y/%m/%d",tz="America/Vancouver"))
ArgTT$FullGroup<-paste(ArgTT$Sex,"_",ArgTT$Group,sep="")
ArgTT$FullGroup<-factor(ArgTT$FullGroup,levels = levels(weekly_monitoring$FullGroup))
ArgTT$Group<-factor(ArgTT$Group,levels = levels(weekly_monitoring$Group))
ArgTT$Days <- sapply(1:nrow(ArgTT),function(i){
  surgDate<-as.POSIXct(strptime("2018/02/20","%Y/%m/%d",tz="America/Vancouver"))
  d<-ArgTT$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
ArgTT$Weeks<-round(ArgTT$Days/7)
ArgTT<-ArgTT[order(ArgTT$Date,
                   ArgTT$Sex,ArgTT$Group,ArgTT$AnimalID),]
ArgTT$AnimalID <- factor(ArgTT$AnimalID, levels = as.character(unique(ArgTT$AnimalID)))

#drop rescued animals
ArgTT<-filter(ArgTT,!(AnimalID=="TJKV06M.26"&Time>60)) |> filter(!(AnimalID=="TJKV06M.36"&Time>60))

#Hormones
ArgTT$GLP1<-ifelse(ArgTT$GLP1.End==ArgTT$GLP1.Start,ArgTT$GLP1.End,NA)
ArgTT$Ins<-ifelse(ArgTT$Ins.End==ArgTT$Ins.Start,ArgTT$Ins.End,NA)
ArgTT$Gcg<-ifelse(ArgTT$Gcg.End==ArgTT$Gcg.Start,ArgTT$Gcg.End,NA)
ArgTT$BG<-ifelse(ArgTT$BG.End==ArgTT$BG.Start,ArgTT$BG.End,NA)

Hormrows<-which(is.na(ArgTT$GLP1.Start))
pr1<-matrix(c(Hormrows,rep(c(18,-5,log(2*40/6),0.9),each=length(Hormrows))),ncol=5) #GLP-1 LLOD from kit instructions (pg/mL)
pr2<-matrix(c(Hormrows,rep(c(19,-5,log(8*40/6/1000),0.95),each=length(Hormrows))),ncol=5) #insulin LLOD from kit instructions (ng/mL)
pr3<-matrix(c(Hormrows,rep(c(20,-5,log(27*40/6),0.9),each=length(Hormrows))),ncol=5) #Glucagon LLOD from kit instructions (pg/mL)
pr<-rbind(pr1,pr2,pr3,c(0,21,0,1.1,0.9)) 

save.image("./preAmelia.RData")

set.seed(12345)
ArgTT.imp <- amelia(ArgTT, m = 54, p2s=1, 
                    idvars = c("Date","Days","Delivery","Ins.End","Ins.Start",
                               "GLP1.End","GLP1.Start","Gcg.End","Gcg.Start",
                               "BG.End","BG.Start","Group","Sex","Weeks"),
                    ts="Time",
                    cs="AnimalID",intercs = F,polytime=1,
                    logs=c("Ins","Gcg","GLP1"),
                    noms = c("FullGroup"), 
                    priors = pr,
                    #bounds=bounds,
                    multicore=8,
                    empri = .001*nrow(ArgTT))

plot(ArgTT.imp,which.vars = 18:21)

#summary(ArgTT.imp.Hormones)

lapply(ArgTT.imp$imputations,function(im){
  im[c(which.max(im$GLP1),which.max(im$Gcg)),]
})

ArgTT.imp <- transform(ArgTT.imp, TimeNom = factor(Time))
imps_ArgTT_Horm<-lapply(ArgTT.imp$imputations, function(im){
  im<-filter(im,Time%in%c(0,7,15,60))
  im$GLP1pM<-im$GLP1*1000/3297.6
  im$InsnM<-im$Ins*1000/5808
  im$GcgpM<-im$Gcg*1000/3485
  im
})

imps_ArgTT<-lapply(ArgTT.imp$imputations, function(im){
  im
})

##### ArgTT: BG Stats #####
#Only one missing value

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)

save(ArgTT,ArgTT.imp,imps_ArgTT,imps_ArgTT_Horm,file="./prebrms_ArgTT.RData")

fitA <- brm(bf(BG~FullGroup*TimeNom+(1|AnimalID)),
                     data=imps_ArgTT[[1]],
                     # iter=1000,
                     # warmup=200,
                     #family="skew_normal",
                     prior=c(set_prior("normal(10,3)", class = "Intercept"),
                             set_prior("normal(0,10)", class = "b")),
                     control = list(adapt_delta=0.8),
                     chains=4,
                     silent=0,
                     backend = "cmdstanr",
                     threads = threading(parallel::detectCores()),
                     save_pars = save_pars(all = TRUE))

fitB <- brm(BG~FullGroup*TimeNom+(TimeNom|AnimalID),
            data=imps_ArgTT[[1]],
            # iter=1000,
            # warmup=200,
            #family="gamma",
            prior=c(set_prior("normal(10,3)", class = "Intercept"),
                    set_prior("normal(0,10)", class = "b")),
            control = list(adapt_delta=0.9),
            chains=4,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo") 
fitB <- add_criterion(fitB, "loo") #much better

loo_compare(fitA,fitB,criterion = "loo")

pp_check(fitA)
pp_check(fitB)

fit_ArgTT_BG <- brm_multiple(BG~FullGroup*TimeNom+(TimeNom|AnimalID),
                             data = imps_ArgTT,
                             iter=8000,
                             thin=5,
                             warmup=1000,
                             prior=c(set_prior("normal(10,3)", class = "Intercept"),
                                     set_prior("normal(0,5)", class = "b")),
                             #family = "skew_normal",
                             control = list(adapt_delta=0.95,
                                            max_treedepth = 12),
                             future=T)

save(fit_ArgTT_BG,file="./fit_ArgTT_BG.RData")

load("./fit_ArgTT_BG.RData")

pp_check(fit_ArgTT_BG) # ~similar plots of observed and predicted values

summary(fit_ArgTT_BG)

newdata4 = data.frame(Group = factor(rep(levels(ArgTT$Group),2*7)|>str_sort()),
                      Sex=factor(c(rep(levels(ArgTT$Sex),7*3))),
                      TimeNom = factor(rep(c(0,7,15,30,60,90,120),each=2)|>rep(3)))

newdata4$FullGroup<-paste(newdata4$Sex, newdata4$Group,sep="_")

ArgTT_BG_mod <- fitted(fit_ArgTT_BG,
               newdata = newdata4, 
               re_formula = NA) # extract the full MCMC

ffArgTT_BG <- ArgTT_BG_mod |>
  as_tibble() |>
  bind_cols(newdata4)

ffArgTT_BG$Time<-rep(c(0,7,15,30,60,90,120),each=2)|>rep(3)
ffArgTT_BG$FullGroup<-factor(ffArgTT_BG$FullGroup,levels=levels(weekly_monitoring$FullGroup))
ffArgTT_BG$Group<-factor(ffArgTT_BG$Group,levels=levels(weekly_monitoring$Group))
ffArgTT_BG$Sex<-factor(ffArgTT_BG$Sex,levels=levels(weekly_monitoring$Sex))

write_csv(ffArgTT_BG,file="./ArgTT_BG.csv")

limits<-c(-5,125)
pd <- position_dodge(width=2)
breaks<-c(0,7,15,30,60,90,120)

p_ArgTT_BG<-ggplot(ArgTT,aes(x=Time,y=BG.End,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "fixed")+
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8,show.legend = F)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),size=0.3,linetype=2,alpha=0.8,show.legend = F)+
  geom_hline(data=data.frame(yint = 1.1),aes(yintercept=yint),colour="#666666",linetype=3)+
  labs(x="Time (minutes)",y="Blood Glucose (mM)")+
  coord_cartesian(ylim=c(0, 14))+
  scale_x_continuous(breaks = breaks)+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(data = ffArgTT_BG,aes(y = Estimate,shape=FullGroup),
             position=pd,size=2,alpha=0.7,show.legend = F)+
  geom_smooth(data = ffArgTT_BG,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1,show.legend = F) 
p_ArgTT_BG
ggsave("BG_ArgTT_Bayes.png",path="./Figures/ArgTT",width = 20, height = 12, units = "cm")

##### Comparisons #####
ArgTTBGtab1<-rbind(hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC)-
                           (Intercept)<0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
                           (Intercept+TimeNom30)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom90+FullGroupFemale_KC:TimeNom90)-
                           (Intercept+TimeNom90)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom120+FullGroupFemale_KC:TimeNom120)-
                           (Intercept+TimeNom120)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC)-
                           (Intercept)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
                           (Intercept+TimeNom30)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom90+FullGroupFemale_TC:TimeNom90)-
                           (Intercept+TimeNom90)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom120+FullGroupFemale_TC:TimeNom120)-
                           (Intercept+TimeNom120)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC)-
                           (Intercept+FullGroupFemale_TC)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
                           (Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom90+FullGroupFemale_KC:TimeNom90)-
                           (Intercept+FullGroupFemale_TC+TimeNom90+FullGroupFemale_TC:TimeNom90)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom120+FullGroupFemale_KC:TimeNom120)-
                           (Intercept+FullGroupFemale_TC+TimeNom120+FullGroupFemale_TC:TimeNom120)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
                           (Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom90+FullGroupMale_KC:TimeNom90)-
                           (Intercept+FullGroupMale_FP+TimeNom90+FullGroupMale_FP:TimeNom90)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom120+FullGroupMale_KC:TimeNom120)-
                           (Intercept+FullGroupMale_FP+TimeNom120+FullGroupMale_FP:TimeNom120)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
                           (Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom90+FullGroupMale_TC:TimeNom90)-
                           (Intercept+FullGroupMale_FP+TimeNom90+FullGroupMale_FP:TimeNom90)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom120+FullGroupMale_TC:TimeNom120)-
                           (Intercept+FullGroupMale_FP+TimeNom120+FullGroupMale_FP:TimeNom120)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_KC)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
                           (Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)=0")$hypothesis,
                hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom90+FullGroupMale_TC:TimeNom90)-
                           (Intercept+FullGroupMale_KC+TimeNom90+FullGroupMale_KC:TimeNom90)=0")$hypothesis,
              hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom120+FullGroupMale_TC:TimeNom120)-
                           (Intercept+FullGroupMale_KC+TimeNom120+FullGroupMale_KC:TimeNom120)=0")$hypothesis)

ArgTTBGtab2<-rbind(         
  hypothesis(fit_ArgTT_BG,"(Intercept+TimeNom7)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+TimeNom15)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+TimeNom30)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+TimeNom60)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+TimeNom90)-
             (Intercept)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+TimeNom120)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom30+FullGroupFemale_KC:TimeNom30)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom90+FullGroupFemale_KC:TimeNom90)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_KC+TimeNom120+FullGroupFemale_KC:TimeNom120)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
             (Intercept+FullGroupFemale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
             (Intercept+FullGroupFemale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom30+FullGroupFemale_TC:TimeNom30)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom90+FullGroupFemale_TC:TimeNom90)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupFemale_TC+TimeNom120+FullGroupFemale_TC:TimeNom120)-
             (Intercept+FullGroupFemale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
             (Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_FP+TimeNom30+FullGroupMale_FP:TimeNom30)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_FP+TimeNom90+FullGroupMale_FP:TimeNom90)-
             (Intercept+FullGroupMale_FP)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_FP+TimeNom120+FullGroupMale_FP:TimeNom120)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
             (Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom30+FullGroupMale_KC:TimeNom30)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom90+FullGroupMale_KC:TimeNom90)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_KC+TimeNom120+FullGroupMale_KC:TimeNom120)-
             (Intercept+FullGroupMale_KC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
             (Intercept+FullGroupMale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom30+FullGroupMale_TC:TimeNom30)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom90+FullGroupMale_TC:TimeNom90)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis,
  hypothesis(fit_ArgTT_BG,"(Intercept+FullGroupMale_TC+TimeNom120+FullGroupMale_TC:TimeNom120)-
             (Intercept+FullGroupMale_TC)<0")$hypothesis)

##### ArgTT: GLP-1 Stats #####
# options(mc.cores = parallel::detectCores())
# plan(multisession)
# rstan::rstan_options(auto_write = TRUE)

set.seed(12345)

load("./prebrms_ArgTT.RData")

fitA <- brm(bf(GLP1pM~FullGroup*TimeNom+(1|AnimalID),
               sigma~FullGroup),
            iter=8000,
            warmup=2000,
            thin=5,
            family="gaussian",
            prior = c(set_prior("normal(0,5)", class = "b"),
                      set_prior("normal(10,5)", class = "Intercept")),
            control = list(adapt_delta=0.9),
            data = imps_ArgTT_Horm[[28]],
            chains=8,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

mix <- mixture(gaussian, gaussian)
prior <- c(
  prior(normal(0, 10), Intercept, dpar = mu1),
  prior(normal(50, 5), Intercept, dpar = mu2),
  prior(normal(0,5), b, dpar=mu1),
  prior(normal(0,5), b, dpar=mu2)
  # prior(normal(0, log(5)), class = Intercept, dpar = sigma1),
  # prior(normal(0, 5),
  #       class = sd, group = Device,
  #       dpar = sigma1),
  # prior(normal(0, log(5)), class = Intercept, dpar = sigma2),
  # prior(normal(0, 5),
  #       class = sd, group = Device,
  #       dpar = sigma2)
)

fitB <- brm(bf(GLP1pM~FullGroup*TimeNom+(TimeNom|AnimalID)),
            iter=8000,
            warmup=2000,
            thin=5,
            family="skew_normal",
            prior = c(set_prior("normal(0,5)", class = "b"),
                      set_prior("normal(10,5)", class = "Intercept")),
            control = list(adapt_delta=0.95),
            data = imps_ArgTT_Horm[[28]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo")
fitB <- add_criterion(fitB, "loo")

loo_compare(fitA,fitB,criterion = "loo")

pp_check(fitA)
pp_check(fitB)

bayestestR::check_prior(fit_ArgTT_GLP1)

#real deal
fit_ArgTT_GLP1 <- brm_multiple(bf(GLP1pM~FullGroup*TimeNom+(1|AnimalID),
                                  sigma~FullGroup),
                               iter=8000,
                               warmup=2000,
                               family="skew_normal",
                               thin=5,
                               prior = c(set_prior("normal(0,5)", class = "b"),
                                         set_prior("normal(10,5)", class = "Intercept")),
                               control = list(adapt_delta=0.95),
                               data = imps_ArgTT_Horm,
                               #chains=8,
                               silent=0,
                               backend = "cmdstanr",
                               threads = threading(parallel::detectCores()),
                               save_pars = save_pars(all = TRUE))

GLP1_summary<-summary(fit_ArgTT_GLP1)

#save(fit_ArgTT_GLP1,GLP1_summary,file="./fit_ArgTT_GLP1.RData")
load("./fit_ArgTT_GLP1.RData")

pp_check(fit_ArgTT_GLP1) # ~similar plots of observed and predicted values
fit_ArgTT_GLP1
max(fit_ArgTT_GLP1$rhats)

newdata4 = data.frame(Group = factor(rep(levels(ArgTT$Group),2*4)|>str_sort()),
                      Sex=factor(c(rep(levels(ArgTT$Sex),4*3))),
                      TimeNom = factor(rep(c(0,7,15,60),each=2)|>rep(3)))

newdata4$FullGroup<-paste(newdata4$Sex, newdata4$Group,sep="_")

ArgTT_GLP1_mod <- fitted(fit_ArgTT_GLP1,
               newdata = newdata4, 
               re_formula = NA) # extract the full MCMC

ff_ArgTT_GLP1 <- ArgTT_GLP1_mod |>
  as_tibble() |>
  bind_cols(newdata4)

ff_ArgTT_GLP1$Time<-rep(c(0,7,15,60),each=2)|>rep(3)
ff_ArgTT_GLP1$FullGroup<-factor(ff_ArgTT_GLP1$FullGroup,levels=levels(weekly_monitoring$FullGroup))
ff_ArgTT_GLP1$Group<-factor(ff_ArgTT_GLP1$Group,levels=levels(weekly_monitoring$Group))
ff_ArgTT_GLP1$Sex<-factor(ff_ArgTT_GLP1$Sex,levels=levels(weekly_monitoring$Sex))

write_csv(ff_ArgTT_GLP1,file="./ArgTT_GLP1.csv")

breaks<-c(0,7,15,60);limits<-c(-5,65)
pd <- position_dodge(width=1)

p_ArgTT_GLP1<-ggplot(ArgTT[ArgTT$Time%in%breaks,],
          aes(x=Time,y=GLP1.End*1000/3297.6,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "free_x")+
  geom_linerange(aes(ymin=GLP1.Start*1000/3297.6,ymax=GLP1.End*1000/3297.6,group=AnimalID),
                 linetype=3,size=0.3,alpha=0.8)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),size=0.3,linetype=2,alpha=0.8)+
  labs(x="Time (minutes)",y="Human or\nMouse GLP-1 (pM)")+
  geom_hline(yintercept = 2*40/6*1000/3297.6,linetype=3,alpha=0.7)+
  scale_x_continuous(breaks = breaks,limits=limits)+ #
  scale_y_continuous(breaks=pretty_breaks(4))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8),
        legend.title = element_blank())+
  geom_smooth(data = ff_ArgTT_GLP1, 
              aes(y = Estimate, ymin = Q2.5, 
                  ymax = Q97.5,
                  fill = FullGroup,colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1)+
  geom_point(data = ff_ArgTT_GLP1, 
             aes(y = Estimate,shape=FullGroup),
             position=pd,size=2,alpha=0.7)
p_ArgTT_GLP1

ggsave("GLP1_ArgTT_Bayes.png", path="./Figures/ArgTT",
       width = 22, height = 14, units = "cm")

##### Comparisons #####
GLP1tab1<-rbind(hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC)-
                           (Intercept)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+TimeNom7)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+TimeNom15)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+TimeNom60)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC)-
                           (Intercept)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC)-
                           (Intercept+FullGroupFemale_TC)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)>0")$hypothesis,
                
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)>0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_KC)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)=0")$hypothesis)

GLP1tab2<-rbind(         
  hypothesis(fit_ArgTT_GLP1,"(Intercept+TimeNom7)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+TimeNom15)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+TimeNom60)-
             (Intercept)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
             (Intercept+FullGroupFemale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
             (Intercept+FullGroupMale_FP)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
             (Intercept+FullGroupMale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
             (Intercept+FullGroupMale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis)

##Sex differences: absolute values
GLP1tab3<-rbind(hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC)-(Intercept+FullGroupFemale_KC)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)<0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)<0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)<0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupFemale_TC)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP)-(Intercept)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_GLP1,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis)

##### ArgTT: Gcg Stats #####
fitA <- brm(bf(GcgpM~FullGroup*TimeNom+(TimeNom|AnimalID)),
            iter=4000,
            warmup=1000,
            family="skew_normal",
            # prior = c(set_prior("normal(0,10)", class = "b"),
            #           set_prior("normal(30,10)", class = "Intercept")),
            control = list(adapt_delta=0.99),
            data = imps_ArgTT_Horm[[1]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitB <- brm(GcgpM~FullGroup*TimeNom+(TimeNom|AnimalID),
            iter=4000,
            warmup=1000,
            family="gamma",
            prior = c(set_prior("normal(-5,10)", class = "b"),
                      set_prior("normal(log(30),10)", class = "Intercept")),
            control = list(adapt_delta=0.8),
            data = imps_ArgTT_Horm[[1]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo")
fitB <- add_criterion(fitB, "loo")

loo_compare(fitA,fitB,criterion = "loo")

pp_check(fitA)
pp_check(fitB)

bayestestR::check_prior(fitA)

fit_ArgTT_Gcg <- brm_multiple(GcgpM~FullGroup*TimeNom+(TimeNom|AnimalID),
                              iter=8000,
                              warmup=1000,
                              thin=5,
                              family="skew_normal",
                              data = imps_ArgTT_Horm,
                              control = list(adapt_delta = 0.99,
                                             max_treedepth=12),
                              chains = 4,
                              silent=0)

load("./fit_ArgTT_Gcg.RData")

pp_check(fit_ArgTT_Gcg) # ~similar plots of observed and predicted values
Gcg_summary

max(fit_ArgTT_Gcg$rhats)

newdata4 = data.frame(Group = factor(rep(levels(ArgTT$Group),2*4)|>str_sort()),
                      Sex=factor(c(rep(levels(ArgTT$Sex),4*3))),
                      TimeNom = factor(rep(c(0,7,15,60),each=2)|>rep(3)))

newdata4$FullGroup<-paste(newdata4$Sex, newdata4$Group,sep="_")

ArgTT_Gcg_mod <- fitted(fit_ArgTT_Gcg,
               newdata = newdata4, 
               re_formula = NA) # extract the full MCMC

ff_ArgTT_Gcg <- ArgTT_Gcg_mod |>
  as_tibble() |>
  bind_cols(newdata4)

ff_ArgTT_Gcg$Time<-rep(c(0,7,15,60),each=2)|>rep(3)
ff_ArgTT_Gcg$FullGroup<-factor(ff_ArgTT_Gcg$FullGroup,levels=levels(weekly_monitoring$FullGroup))
ff_ArgTT_Gcg$Group<-factor(ff_ArgTT_Gcg$Group,levels=levels(weekly_monitoring$Group))
ff_ArgTT_Gcg$Sex<-factor(ff_ArgTT_Gcg$Sex,levels=levels(weekly_monitoring$Sex))

write_csv(ff_ArgTT_Gcg,file="./ArgTT_Gcg.csv")

breaks<-c(0,7,15,60);limits<-c(-5,65)

p_ArgTT_Gcg<-ggplot(ArgTT[ArgTT$Time%in%breaks,],
          aes(x=Time,y=Gcg.End*1000/3485,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "free_x")+
  geom_linerange(aes(ymin=Gcg.Start*1000/3485,ymax=Gcg.End*1000/3485,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8,show.legend = F)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),size=0.3,
            linetype=2,alpha=0.8,show.legend = F)+
  labs(x="Time (minutes)",y="Human or\nMouse Glucagon (pM)")+
  geom_hline(yintercept = 27*40/6*1000/3485,linetype=3,alpha=0.7)+
  scale_x_continuous(breaks = breaks,limits=limits)+ #
  scale_y_continuous(breaks=pretty_breaks(4))+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8),
        legend.title = element_blank())+
  geom_smooth(data = ff_ArgTT_Gcg, 
              aes(y = Estimate, 
                  ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup,colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1,show.legend = F)+
  geom_point(data = ff_ArgTT_Gcg, 
             aes(y = Estimate,shape=FullGroup),
             position=pd,size=2,alpha=0.7,show.legend = F)
p_ArgTT_Gcg

ggsave("Gcg_ArgTT_Bayes.png", path="./Figures/GSIS",
       width = 22, height = 14, units = "cm")

p<-ggplot(ArgTT[ArgTT$Time%in%breaks,],
          aes(x=Time,y=Gcg.End*1000/3485,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "free_x")+
  geom_linerange(aes(ymin=Gcg.Start*1000/3485,ymax=Gcg.End*1000/3485,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),position=pd,size=0.3,linetype=2,alpha=0.5)+
  labs(x="Time (minutes)",y="Human or Mouse Glucagon (pM)")+
  geom_hline(yintercept = 27*20/3*1000/3485,linetype=3,alpha=0.7)+
  scale_x_continuous(breaks = breaks,limits=limits)+ #
  scale_y_continuous(breaks=pretty_breaks(6))+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8),
        legend.title = element_blank())+
  geom_smooth(data = ff, 
              aes(y = Estimate, 
                  ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup,colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1/2)+
  geom_point(data = ff, 
             aes(y = Estimate,shape=FullGroup),
             position=pd,size=1.75,alpha=0.5)
p

ggsave("Gcg_ArgTT_Bayes.png", path="./Figures/ArgTT",
       width = 22, height = 14, units = "cm")

##### Comparisons #####
Gcgtab1<-rbind(hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC)-
                           (Intercept)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+TimeNom7)>0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+TimeNom15)>0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+TimeNom60)>0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC)-
                           (Intercept)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC)-
                           (Intercept+FullGroupFemale_TC)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)>0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)>0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)>0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_KC)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)=0")$hypothesis)

Gcgtab2<-rbind(         
  hypothesis(fit_ArgTT_Gcg,"(Intercept+TimeNom7)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+TimeNom15)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+TimeNom60)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
             (Intercept+FullGroupFemale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
             (Intercept+FullGroupMale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_TC)=0")$hypothesis)

##Sex differences: absolute values
Gcgtab3<-rbind(hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC)-(Intercept+FullGroupFemale_KC)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)<0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)<0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)<0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupFemale_TC)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP)-(Intercept)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
                hypothesis(fit_ArgTT_Gcg,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis)

##### ArgTT: Ins Stats #####
options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)


fitA <- brm(bf(InsnM~FullGroup*TimeNom+(TimeNom|AnimalID),
               sigma~FullGroup),
            iter=8000,
            warmup=1000,
            family="skew_normal",
            # prior = c(set_prior("normal(0,10)", class = "b"),
            #           set_prior("normal(30,10)", class = "Intercept")),
            control = list(adapt_delta=0.8),
            data = imps_ArgTT_Horm[[9]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitB <- brm(bf(InsnM~FullGroup*TimeNom+(1|AnimalID),
               shape~Group),
            iter=8000,
            warmup=1000,
            family="gamma",
            # prior = c(set_prior("normal(-5,10)", class = "b"),
            #           set_prior("normal(log(30),10)", class = "Intercept")),
            control = list(adapt_delta=0.8),
            data = imps_ArgTT_Horm[[9]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitC <- brm(InsnM~FullGroup*TimeNom+(1|AnimalID),
            iter=8000,
            warmup=1000,
            family="gamma",
            # prior = c(set_prior("normal(-5,10)", class = "b"),
            #           set_prior("normal(log(30),10)", class = "Intercept")),
            control = list(adapt_delta=0.8),
            data = imps_ArgTT_Horm[[9]],
            # chains=1,
            silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

fitA <- add_criterion(fitA, "loo")
fitB <- add_criterion(fitB, "loo")
fitC <- add_criterion(fitC, "loo")

loo_compare(fitA,fitB,fitC,criterion = "loo")

pp_check(fitA)
pp_check(fitB)
pp_check(fitC)

bayestestR::check_prior(fitA)

fit_ArgTT_Ins <- brm_multiple(InsnM~FullGroup*TimeNom+(1|AnimalID),
                              iter=8000,
                              warmup=1000,
                              thin=5,
                              data = imps_ArgTT_Horm,
                              family="gamma",
                              #prior = set_prior("normal(0,1)", class = "b"),
                              control = list(adapt_delta = 0.9,
                                             max_treedepth=12),
                              #chains = 4,
                              silent=0)

load("./fit_ArgTT_Ins.RData")

pp_check(fit_ArgTT_Ins) # ~similar plots of observed and predicted values
summary(fit_ArgTT_Ins)
max(fit_ArgTT_Ins$rhats) #no issues

newdata4 = data.frame(Group = factor(rep(levels(ArgTT$Group),2*4)|>str_sort()),
                      Sex=factor(c(rep(levels(ArgTT$Sex),4*3))),
                      TimeNom = factor(rep(c(0,7,15,60),each=2)|>rep(3)))

newdata4$FullGroup<-paste(newdata4$Sex, newdata4$Group,sep="_")

ArgTT_Ins_mod <- fitted(fit_ArgTT_Ins,
               newdata = newdata4, 
               re_formula = NA) # extract the full MCMC

ff_ArgTT_Ins <- ArgTT_Ins_mod |>
  as_tibble() |>
  bind_cols(newdata4)

ff_ArgTT_Ins$Time<-rep(c(0,7,15,60),each=2)|>rep(3)
ff_ArgTT_Ins$FullGroup<-factor(ff_ArgTT_Ins$FullGroup,levels=levels(weekly_monitoring$FullGroup))
ff_ArgTT_Ins$Group<-factor(ff_ArgTT_Ins$Group,levels=levels(weekly_monitoring$Group))
ff_ArgTT_Ins$Sex<-factor(ff_ArgTT_Ins$Sex,levels=levels(weekly_monitoring$Sex))

write_csv(ff_ArgTT_Ins,file="./ArgTT_Ins.csv")

breaks<-c(0,7,15,60);limits<-c(-5,65)

p_ArgTT_Ins<-ggplot(ArgTT[ArgTT$Time%in%breaks,],
          aes(x=Time,y=Ins.End*1000/5808,colour=FullGroup,group=FullGroup,fill=FullGroup))+
  facet_grid(Sex~.,scales = "free_x")+
  geom_linerange(aes(ymin=Ins.Start*1000/5808,ymax=Ins.End*1000/5808,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8,show.legend = F)+
  scale_colour_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                     "Female: Kidney Capsule",
                                                                     "Female: TheraCyte",
                                                                     "Male: Fat Pad",
                                                                     "Male: Kidney Capsule",
                                                                     "Male: TheraCyte"))+
  scale_shape_manual(name="FullGroup", values=GroupShapes,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  scale_fill_manual(name="FullGroup", values=GroupPalette,labels=c("Female: Fat Pad",
                                                                   "Female: Kidney Capsule",
                                                                   "Female: TheraCyte",
                                                                   "Male: Fat Pad",
                                                                   "Male: Kidney Capsule",
                                                                   "Male: TheraCyte"))+
  geom_line(aes(colour=FullGroup,group=AnimalID),
            size=0.3,linetype=2,alpha=0.8,show.legend = F)+
  labs(x="Time (minutes)",y="\nHuman Insulin (nM)")+
  geom_hline(yintercept = 8*40/6/1000/5808,
             linetype=3,alpha=0.7)+
  scale_x_continuous(limits=limits,breaks = breaks)+
  facet_grid(Sex~.,scales = "free_x")+
  # guides(fill=guide_legend(ncol=2),
  #        colour=guide_legend(ncol=2),
  #        shape=guide_legend(ncol=2))+
  theme(axis.title = element_text(family = "Arial", color="black", size=8))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=8))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=8))+
  theme(legend.text = element_text(family = "Arial",color="black",size=8),
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=8))+
  # theme(strip.text = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_smooth(data = ff_ArgTT_Ins, 
              aes(y = Estimate, 
                  ymin = Q2.5, ymax = Q97.5,
                  fill = FullGroup,colour=FullGroup),
              stat = "identity", 
              alpha = 1/4, size = 1,show.legend = F)+
  geom_point(data = ff_ArgTT_Ins, 
             aes(y = Estimate,shape=FullGroup),
             position=pd,size=2,alpha=0.8,show.legend = F)
p_ArgTT_Ins

ggsave("Ins_ArgTT_Bayes.png", path="./Figures/ArgTT",width = 22, height = 14, units = "cm")

##### Comparisons #####
Instab1<-rbind(hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC)-
                           (Intercept)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+TimeNom7)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC)-
                           (Intercept)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
                           (Intercept+TimeNom7)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
                           (Intercept+TimeNom15)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
                           (Intercept+TimeNom60)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC)-
                           (Intercept+FullGroupFemale_TC)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)<0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC)-
                           (Intercept+FullGroupMale_FP)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_FP)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC)-
                           (Intercept+FullGroupMale_KC)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)=0")$hypothesis)

Instab2<-rbind(         
  hypothesis(fit_ArgTT_Ins,"(Intercept+TimeNom7)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+TimeNom15)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+TimeNom60)-
             (Intercept)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)-
             (Intercept+FullGroupFemale_KC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)-
             (Intercept+FullGroupFemale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)-
             (Intercept+FullGroupFemale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)-
             (Intercept+FullGroupFemale_TC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
             (Intercept+FullGroupMale_FP)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
             (Intercept+FullGroupMale_KC)=0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
             (Intercept+FullGroupMale_TC)>0")$hypothesis,
  hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
             (Intercept+FullGroupMale_TC)=0")$hypothesis)

##Sex differences: absolute values
Instab3<-rbind(hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC)-(Intercept+FullGroupFemale_KC)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom7+FullGroupMale_KC:TimeNom7)-
                           (Intercept+FullGroupFemale_KC+TimeNom7+FullGroupFemale_KC:TimeNom7)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom15+FullGroupMale_KC:TimeNom15)-
                           (Intercept+FullGroupFemale_KC+TimeNom15+FullGroupFemale_KC:TimeNom15)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_KC+TimeNom60+FullGroupMale_KC:TimeNom60)-
                           (Intercept+FullGroupFemale_KC+TimeNom60+FullGroupFemale_KC:TimeNom60)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC)-(Intercept+FullGroupFemale_TC)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom7+FullGroupMale_TC:TimeNom7)-
                           (Intercept+FullGroupFemale_TC+TimeNom7+FullGroupFemale_TC:TimeNom7)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom15+FullGroupMale_TC:TimeNom15)-
                           (Intercept+FullGroupFemale_TC+TimeNom15+FullGroupFemale_TC:TimeNom15)>0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_TC+TimeNom60+FullGroupMale_TC:TimeNom60)-
                           (Intercept+FullGroupFemale_TC+TimeNom60+FullGroupFemale_TC:TimeNom60)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP)-(Intercept)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP+TimeNom7+FullGroupMale_FP:TimeNom7)-
                           (Intercept+TimeNom7)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP+TimeNom15+FullGroupMale_FP:TimeNom15)-
                           (Intercept+TimeNom15)=0")$hypothesis,
               hypothesis(fit_ArgTT_Ins,"(Intercept+FullGroupMale_FP+TimeNom60+FullGroupMale_FP:TimeNom60)-
                           (Intercept+TimeNom60)=0")$hypothesis)

##### Results Tables #####
######Supplementary Table 2 - BW#####
BWtab1<-filter(BWtab1,Star=="*")
BWtab1$Sex<-ifelse(grepl("Female",BWtab1$Hypothesis),"Female","Male") 
BWtab1$ContrSite<-str_split_fixed(BWtab1$Hypothesis, boundary("word"),6)[,1]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
BWtab1$RefSite<-str_split_fixed(BWtab1$Hypothesis, boundary("word"),6)[,3]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
BWtab1$RefSite[BWtab1$RefSite=="0"]<-"FP"

BWtab1$FullGroupRef=paste(BWtab1$Sex,BWtab1$RefSite,sep="_")
BWtab1$FullGroupContr=paste(BWtab1$Sex,BWtab1$ContrSite,sep="_")

BWtab1$Weeks<-str_split_fixed(BWtab1$Hypothesis,"[:)-]",3)[,2]|>
  str_replace_all("DaysNom","")|>
  as.numeric()/7

BWtab1$Weeks<-round(BWtab1$Weeks)

SuppTab2<-BWtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                             Evid.Ratio,Sex,Weeks,RefSite,ContrSite,
                             FullGroupRef,FullGroupContr))|>
  mutate(Evid.Ratio=case_when(BWtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              BWtab1$Evid.Ratio<=10^(2)~"Very strong",
                              BWtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Weeks,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 2. Weekly Body Weight Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (g)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,RefSite,ContrSite,Weeks)
  )|>
  cols_hide(columns=c(FullGroupRef,FullGroupContr))|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Body Weight (g)"),
             Weeks = "Weeks Post-Implant")|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab2

gtsave(SuppTab2,"SuppTab2.pdf",path="./Tables")

######Supplementary Table 3&4 - BG#####
BGtab1<-filter(BGtab1,Star=="*")
BGtab1$Sex<-ifelse(grepl("Female",BGtab1$Hypothesis),"Female","Male") #wrong - need female FP group
BGtab1$Site<-str_split_fixed(BGtab1$Hypothesis,"[:+]",3)[,2]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
BGtab1$Site[!grepl("FullGroup",BGtab1$Hypothesis)]<-"FP"
BGtab1$Sex[!grepl("FullGroup",BGtab1$Hypothesis)]<-"Female"

BGtab1$Weeks<-str_extract(BGtab1$Hypothesis,"DaysNom[0-9]{1,3}")|>
  str_replace_all("DaysNom","")|>
  as.numeric()/7
BGtab1$Weeks<-round(BGtab1$Weeks)

SuppTab3<-BGtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                             Evid.Ratio,Sex,Weeks,Site))|>
  mutate(Evid.Ratio=case_when(BGtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              BGtab1$Evid.Ratio<=10^(2)~"Very strong",
                              BGtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Weeks)|>
  gt() |>
  row_group_order(c("Female - FP","Female - KC","Female - TC",
                    "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 3. Weekly Blood Glucose Compared to Pre-implant")) |>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Weeks)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Weeks = html("Weeks<br>Post-Implant"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab3

gtsave(SuppTab3,"SuppTab3.pdf",path="./Tables")


BGtab2<-filter(BGtab2,Star=="*")
BGtab2$Sex<-ifelse(grepl("Female",BGtab2$Hypothesis),"Female","Male") 
BGtab2$ContrSite<-str_split_fixed(BGtab2$Hypothesis, boundary("word"),6)[,1]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
BGtab2$RefSite<-str_split_fixed(BGtab2$Hypothesis, boundary("word"),6)[,3]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
BGtab2$RefSite[BGtab2$RefSite=="0"]<-"FP"

BGtab2$FullGroupRef=paste(BGtab2$Sex,BGtab2$RefSite,sep="_")
BGtab2$FullGroupContr=paste(BGtab2$Sex,BGtab2$ContrSite,sep="_")


BGtab2$Weeks<-str_split_fixed(BGtab2$Hypothesis,"[:)-]",3)[,2]|>
  str_replace_all("DaysNom","")|>
  as.numeric()/7

BGtab2$Weeks<-round(BGtab2$Weeks)

SuppTab4<-BGtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                             Evid.Ratio,Sex,Weeks,RefSite,ContrSite,
                             FullGroupRef,FullGroupContr))|>
  mutate(Evid.Ratio=case_when(BGtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              BGtab2$Evid.Ratio<=10^(2)~"Very strong",
                              BGtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Weeks,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 4. Weekly Blood Glucose Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Weeks)
  )|>
  cols_hide(columns=c(FullGroupRef,FullGroupContr))|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Weeks = html("Weeks<br>Post-Implant"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab4

gtsave(SuppTab4,"SuppTab4.pdf",path="./Tables")

######Supplementary Table 5: GSIS BG#####
GSISBGtab1<-filter(GSISBGtab1,Star=="*")
h<-str_replace(GSISBGtab1$Hypothesis,"[(]exp[(]","")|>str_replace("[)]{2}","")
GSISBGtab1$Sex<-ifelse(grepl("Female",h),"Female","Male") 
GSISBGtab1$ContrSite<-str_split_fixed(h, boundary("word"),18)[,2]
GSISBGtab1$ContrSite<-str_replace_all(GSISBGtab1$ContrSite,"FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GSISBGtab1$RefSite<-str_split_fixed(h, boundary("word"),18)[,11]
GSISBGtab1$RefSite<-str_replace_all(GSISBGtab1$RefSite,"FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GSISBGtab1$RefSite<-case_when(grepl("Weeks",GSISBGtab1$RefSite)~"FP",
                              TRUE~GSISBGtab1$RefSite)

GSISBGtab1$FullGroupRef=paste(GSISBGtab1$Sex,GSISBGtab1$RefSite,sep="_")
GSISBGtab1$FullGroupContr=paste(GSISBGtab1$Sex,GSISBGtab1$ContrSite,sep="_")

GSISBGtab1$Weeks<-str_split_fixed(h, boundary("word"),18)[,3]
GSISBGtab1$Weeks<-str_replace_all(GSISBGtab1$Weeks,"WeeksNom","")|>paste(" Weeks")

GSISBGtab1$Time<-str_split_fixed(h, boundary("word"),18)[,5]|>
  str_replace_all("TimeNom","")|>as.numeric()

SuppTab5<-GSISBGtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Weeks,Time,RefSite,ContrSite,
                                 FullGroupRef,FullGroupContr))|>
  mutate(Evid.Ratio=case_when(GSISBGtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              GSISBGtab1$Evid.Ratio<=10^(2)~"Very strong",
                              GSISBGtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Weeks) |> 
  arrange(Weeks,Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 5. Glucose-Stimulated Blood Glucose Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Weeks,Time)
  )|>
  cols_hide(columns=c(FullGroupRef,FullGroupContr))|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Time = html("Time<br>Post-Glucose (min)"),
             Weeks = html("Weeks<br>Post-Implant"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab5

gtsave(SuppTab5,"SuppTab5.pdf",path="./Tables")


GSISBGtab2<-filter(GSISBGtab2,Star=="*")
GSISBGtab2$Sex<-ifelse(grepl("Female",GSISBGtab2$Hypothesis),"Female","Male")
GSISBGtab2$Sex[!grepl("FullGroup",GSISBGtab2$Hypothesis)]<-"Female"

GSISBGtab2$Site<-str_extract(GSISBGtab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GSISBGtab2$Site[is.na(GSISBGtab2$Site)]<-"FP"

GSISBGtab2$Time<-str_extract(GSISBGtab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
GSISBGtab2$Time[is.na(GSISBGtab2$Time)]<-0

GSISBGtab2$Weeks<-str_extract(GSISBGtab2$Hypothesis,"WeeksNom[0-9]{1,2}")|>
  str_replace_all("sNom"," ")

GSISBGtab2$Weeks<-factor(GSISBGtab2$Weeks,levels=c("Week 4","Week 8","Week 10","Week 12","Week 16","Week 21"))

SuppTab5B<-GSISBGtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Weeks,Site,Time))|>
  mutate(Evid.Ratio=case_when(GSISBGtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              GSISBGtab2$Evid.Ratio<=10^(2)~"Very strong",
                              GSISBGtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Weeks,Time)|>
  gt() |>
  row_group_order(c("Female - FP","Female - KC","Female - TC",
                    "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 5B. Glucose-Stimulated Blood Glucose Compared to Fasted")) |>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Weeks,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Weeks = html("Weeks<br>Post-Implant"),
             Time = html("Time<br>Post-Glucose (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab5B

gtsave(SuppTab5B,"SuppTab5B.pdf",path="./Tables")


######Supplementary Table 6-8: GSIS CP#####
GSISCPtab1<-filter(GSISCPtab1,Star=="*")
GSISCPtab1$Sex<-ifelse(grepl("Female",GSISCPtab1$Hypothesis),"Female","Male")

h<-str_replace_all(GSISCPtab1$Hypothesis,"[(]exp[(]","")|>str_replace("[)]{2}","")|>str_replace("Intercept","")

GSISCPtab1$ContrSite<-str_split_fixed(h, "[+]",3)[,2]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GSISCPtab1$RefSite<-str_split_fixed(h, "Intercept[+]",2)[,2]
GSISCPtab1$RefSite<-str_split_fixed(GSISCPtab1$RefSite,"[+]",2)[,1]|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GSISCPtab1$RefSite<-case_when(grepl("Weeks",GSISCPtab1$RefSite)~"FP",
                              TRUE~GSISCPtab1$RefSite)

GSISCPtab1$Weeks<-str_extract(h,"WeeksNom[0-9]{1,2}")|>
  str_replace_all("sNom"," ")

GSISCPtab1$Time<-str_extract(h,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
GSISCPtab1$Time[is.na(GSISCPtab1$Time)]<-0

SuppTab6<-GSISCPtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Weeks,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(GSISCPtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              GSISCPtab1$Evid.Ratio<=10^(2)~"Very strong",
                              GSISCPtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Weeks) |> 
  arrange(Weeks,Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 6. Glucose-Stimulated C-peptide Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Weeks,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>C-peptide (nM)"),
             Time = html("Time<br>Post-Glucose (min)"),
             Weeks = html("Weeks<br>Post-Implant"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab6

gtsave(SuppTab6,"SuppTab6.pdf",path="./Tables")

GSISCPtab2<-filter(GSISCPtab2,Star=="*")
GSISCPtab2$Sex<-ifelse(grepl("Female",GSISCPtab2$Hypothesis),"Female","Male")
GSISCPtab2$Sex[!grepl("FullGroup",GSISCPtab2$Hypothesis)]<-"Female"

GSISCPtab2$Site<-str_extract(GSISCPtab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
    str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GSISCPtab2$Site[is.na(GSISCPtab2$Site)]<-"FP"

GSISCPtab2$Time<-str_extract(GSISCPtab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
GSISCPtab2$Time[is.na(GSISCPtab2$Time)]<-0

GSISCPtab2$Weeks<-str_extract(GSISCPtab2$Hypothesis,"WeeksNom[0-9]{1,2}")|>
  str_replace_all("sNom"," ")

GSISCPtab2$Weeks<-factor(GSISCPtab2$Weeks,levels=c("Week 8","Week 10","Week 12","Week 16","Week 21"))

SuppTab7<-GSISCPtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Weeks,Site,Time))|>
  mutate(Evid.Ratio=case_when(GSISCPtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              GSISCPtab2$Evid.Ratio<=10^(2)~"Very strong",
                              GSISCPtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Weeks,Time)|>
  gt() |>
  row_group_order(c("Female - FP","Female - KC","Female - TC",
                    "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 7. Glucose-Stimulated C-peptide Compared to Fasted")) |>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Weeks,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>C-peptide (nM)"),
             Weeks = html("Weeks<br>Post-Implant"),
             Time = html("Time<br>Post-Glucose (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab7

gtsave(SuppTab7,"SuppTab7.pdf",path="./Tables")

GSISCPtab3<-filter(GSISCPtab3,Star=="*")
GSISCPtab3$Site<-str_extract(GSISCPtab3$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")

GSISCPtab3$ContrSex<-"Female"
GSISCPtab3$RefSex<-"Male"

GSISCPtab3$Time<-str_extract(GSISCPtab3$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
GSISCPtab3$Time[is.na(GSISCPtab3$Time)]<-0

GSISCPtab3$Weeks<-str_extract(GSISCPtab3$Hypothesis,"WeeksNom[0-9]{1,2}")|>
  str_replace_all("sNom"," ")

SuppTab8<-GSISCPtab3 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Weeks,Time,Site))|>
  mutate(Evid.Ratio=case_when(GSISCPtab3$Evid.Ratio<=10^(3/2)~"Strong",
                              GSISCPtab3$Evid.Ratio<=10^(2)~"Very strong",
                              GSISCPtab3$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Site,Weeks) |> 
  arrange(Site,Weeks,Time)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 8. Glucose-Stimulated C-peptide Compared Within Group; Female Compared to Male")) |>
  # tab_spanner(
  #   label = "Sex Comparisons",
  #   columns = c(ContrSex,RefSex)
  # )|>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Site,Weeks,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             Estimate = html("Difference in<br>C-peptide (nM)"),
             Time = html("Time<br>Post-Glucose (min)"),
             Weeks = html("Weeks<br>Post-Implant"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab8

gtsave(SuppTab8,"SuppTab8.pdf",path="./Tables")

######Supplementary Table 9-11: ITT BG#####
ITTBGtab1<-filter(ITTBGtab1,Star=="*")
ITTBGtab1$Sex<-ifelse(grepl("Female",ITTBGtab1$Hypothesis),"Female","Male") 
ITTBGtab1$ContrSite<-str_extract(ITTBGtab1$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ITTBGtab1$ContrSite[is.na(ITTBGtab1$ContrSite)]<-"FP"

ITTBGtab1$RefSite<-str_split_fixed(ITTBGtab1$Hypothesis, "-\\(Intercept",2)[,2]|>
  str_extract("FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ITTBGtab1$RefSite[is.na(ITTBGtab1$RefSite)]<-"FP"

ITTBGtab1$Time<-str_extract(ITTBGtab1$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
ITTBGtab1$Time[is.na(ITTBGtab1$Time)]<-0

SuppTab9<-ITTBGtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                                Evid.Ratio,Sex,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(ITTBGtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              ITTBGtab1$Evid.Ratio<=10^(2)~"Very strong",
                              ITTBGtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(ContrSite,RefSite,Time)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 9. Blood Glucose Changes During ITT Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"), #Δ 
             Time = html("Time<br>Post-Insulin (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ITT=Insulin Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  #<br>Results reported at time 0 are absolute differences in fasted blood glucose between groups;<br>
  #other results are the differences between groups in the change in blood glucose from the indicated<br>time compared to time 0.
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )) 

SuppTab9

gtsave(SuppTab9,"SuppTab9.pdf",path="./Tables")

ITTBGtab2<-filter(ITTBGtab2,Star=="*")
ITTBGtab2$Sex<-str_extract(ITTBGtab2$Hypothesis,"FullGroup(Female|Male)")|>
  str_replace_all("FullGroup","")
ITTBGtab2$Sex[is.na(ITTBGtab2$Sex)]<-"Female"

ITTBGtab2$Site<-str_extract(ITTBGtab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ITTBGtab2$Site[is.na(ITTBGtab2$Site)]<-"FP"

ITTBGtab2$Time<-str_extract(ITTBGtab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()

SuppTab10<-ITTBGtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Site,Time))|>
  mutate(Evid.Ratio=case_when(ITTBGtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              ITTBGtab2$Evid.Ratio<=10^(2)~"Very strong",
                              ITTBGtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Sex,Time)|>
  gt() |>
  row_group_order(c("Female - KC","Female - TC",
                    "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 10. Blood Glucose During ITT Compared to Fasted")) |>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Time = html("Time<br>Post-Insulin (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("ITT=Insulin Tolerance Test.<br>If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab10

gtsave(SuppTab10,"SuppTab10.pdf",path="./Tables")

# ITTBGtab3<-filter(ITTBGtab3,Star=="*")
# ITTBGtab3$ContrSex<-str_extract(ITTBGtab3$Hypothesis,"FullGroup(Female|Male)")|>
#   str_replace_all("FullGroup","")
# ITTBGtab3$ContrSex[is.na(ITTBGtab3$ContrSex)]<-"Female"
# 
# ITTBGtab3$RefSex<-str_split_fixed(ITTBGtab3$Hypothesis, "-\\(Intercept",2)[,2]|>
#   str_extract("FullGroup(Female|Male)")|>
#   str_replace_all("FullGroup","")
# ITTBGtab3$RefSex[is.na(ITTBGtab3$RefSex)]<-"Male"
# 
# ITTBGtab3$Time<-str_extract(ITTBGtab3$Hypothesis,"TimeNom[0-9]{1,3}")|>
#   str_replace_all("TimeNom","")|>as.numeric()
# ITTBGtab3$Time[is.na(ITTBGtab3$Time)]<-0
# 
# ITTBGtab3<-filter(ITTBGtab3,Star=="*")
# ITTBGtab3$Site<-str_split_fixed(ITTBGtab3$Hypothesis, boundary("word"),5)[,1]|>
#   str_split_fixed("[:]",2)
# ITTBGtab3$Site<-ITTBGtab3$Site[,1]|>
#   str_replace_all("FullGroupFemale_","")|>
#   str_replace_all("FullGroupMale_","")
# 
# ITTBGtab3$ContrSex<-str_split_fixed(ITTBGtab3$Hypothesis, boundary("word"),6)[,1]|>
#   str_split_fixed("[_]",2)
# ITTBGtab3$ContrSex<-ITTBGtab3$ContrSex[,1]|>str_replace_all("FullGroup","")
# 
# ITTBGtab3$RefSex<-str_split_fixed(ITTBGtab3$Hypothesis, boundary("word"),5)[,3]|>
#   str_split_fixed("[_]",2)
# ITTBGtab3$RefSex<-ITTBGtab3$RefSex[,1]|>str_replace_all("FullGroup","")
# ITTBGtab3$RefSex[ITTBGtab3$RefSex=="0"]<-"Female"
# ITTBGtab3$RefSex[ITTBGtab3$RefSex==""]<-"Female"
# 
# ITTBGtab3$FullGroupRef=paste(ITTBGtab3$RefSex,ITTBGtab3$Site,sep="_")
# ITTBGtab3$FullGroupContr=paste(ITTBGtab3$ContrSex,ITTBGtab3$Site,sep="_")
# 
# ITTBGtab3$Time<-str_split_fixed(ITTBGtab3$Hypothesis, boundary("word"),5)[,1]|>
#   str_split_fixed("[:]",2)
# ITTBGtab3$Time<-ITTBGtab3$Time[,2]|>str_replace_all("TimeNom","")
# ITTBGtab3$Time[ITTBGtab3$Time==""]<-"0"
# ITTBGtab3$Time<-as.numeric(ITTBGtab3$Time)
# 
# SuppTab11<-ITTBGtab3 |> select(c(Estimate,CI.Lower,CI.Upper,
#                                  Evid.Ratio,RefSex,ContrSex,Time,Site,
#                                  FullGroupRef,FullGroupContr))|>
#   mutate(Evid.Ratio=case_when(ITTBGtab3$Evid.Ratio<=10^(3/2)~"Strong",
#                               ITTBGtab3$Evid.Ratio<=10^(2)~"Very strong",
#                               ITTBGtab3$Evid.Ratio>10^(2)~"Decisive")) |>
#   group_by(Site) |> 
#   arrange(Time,ContrSex,RefSex)|>
#   gt() |>
#   # row_group_order(c("Female - FP","Female - KC","Female - TC",
#   #                   "Male - FP","Male - KC", "Male - TC"))|>
#   tab_header(title = md("Supplementary Table 11. Blood Glucose Changes During ITT Compared Within Group, Between Sex")) |>
#   tab_spanner(
#     label = "Sex Comparisons",
#     columns = c(ContrSex,RefSex)
#   )|>
#   tab_spanner(
#     label = "95% Credibility Interval (mM)",
#     columns = c(CI.Lower,CI.Upper)
#   )|>
#   cols_move_to_start(
#     columns = c(Site,ContrSex,RefSex,Time)
#   )|>
#   cols_hide(columns=c(FullGroupRef,FullGroupContr))|>
#   cols_label(CI.Lower = "Lower Bound",
#              CI.Upper = "Upper Bound",
#              Evid.Ratio = "Strength of Evidence",
#              RefSex = "Reference",
#              ContrSex = "Contrast",
#              Estimate = html("Difference in<br>Δ Blood Glucose (mM)"),
#              Time = html("Time<br>Post-Insulin (min)")
#   )|>
#   fmt_number(
#     columns = c(Estimate,CI.Lower,CI.Upper),
#     decimals = 1,
#     use_seps = TRUE
#   ) |>
#   # data_color( 
#   #   columns = c(RefSite), 
#   #   colors = scales::col_factor( 
#   #     palette = GroupPalette,
#   #     domain = levels(weekly_monitoring$FullGroup)
#   #   )) |>
#   # data_color(
#   #   columns = c(ContrSite),
#   #   colors = scales::col_factor(
#   #     palette = GroupPalette,
#   #     domain =  levels(weekly_monitoring$FullGroup)
# #   )) |>
# tab_source_note(md("ITT=Insulin Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.<br>Results reported are the differences between sexes in the change in blood glucose from the indicated<br>time compared to time 0.")) |>
#   tab_style(
#     locations = cells_column_labels(columns = everything()),
#     style     = list(
#       cell_borders(sides = "bottom", weight = px(3)),
#       cell_text(weight = "bold")
#     )
#   ) |>
#   tab_style(
#     locations = cells_title(groups = "title"),
#     style     = list(
#       cell_text(weight = "bold", size = 24)
#     )
#   ) |>
#   tab_options(
#     table.font.style = "Arial",
#     table.font.color = "black"
#   )
# SuppTab11
# 
# gtsave(SuppTab11,"SuppTab11.pdf",path="./Tables")

#######Supplementary Table 11-13: ITT CP######
ITTCPtab1<-filter(ITTCPtab1,Star=="*")
ITTCPtab1$Sex<-ifelse(grepl("Female",ITTCPtab1$Hypothesis),"Female","Male") 
ITTCPtab1$ContrSite<-str_extract(ITTCPtab1$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ITTCPtab1$ContrSite[is.na(ITTCPtab1$ContrSite)]<-"FP"

ITTCPtab1$RefSite<-str_split_fixed(ITTCPtab1$Hypothesis, "-\\(Intercept",2)[,2]|>
  str_extract("FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ITTCPtab1$RefSite[is.na(ITTCPtab1$RefSite)]<-"FP"

ITTCPtab1$Time<-str_extract(ITTCPtab1$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
ITTCPtab1$Time[is.na(ITTCPtab1$Time)]<-0

SuppTab11<-ITTCPtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(ITTCPtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              ITTCPtab1$Evid.Ratio<=10^(2)~"Very strong",
                              ITTCPtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 11. C-peptide During ITT Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>C-peptide (nM)"),
             Time = html("Time<br>Post-Insulin (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ITT=Insulin Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.<br>Results are absolute differences in C-peptide between groups.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )) 

SuppTab11
#is there a text wrapping option??

gtsave(SuppTab11,"SuppTab11.pdf",path="./Tables")

ITTCPtab2<-filter(ITTCPtab2,Star=="*")
ITTCPtab2$Sex<-str_extract(ITTCPtab2$Hypothesis,"FullGroup(Female|Male)")|>
  str_replace_all("FullGroup","")
ITTCPtab2$Sex[is.na(ITTCPtab2$Sex)]<-"Female"

ITTCPtab2$Site<-str_extract(ITTCPtab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ITTCPtab2$Site[is.na(ITTCPtab2$Site)]<-"FP"

ITTCPtab2$Time<-str_extract(ITTCPtab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()


SuppTab12<-ITTCPtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,Sex,Site,Time))|>
  mutate(Evid.Ratio=case_when(ITTCPtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              ITTCPtab2$Evid.Ratio<=10^(2)~"Very strong",
                              ITTCPtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Time)|>
  gt() |>
  row_group_order(c("Female - FP", "Female - TC", "Female - KC", "Male - FP", 
                    "Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 12. C-peptide During ITT Compared to Fasted")) |>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>C-peptide (nM)"),
             Time = html("Time<br>Post-Insulin (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("ITT=Insulin Tolerance Test.<br>If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab12

gtsave(SuppTab12,"SuppTab12.pdf",path="./Tables")

ITTCPtab3<-filter(ITTCPtab3,Star=="*")
ITTCPtab3$Site<-str_extract(ITTCPtab3$Hypothesis, "FullGroup[Fem|M]ale[_](FP|TC|KC)")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")

ITTCPtab3$ContrSex<-str_extract_all(ITTCPtab3$Hypothesis, "FullGroup[Fem|M]ale[_](FP|TC|KC)",simplify = T)[,1]|>
  str_replace_all("FullGroup","")|>
  str_replace_all("[_](FP|TC|KC)","")

ITTCPtab3$RefSex<-case_when(ITTCPtab3$ContrSex=="Male"~"Female",
                            TRUE~ITTCPtab3$ContrSex)

ITTCPtab3$Time<-str_extract(ITTCPtab3$Hypothesis, "TimeNom[0-9]{2}")|>
  str_replace_all("TimeNom","")
ITTCPtab3$Time[is.na(ITTCPtab3$Time)]<-"0"|>as.numeric()

SuppTab13<-ITTCPtab3 |> select(c(Estimate,CI.Lower,CI.Upper,
                                 Evid.Ratio,RefSex,ContrSex,Time,Site))|>
  mutate(Evid.Ratio=case_when(ITTCPtab3$Evid.Ratio<=10^(3/2)~"Strong",
                              ITTCPtab3$Evid.Ratio<=10^(2)~"Very strong",
                              ITTCPtab3$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Site) |> 
  arrange(Time,ContrSex,RefSex)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 13. C-peptide During ITT Compared Within Group, Between Sex")) |>
  tab_spanner(
    label = "Sex Comparisons",
    columns = c(ContrSex,RefSex)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  cols_move_to_start(
    columns = c(Site,ContrSex,RefSex,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSex = "Reference",
             ContrSex = "Contrast",
             Estimate = html("Difference in<br>C-peptide (nM)"),
             Time = html("Time<br>Post-Insulin (min)")
  )|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ITT=Insulin Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab13

gtsave(SuppTab13,"SuppTab13.pdf",path="./Tables")

#######Supplementary Table 14-15: ArgTT BG######
ArgTTBGtab1<-filter(ArgTTBGtab1,Star=="*")
ArgTTBGtab1$Sex<-ifelse(grepl("Female",ArgTTBGtab1$Hypothesis),"Female","Male") 
ArgTTBGtab1$ContrSite<-str_extract(ArgTTBGtab1$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ArgTTBGtab1$ContrSite[is.na(ArgTTBGtab1$ContrSite)]<-"FP"

ArgTTBGtab1$RefSite<-str_split_fixed(ArgTTBGtab1$Hypothesis, "-\\(Intercept",2)[,2]|>
  str_extract("FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ArgTTBGtab1$RefSite[is.na(ArgTTBGtab1$RefSite)]<-"FP"

ArgTTBGtab1$Time<-str_extract(ArgTTBGtab1$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
ArgTTBGtab1$Time[is.na(ArgTTBGtab1$Time)]<-0

SuppTab14<-ArgTTBGtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                                   Evid.Ratio,Sex,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(ArgTTBGtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              ArgTTBGtab1$Evid.Ratio<=10^(2)~"Very strong",
                              ArgTTBGtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 14. Blood Glucose Changes During ArgTT Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.<br>")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )) 

SuppTab14
#is there a text wrapping option??

gtsave(SuppTab14,"SuppTab14.pdf",path="./Tables")


ArgTTBGtab2<-filter(ArgTTBGtab2,Star=="*")
ArgTTBGtab2$Sex<-str_extract(ArgTTBGtab2$Hypothesis,"FullGroup(Female|Male)")|>
  str_replace_all("FullGroup","")
ArgTTBGtab2$Sex[is.na(ArgTTBGtab2$Sex)]<-"Female"

ArgTTBGtab2$Site<-str_extract(ArgTTBGtab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
ArgTTBGtab2$Site[is.na(ArgTTBGtab2$Site)]<-"FP"

ArgTTBGtab2$Time<-str_extract(ArgTTBGtab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()

SuppTab15<-ArgTTBGtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                                   Evid.Ratio,Sex,Site,Time))|>
  mutate(Evid.Ratio=case_when(ArgTTBGtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              ArgTTBGtab2$Evid.Ratio<=10^(2)~"Very strong",
                              ArgTTBGtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Time,Site,Sex)|>
  gt() |>
  row_group_order(c("Female - FP","Female - KC","Female - TC",
                    "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 15. Blood Glucose During ArgTT Compared to Fasted")) |>
  tab_spanner(
    label = "95% Credibility Interval (mM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>Blood Glucose (mM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("ArgTT=Arginine Tolerance Test.<br>If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab15

gtsave(SuppTab15,"SuppTab15.pdf",path="./Tables")

# ArgTTBGtab3<-filter(ArgTTBGtab3,Star=="*")
# ArgTTBGtab3$Site<-str_split_fixed(ArgTTBGtab3$Hypothesis, boundary("word"),5)[,1]|>
#   str_split_fixed("[:]",2)
# ArgTTBGtab3$Site<-ArgTTBGtab3$Site[,1]|>
#   str_replace_all("FullGroupFemale_","")|>
#   str_replace_all("FullGroupMale_","")
# 
# ArgTTBGtab3$ContrSex<-str_split_fixed(ArgTTBGtab3$Hypothesis, boundary("word"),6)[,1]|>
#   str_split_fixed("[_]",2)
# ArgTTBGtab3$ContrSex<-ArgTTBGtab3$ContrSex[,1]|>str_replace_all("FullGroup","")
# 
# ArgTTBGtab3$RefSex<-str_split_fixed(ArgTTBGtab3$Hypothesis, boundary("word"),5)[,3]|>
#   str_split_fixed("[_]",2)
# ArgTTBGtab3$RefSex<-ArgTTBGtab3$RefSex[,1]|>str_replace_all("FullGroup","")
# ArgTTBGtab3$RefSex[ArgTTBGtab3$RefSex=="0"]<-"Female"
# 
# ArgTTBGtab3$FullGroupRef=paste(ArgTTBGtab3$RefSex,ArgTTBGtab3$Site,sep="_")
# ArgTTBGtab3$FullGroupContr=paste(ArgTTBGtab3$ContrSex,ArgTTBGtab3$Site,sep="_")
# 
# ArgTTBGtab3$Time<-str_split_fixed(ArgTTBGtab3$Hypothesis, boundary("word"),5)[,1]|>
#   str_split_fixed("[:]",2)
# ArgTTBGtab3$Time<-ArgTTBGtab3$Time[,2]|>str_replace_all("TimeNom","")
# ArgTTBGtab3$Time<-as.numeric(ArgTTBGtab3$Time)
# 
# SuppTab17<-ArgTTBGtab3 |> select(c(Estimate,CI.Lower,CI.Upper,
#                                    Evid.Ratio,RefSex,ContrSex,Time,Site,
#                                    FullGroupRef,FullGroupContr))|>
#   mutate(Evid.Ratio=case_when(ArgTTBGtab3$Evid.Ratio<=10^(3/2)~"Strong",
#                               ArgTTBGtab3$Evid.Ratio<=10^(2)~"Very strong",
#                               ArgTTBGtab3$Evid.Ratio>10^(2)~"Decisive")) |>
#   group_by(Site) |> 
#   arrange(Time,ContrSex,RefSex)|>
#   gt() |>
#   # row_group_order(c("Female - FP","Female - KC","Female - TC",
#   #                   "Male - FP","Male - KC", "Male - TC"))|>
#   tab_header(title = md("Supplementary Table 17. Blood Glucose Changes During ArgTT Compared Within Group, Between Sex")) |>
#   tab_spanner(
#     label = "Sex Comparisons",
#     columns = c(ContrSex,RefSex)
#   )|>
#   tab_spanner(
#     label = "95% Credibility Interval (mM)",
#     columns = c(CI.Lower,CI.Upper)
#   )|>
#   cols_move_to_start(
#     columns = c(Site,ContrSex,RefSex,Time)
#   )|>
#   cols_hide(columns=c(FullGroupRef,FullGroupContr))|>
#   cols_label(CI.Lower = "Lower Bound",
#              CI.Upper = "Upper Bound",
#              Evid.Ratio = "Strength of Evidence",
#              RefSex = "Reference",
#              ContrSex = "Contrast",
#              Estimate = html("Difference in<br>Δ Blood Glucose (mM)"),
#              Time = html("Time<br>Post-Arginine (min)")
#   )|>
#   fmt_number(
#     columns = c(Estimate,CI.Lower,CI.Upper),
#     decimals = 1,
#     use_seps = TRUE
#   ) |>
#   # data_color( 
#   #   columns = c(RefSite), 
#   #   colors = scales::col_factor( 
#   #     palette = GroupPalette,
#   #     domain = levels(weekly_monitoring$FullGroup)
#   #   )) |>
#   # data_color(
#   #   columns = c(ContrSite),
#   #   colors = scales::col_factor(
#   #     palette = GroupPalette,
#   #     domain =  levels(weekly_monitoring$FullGroup)
# #   )) |>
# tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.<br>Results reported at time 0 are absolute differences in fasted blood glucose between groups;<br>other results are the differences between groups in the change in blood glucose from the indicated<br>time compared to time 0.")) |>
#   tab_style(
#     locations = cells_column_labels(columns = everything()),
#     style     = list(
#       cell_borders(sides = "bottom", weight = px(3)),
#       cell_text(weight = "bold")
#     )
#   ) |>
#   tab_style(
#     locations = cells_title(groups = "title"),
#     style     = list(
#       cell_text(weight = "bold", size = 24)
#     )
#   ) |>
#   tab_options(
#     table.font.style = "Arial",
#     table.font.color = "black"
#   )
# SuppTab17
# 
# gtsave(SuppTab17,"SuppTab17.pdf",path="./Tables")

#######Supplementary Table 16-18: ArgTT Ins######
Instab1<-filter(Instab1,Star=="*")
Instab1$Sex<-ifelse(grepl("Female",Instab1$Hypothesis),"Female","Male") 
Instab1$ContrSite<-str_extract(Instab1$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
Instab1$ContrSite[is.na(Instab1$ContrSite)]<-"FP"

Instab1$RefSite<-str_split_fixed(Instab1$Hypothesis, "-\\(Intercept",2)[,2]|>
  str_extract("FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
Instab1$RefSite[is.na(Instab1$RefSite)]<-"FP"

Instab1$Time<-str_extract(Instab1$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
Instab1$Time[is.na(Instab1$Time)]<-0

SuppTab16<-Instab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                               Evid.Ratio,Sex,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(Instab1$Evid.Ratio<=10^(3/2)~"Strong",
                              Instab1$Evid.Ratio<=10^(2)~"Very strong",
                              Instab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 16. Human Insulin During ArgTT Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Insulin (nM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )) 

SuppTab16
#is there a text wrapping option??

gtsave(SuppTab16,"SuppTab16.pdf",path="./Tables")


Instab2<-filter(Instab2,Star=="*")
Instab2$Sex<-str_extract(Instab2$Hypothesis,"FullGroup(Female|Male)")|>
  str_replace_all("FullGroup","")
Instab2$Sex[is.na(Instab2$Sex)]<-"Female"

Instab2$Site<-str_extract(Instab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
Instab2$Site[is.na(Instab2$Site)]<-"FP"

Instab2$Time<-str_extract(Instab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()


SuppTab17<-Instab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                               Evid.Ratio,Sex,Site,Time))|>
  mutate(Evid.Ratio=case_when(Instab2$Evid.Ratio<=10^(3/2)~"Strong",
                              Instab2$Evid.Ratio<=10^(2)~"Very strong",
                              Instab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Time)|>
  gt() |>
  row_group_order(c("Female - KC","Female - TC", #"Male - FP",#"Female - FP","Male - KC",
                    "Male - TC"))|>
  tab_header(title = md("Supplementary Table 17. Human Insulin During ArgTT Compared to Fasting")) |>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>Insulin (nM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("ArgTT=Arginine Tolerance Test.<br>If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab17

gtsave(SuppTab17,"SuppTab17.pdf",path="./Tables")

Instab3<-filter(Instab3,Star=="*")
Instab3$Site<-str_extract(Instab3$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")

Instab3$ContrSex<-"Female"
Instab3$RefSex<-"Male"

Instab3$Time<-str_extract(Instab3$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
Instab3$Time[is.na(Instab3$Time)]<-0

SuppTab18<-Instab3 |> select(c(Estimate,CI.Lower,CI.Upper,
                               Evid.Ratio,RefSex,ContrSex,Time,Site))|>
  mutate(Evid.Ratio=case_when(Instab3$Evid.Ratio<=10^(3/2)~"Strong",
                              Instab3$Evid.Ratio<=10^(2)~"Very strong",
                              Instab3$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Site) |> 
  arrange(Time,ContrSex,RefSex)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 18. Human Insulin During ArgTT Compared Within Group, Between Sex")) |>
  tab_spanner(
    label = "Sex Comparisons",
    columns = c(ContrSex,RefSex)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (nM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  cols_move_to_start(
    columns = c(Site,ContrSex,RefSex,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSex = "Reference",
             ContrSex = "Contrast",
             Estimate = html("Difference in<br>Insulin (nM)"),
             Time = html("Time<br>Post-Arginine (min)")
  )|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab18

gtsave(SuppTab18,"SuppTab18.pdf",path="./Tables")

#######Supplementary Table 19-21: ArgTT GLP1######
GLP1tab1<-filter(GLP1tab1,Star=="*")
GLP1tab1$Sex<-ifelse(grepl("Female",GLP1tab1$Hypothesis),"Female","Male") 
GLP1tab1$ContrSite<-str_extract(GLP1tab1$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GLP1tab1$ContrSite[is.na(GLP1tab1$ContrSite)]<-"FP"

GLP1tab1$RefSite<-str_split_fixed(GLP1tab1$Hypothesis, "-\\(Intercept",2)[,2]|>
  str_extract("FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GLP1tab1$RefSite[is.na(GLP1tab1$RefSite)]<-"FP"

GLP1tab1$Time<-str_extract(GLP1tab1$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
GLP1tab1$Time[is.na(GLP1tab1$Time)]<-0

SuppTab19<-GLP1tab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                               Evid.Ratio,Sex,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(GLP1tab1$Evid.Ratio<=10^(3/2)~"Strong",
                              GLP1tab1$Evid.Ratio<=10^(2)~"Very strong",
                              GLP1tab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 19. Human or Mouse GLP-1 During ArgTT Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (pM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>GLP-1 (pM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )) 

SuppTab19
#is there a text wrapping option??

gtsave(SuppTab19,"SuppTab19.pdf",path="./Tables")


GLP1tab2<-filter(GLP1tab2,Star=="*")
GLP1tab2$Sex<-str_extract(GLP1tab2$Hypothesis,"FullGroup(Female|Male)")|>
  str_replace_all("FullGroup","")
GLP1tab2$Sex[is.na(GLP1tab2$Sex)]<-"Female"

GLP1tab2$Site<-str_extract(GLP1tab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
GLP1tab2$Site[is.na(GLP1tab2$Site)]<-"FP"

GLP1tab2$Time<-str_extract(GLP1tab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()


SuppTab20<-GLP1tab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                               Evid.Ratio,Sex,Site,Time))|>
  mutate(Evid.Ratio=case_when(GLP1tab2$Evid.Ratio<=10^(3/2)~"Strong",
                              GLP1tab2$Evid.Ratio<=10^(2)~"Very strong",
                              GLP1tab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Time)|>
  gt() |>
  row_group_order(c("Female - KC","Female - TC", #"Male - FP",#"Female - FP","Male - KC",
                    "Male - TC"))|>
  tab_header(title = md("Supplementary Table 20. Human or Mouse GLP-1 During ArgTT Compared to Fasting")) |>
  tab_spanner(
    label = "95% Credibility Interval (pM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>GLP-1 (pM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("ArgTT=Arginine Tolerance Test.<br>If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab20

gtsave(SuppTab20,"SuppTab20.pdf",path="./Tables")

GLP1tab3<-filter(GLP1tab3,Star=="*")
GLP1tab3$Site<-str_extract(GLP1tab3$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")

GLP1tab3$ContrSex<-"Female"
GLP1tab3$RefSex<-"Male"

GLP1tab3$Time<-str_extract(GLP1tab3$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
GLP1tab3$Time[is.na(GLP1tab3$Time)]<-0

SuppTab21<-GLP1tab3 |> select(c(Estimate,CI.Lower,CI.Upper,
                               Evid.Ratio,RefSex,ContrSex,Time,Site))|>
  mutate(Evid.Ratio=case_when(GLP1tab3$Evid.Ratio<=10^(3/2)~"Strong",
                              GLP1tab3$Evid.Ratio<=10^(2)~"Very strong",
                              GLP1tab3$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Site) |> 
  arrange(Time,ContrSex,RefSex)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 21. Human GLP-1 During ArgTT Compared Within Group, Between Sex")) |>
  tab_spanner(
    label = "Sex Comparisons",
    columns = c(ContrSex,RefSex)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (pM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  cols_move_to_start(
    columns = c(Site,ContrSex,RefSex,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSex = "Reference",
             ContrSex = "Contrast",
             Estimate = html("Difference in<br>GLP-1 (pM)"),
             Time = html("Time<br>Post-Arginine (min)")
  )|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab21

gtsave(SuppTab21,"SuppTab21.pdf",path="./Tables")

#######Supplementary Table 22-25: ArgTT Gcg######
Gcgtab1<-filter(Gcgtab1,Star=="*")
Gcgtab1$Sex<-ifelse(grepl("Female",Gcgtab1$Hypothesis),"Female","Male") 
Gcgtab1$ContrSite<-str_extract(Gcgtab1$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
Gcgtab1$ContrSite[is.na(Gcgtab1$ContrSite)]<-"FP"

Gcgtab1$RefSite<-str_split_fixed(Gcgtab1$Hypothesis, "-\\(Intercept",2)[,2]|>
  str_extract("FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
Gcgtab1$RefSite[is.na(Gcgtab1$RefSite)]<-"FP"

Gcgtab1$Time<-str_extract(Gcgtab1$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
Gcgtab1$Time[is.na(Gcgtab1$Time)]<-0

SuppTab22<-Gcgtab1 |> select(c(Estimate,CI.Lower,CI.Upper,
                                Evid.Ratio,Sex,Time,RefSite,ContrSite))|>
  mutate(Evid.Ratio=case_when(Gcgtab1$Evid.Ratio<=10^(3/2)~"Strong",
                              Gcgtab1$Evid.Ratio<=10^(2)~"Very strong",
                              Gcgtab1$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex) |> 
  arrange(Time,ContrSite,RefSite)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 22. Human or Mouse Glucagon During ArgTT Compared Within Sex, Between Groups")) |>
  tab_spanner(
    label = "Site Comparisons",
    columns = c(ContrSite,RefSite)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (pM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,ContrSite,RefSite,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSite = "Reference",
             ContrSite = "Contrast",
             Estimate = html("Difference in<br>Glucagon (pM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )) 

SuppTab22
#is there a text wrapping option??

gtsave(SuppTab22,"SuppTab22.pdf",path="./Tables")


Gcgtab2<-filter(Gcgtab2,Star=="*")
Gcgtab2$Sex<-str_extract(Gcgtab2$Hypothesis,"FullGroup(Female|Male)")|>
  str_replace_all("FullGroup","")
Gcgtab2$Sex[is.na(Gcgtab2$Sex)]<-"Female"

Gcgtab2$Site<-str_extract(Gcgtab2$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")
Gcgtab2$Site[is.na(Gcgtab2$Site)]<-"FP"

Gcgtab2$Time<-str_extract(Gcgtab2$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()


SuppTab23<-Gcgtab2 |> select(c(Estimate,CI.Lower,CI.Upper,
                                Evid.Ratio,Sex,Site,Time))|>
  mutate(Evid.Ratio=case_when(Gcgtab2$Evid.Ratio<=10^(3/2)~"Strong",
                              Gcgtab2$Evid.Ratio<=10^(2)~"Very strong",
                              Gcgtab2$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Sex,Site) |> 
  arrange(Time)|>
  gt() |>
  row_group_order(c("Female - KC","Female - TC", #"Male - FP",#"Female - FP","Male - KC",
                    "Male - TC"))|>
  tab_header(title = md("Supplementary Table 23. Human or Mouse Glucagon During ArgTT Compared to Fasting")) |>
  tab_spanner(
    label = "95% Credibility Interval (pM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  # tab_spanner(
  #   label = "Genotype Comparisons",
  #   columns = c(Comp2,Comp1)
  # )|>
  cols_move_to_start(
    columns = c(Sex,Site,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             # Comp2 = "Reference",
             # Comp1 = "Contrast",
             Estimate = html("Difference in<br>Glucagon (pM)"),
             Time = html("Time<br>Post-Arginine (min)"))|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  tab_source_note(md("ArgTT=Arginine Tolerance Test.<br>If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab23

gtsave(SuppTab23,"SuppTab23.pdf",path="./Tables")

Gcgtab3<-filter(Gcgtab3,Star=="*")
Gcgtab3$Site<-str_extract(Gcgtab3$Hypothesis,"FullGroup(Female|Male)[_][A-Z]{2}")|>
  str_replace_all("FullGroupFemale_","")|>
  str_replace_all("FullGroupMale_","")

Gcgtab3$ContrSex<-"Female"
Gcgtab3$RefSex<-"Male"

Gcgtab3$Time<-str_extract(Gcgtab3$Hypothesis,"TimeNom[0-9]{1,3}")|>
  str_replace_all("TimeNom","")|>as.numeric()
Gcgtab3$Time[is.na(Gcgtab3$Time)]<-0

SuppTab24<-Gcgtab3 |> select(c(Estimate,CI.Lower,CI.Upper,
                                Evid.Ratio,RefSex,ContrSex,Time,Site))|>
  mutate(Evid.Ratio=case_when(Gcgtab3$Evid.Ratio<=10^(3/2)~"Strong",
                              Gcgtab3$Evid.Ratio<=10^(2)~"Very strong",
                              Gcgtab3$Evid.Ratio>10^(2)~"Decisive")) |>
  group_by(Site) |> 
  arrange(Time,ContrSex,RefSex)|>
  gt() |>
  # row_group_order(c("Female - FP","Female - KC","Female - TC",
  #                   "Male - FP","Male - KC", "Male - TC"))|>
  tab_header(title = md("Supplementary Table 24. Human Glucagon During ArgTT Compared Within Group, Between Sex")) |>
  tab_spanner(
    label = "Sex Comparisons",
    columns = c(ContrSex,RefSex)
  )|>
  tab_spanner(
    label = "95% Credibility Interval (pM)",
    columns = c(CI.Lower,CI.Upper)
  )|>
  cols_move_to_start(
    columns = c(Site,ContrSex,RefSex,Time)
  )|>
  cols_label(CI.Lower = "Lower Bound",
             CI.Upper = "Upper Bound",
             Evid.Ratio = "Strength of Evidence",
             RefSex = "Reference",
             ContrSex = "Contrast",
             Estimate = html("Difference in<br>Glucagon (pM)"),
             Time = html("Time<br>Post-Arginine (min)")
  )|>
  fmt_number(
    columns = c(Estimate,CI.Lower,CI.Upper),
    decimals = 1,
    use_seps = TRUE
  ) |>
  # data_color( 
  #   columns = c(RefSite), 
  #   colors = scales::col_factor( 
  #     palette = GroupPalette,
  #     domain = levels(weekly_monitoring$FullGroup)
  #   )) |>
  # data_color(
  #   columns = c(ContrSite),
  #   colors = scales::col_factor(
  #     palette = GroupPalette,
  #     domain =  levels(weekly_monitoring$FullGroup)
#   )) |>
tab_source_note(md("ArgTT=Arginine Tolerance Test. If no result is reported, the posterior probability did not exceed 95%.")) |>
  tab_style(
    locations = cells_column_labels(columns = everything()),
    style     = list(
      cell_borders(sides = "bottom", weight = px(3)),
      cell_text(weight = "bold")
    )
  ) |>
  tab_style(
    locations = cells_title(groups = "title"),
    style     = list(
      cell_text(weight = "bold", size = 24)
    )
  ) |>
  tab_options(
    table.font.style = "Arial",
    table.font.color = "black"
  )
SuppTab24

gtsave(SuppTab24,"SuppTab24.pdf",path="./Tables")
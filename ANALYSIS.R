############## SUBSISTENCE MODEL ################

## RELEASE 1.0 NOTES ###


### Step 2-4 include processes used to run paleoproductivity model.
 ## Spatial data used in these steps (catchment.shp, drainage_area_high.shp, drainage_area_moderate.shp, drinage_area_low.shp, and mesa_top_area.shp)
 ## do not represent actual locations of farming areas and catchments used in the study, due to 
 ## site sensitivity data. Instead, these data are "invented" to show how the model functioned and were 
 ## are located in publicly accessible areas of the Far View community in Mesa Verde National Park. 
 ## Climate data used in these steps were downloaded at https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/usa/bocinsky2016/
 ## and are derived from Paleocar (https://github.com/bocinsky/paleocar)

### Step 5-6 include processes for determining surpluses/deficits.
 ## Data needed for these steps are complete tabular results from the model. 

### Step 7-8 include processes for producing figures 5-8 from the associated publication.
 ## Data needed for these steps are complete tabular results from the model.

## END RELEASE 1.0 NOTES ###



######### STEP 0.0: # ENVIRONMENTS & DATA ######################################
packages <-c('sf','terra','ggplot2','dplyr',
             'gdistance','tidyr','zoo','ggpubr','ggridges','tidyterra')
for(p in packages) if(p %in% rownames(installed.packages()) == F) { install.packages(p) }
for(p in packages) suppressPackageStartupMessages(library(p,quietly=T,character.only=T))

setwd("C:/Users/sfield3/OneDrive - University of Wyoming/RESEARCH/PROJECT/PALEOPROD_FARVIEW")
theme_set(theme_bw())



######### STEP 1.0: IMPORT FUNCTIONS & DATA #####################################
## Import paleoclimate rasters
precip_full <- terra::rast("./DATA/GITHUB DATA/PPT_1-2000.tif")
crs(precip_full)<- "epsg:4326"

niche_full <- terra::rast("./DATA/GITHUB DATA/niche_1-2000.tif")
crs(niche_full) <- "epsg:4326"

## Import community catchment boundary
bound <- sf::read_sf("./DATA/GITHUB DATA/catchment.shp")%>%
  st_transform(4326)

## Import potential farming areas
mt <- sf::read_sf("./DATA/GITHUB DATA/mesa_top_area.shp")%>%
  st_transform(4326)

d_high <- sf::read_sf("./DATA/GITHUB DATA/drainage_area_high.shp")%>%
  st_transform(4326)

d_mod <- sf::read_sf("./DATA/GITHUB DATA/drainage_area_moderate.shp")%>%
  st_transform(4326)

d_low <- sf::read_sf("./DATA/GITHUB DATA/drainage_area_low.shp")%>%
  st_transform(4326)


## Import productivity and yield rates 
norm_cd <- read.csv("./DATA/GITHUB DATA/check_dam_success_distribution.csv",header=T,fileEncoding = 'UTF-8-BOM')
norm_mt <- read.csv("./DATA/GITHUB DATA/mesatop_success_distributiontotalarea.csv",header=T,fileEncoding = 'UTF-8-BOM')
yields <- read.csv("./DATA/GITHUB DATA/yield_summary_table_pueblofarmingproject.csv",header=T,fileEncoding = 'UTF-8-BOM')


## Import demographic reconstruction
demographic <- read.csv("https://raw.githubusercontent.com/sfield2/Acute-Climate-Stress/main/DATA/all_community_occupations.csv")%>%
  .[,1:4]



############# STEP 2: Structure data for model #################
yields <- subset (yields, PFP.experimental.yield..kg.ha. > 0)
yields_cd <- subset(yields, Garden =="CDG")
yields_mt <- subset(yields, Garden !="CDG")


# smooth population curve for better idea of pop needs
pop.lo <- loess(fv_pop_av~Year, data=demographic, span=0.2)
lo.val <- predict(pop.lo)

demographic<-demographic %>%
  mutate(fv_pop_smooth = round(lo.val))%>%
  mutate(fv_pop_smooth = ifelse(fv_pop_smooth < 0 , 0, fv_pop_smooth))

rm(pop.lo,lo.val)



############# STEP 3: Build output for the productivity model #################
n <- nrow(demographic)
firstyear <- min(demographic$Year)
lastyear <- max(demographic$Year)

results_simulation <- as.data.frame(matrix(NA,n,3))%>%
  `colnames<-`(c("Year", "mesa_top_productivity", "drainage_productivity"))%>%
  mutate(Year=firstyear:lastyear)


############# STEP 4: Run the model #################
for (i in firstyear:lastyear) {
  
  # read in that years paleo-climate reconstructions
  niche <- niche_full[[i]]
  precip <- precip_full[[i]]

  # clip climate reconstructions by community boundary
  niche_comm <- terra::crop(niche, bound)

  # build polygon that encompasses in-niche areas
  niche_t <- as.polygons(niche_comm)
  names(niche_t) <- c("inorout")
  
  in_niche <- subset(niche_t, niche_t$inorout == 1)%>%
    sf::st_as_sf()
  
  # determine how much of the catchment is in the niche 
  comm_niche <- st_intersection(bound,in_niche)

  # determine how much of the mesatop area is in the niche 
  mt_tot<- st_intersection(comm_niche,mt)
  
  # find mesa top area (in HA) and assess if its zero
  mt_area<- as.numeric(st_area(mt_tot)/10000)
  
  if(length(mt_area)==0){
    mt_area = 0
  } else {
    mt_area <- mt_area
  }
  
  # sample a productivity rate for mesa top areas
  samp_area <- sample((norm_mt$range),1)
  
  #sample a yield rate for productive mesa top areas
  samp_yield <- sample((yields_mt$PFP.experimental.yield..kg.ha.),1)
  
  # downscale mesa top area, multiply by yield rate, and transform into cal/kg
  results_simulation[(i-699),2]<- (((mt_area * samp_area) * samp_yield) * 3500)
  
  
  # calculate how much precipitation in community
  precip_comm <- terra::crop(precip, bound)
  precip_ave <- as.numeric(terra::global(precip_comm,fun="mean"))/25.4
  
  #determine which check dam network to use 
  if (precip_ave>20) {
    drainage <- d_high
  } else if (precip_ave <= 20 & precip_ave > 15) {
    drainage <- d_mod
  } else {
    drainage <- d_low
  }
  
  
  drainage <- st_intersection(drainage,comm_niche)
  drainage_df <- as.data.frame(drainage)
  drainage_area <-(nrow(drainage_df)*30)/10000
    
  # sample a productivity rate for drainage areas
  samp_area <- sample((norm_cd$range),1)
  
  #sample a yield rate for productive drainage areas
  samp_yield <- sample((yields_cd$PFP.experimental.yield..kg.ha.),1)
  
  # downscale drainage area, multiply by yield rate, and transform into cal/kg
  results_simulation[(i-699),3]<- (((drainage_area * samp_area) * samp_yield) * 3500)

}





############# STEP 5: Determine surplus or deficits #################
# read in model results
results_total_summary <- read.csv("./DATA/GITHUB DATA/results_total_summary.csv",header=T,fileEncoding = 'UTF-8-BOM')%>%
  mutate(total_productivity = mesa_top_productivity + drainage_productivity)%>%
  # remove what was used for seeds
  mutate(difference = (total_productivity * 0.75) - community_needs)

## build accumulation series
accumulated_cal <- as.data.frame(results_total_summary[,c("Year","difference")])%>%
  mutate(cal_remain_lag_one = ifelse(difference > 0, difference*0.90 , 0)) %>%
  mutate(cal_remain_lag_two = ifelse(cal_remain_lag_one > 0, difference*0.65, 0))%>%
  mutate(cal_remain_lag_three = ifelse(cal_remain_lag_two > 0, difference*0.3,0))%>%
  mutate(accumulated_cal = 0)


# run an accumulation model 
years <- nrow(results_total_summary)

for (i in 1:years){
  
  if (i == 1){
    accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")]
    
  } else if (i == 2){
    accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")]
    
  } else if (i == 3){
    accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")] + accumulated_cal[i-2,c("cal_remain_lag_two")]
    
  } else {
    accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")] + accumulated_cal[i-2,c("cal_remain_lag_two")]  + accumulated_cal[i-3,c("cal_remain_lag_three")] 
  }
}

deficit_years <- subset(accumulated_cal, accumulated_cal < 0 )




############# STEP 6: Determine how much production to offset deficit #################
deficit_years <- read.csv("./DATA/GITHUB DATA/Full Model Results/deficit_years.csv",header=T,fileEncoding = 'UTF-8-BOM')


increase_sim <- as.data.frame(matrix(0,50,3))%>%
  setNames(c("increase_amount","in_deficit_count","out_deficit_count"))

n_sim <- nrow(increase_sim)
for (a in 1:n_sim){
  
  increase_amt <- 1 +(a/n_sim)
  
  results_total_summary_increase <- read.csv("./DATA/GITHUB DATA/results_total_summary_intensification.csv",header=T,fileEncoding = 'UTF-8-BOM')%>%
    mutate(total_cal_up = ((mesatop_area * increase_amt) * (cal_per_mt_ha * increase_amt)) + 
             ((drainage_area * increase_amt) * (cal_per_drainage_ha * increase_amt)))%>%
    # remove what was used for seeds
    mutate(difference = (total_cal_up * 0.75) - community_needs)
  
  
  accumulated_cal <- as.data.frame(results_total_summary_increase[,c("Year","difference")])%>%
    mutate(cal_remain_lag_one = ifelse(difference > 0, difference*0.90 , 0)) %>%
    mutate(cal_remain_lag_two = ifelse(cal_remain_lag_one > 0, difference*0.65, 0))%>%
    mutate(cal_remain_lag_three = ifelse(cal_remain_lag_two > 0, difference*0.3,0))%>%
    mutate(accumulated_cal = 0)
  
  # run an accumulation model 
  years <- nrow(results_total_summary_increase)
  
  for (i in 1:years){
    
    if (i == 1){
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")]
      
    } else if (i == 2){
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")]
      
    } else if (i == 3){
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")] + accumulated_cal[i-2,c("cal_remain_lag_two")]
      
    } else {
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")] + accumulated_cal[i-2,c("cal_remain_lag_two")]  + accumulated_cal[i-3,c("cal_remain_lag_three")] 
    }
  }
  
  deficit_up <- subset(accumulated_cal, accumulated_cal < 0 )
  
  # add to deficit_years
  deficit_years_summary <- deficit_years%>%
    .[,c("Year","accumulated_cal")]%>%
    mutate(deficit = T)%>%
    mutate(up_deficit = Year %in% deficit_up$Year)
  
  
  results <-as.data.frame(table(deficit_years_summary$up_deficit))
  increase_sim[a,c("increase_amount")] <-increase_amt
  increase_sim[a,c("in_deficit_count")] <- results[2,2]
  increase_sim[a,c("out_deficit_count")] <- results[1,2]
  
}
  



############# STEP 6: Identify persistent deficits  #################

# years with marginal and/or persistent deficits (persistent could not be avoided with 20% increase in farming area & yields)
increase_amt <- 1.2

results_total_summary_increase <- read.csv("./DATA/GITHUB DATA/results_total_summary_intensification.csv",header=T,fileEncoding = 'UTF-8-BOM')%>%
  mutate(total_cal_up = ((mesatop_area * increase_amt) * (cal_per_mt_ha * increase_amt)) + 
           ((drainage_area * increase_amt) * (cal_per_drainage_ha * increase_amt)))%>%
  # remove what was used for seeds
  mutate(difference = (total_cal_up * 0.75) - community_needs)


  accumulated_cal <- as.data.frame(results_total_summary_increase[,c("Year","difference")])%>%
    mutate(cal_remain_lag_one = ifelse(difference > 0, difference*0.90 , 0)) %>%
    mutate(cal_remain_lag_two = ifelse(cal_remain_lag_one > 0, difference*0.65, 0))%>%
    mutate(cal_remain_lag_three = ifelse(cal_remain_lag_two > 0, difference*0.3,0))%>%
    mutate(accumulated_cal = 0)
  
  # run an accumulation model 
  years <- nrow(results_total_summary_increase)
  
  for (i in 1:years){
    
    if (i == 1){
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")]
      
    } else if (i == 2){
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")]
      
    } else if (i == 3){
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")] + accumulated_cal[i-2,c("cal_remain_lag_two")]
      
    } else {
      accumulated_cal[i,c("accumulated_cal")] <- accumulated_cal[i,c("difference")] + accumulated_cal[i-1,c("cal_remain_lag_one")] + accumulated_cal[i-2,c("cal_remain_lag_two")]  + accumulated_cal[i-3,c("cal_remain_lag_three")] 
    }
  }

  
persistent_deficit <- subset(accumulated_cal, accumulated_cal < 0 )


deficit_years <- deficit_years%>%
  mutate(persistent_deficit = Year %in% persistent_deficit$Year)

rm(list=ls())






### STEP 7: Read in data for manuscript figures
results_total_summary <- read.csv("./DATA/GITHUB DATA/Full Model Results/results_total_summary.csv",header=T,fileEncoding = 'UTF-8-BOM')
results_total <- read.csv("./DATA/GITHUB DATA/Full Model Results/results_total.csv",header=T,fileEncoding = 'UTF-8-BOM')
accumulated <- read.csv("./DATA/GITHUB DATA/Full Model Results/results_model_accumulated.csv",header=T,fileEncoding = 'UTF-8-BOM')
results_990 <- read.csv("./DATA/GITHUB DATA/Full Model Results/results_990-999.csv",header=T,fileEncoding = 'UTF-8-BOM')
climate <- read.csv("./DATA/GITHUB DATA/Full Model Results/results_climate.csv",header=T,fileEncoding = 'UTF-8-BOM')
deficit_years <- read.csv("./DATA/GITHUB DATA/Full Model Results/deficit_years.csv",header=T,fileEncoding = 'UTF-8-BOM')
increase_sim <- read.csv("./DATA/GITHUB DATA/Full Model Results/avoiding_deficits.csv",header=T,fileEncoding = 'UTF-8-BOM')

# a little summarizing
not_in_niche <- subset(climate,In_Niche == 0)
accumulated <- accumulated%>%
  mutate(in_out = ifelse(accumulated_cal >= 0,"Surplus","Deficit"))
persistent_deficits <- subset(deficit_years, persistent_deficit == "TRUE")




increase_sim_for_plot <- gather(increase_sim,"type","count",2:3)%>%
  mutate(pct_plot = as.numeric(round(count/129,2)))%>%
  mutate(pct_in = ifelse(type =="in_deficit_count",as.numeric(round(count/129,2)),NA))%>%
  mutate(pct_out = ifelse(type =="out_deficit_count",as.numeric(round(count/129,2)),NA))



### STEP 8: Figure 4
results_total_ex_non <- subset(results_total, Year %in% c("990","991","992","993","994","995","996","997","998","999"))


# niche percentages
p1 <- ggplot(results_990)+
  geom_bar(aes(Year,V1.1),stat="identity",fill="darkolivegreen4")+
  scale_x_continuous(limits=c(989.5,1000.8),breaks=c(990:999))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=0.5))+
  labs(title="[A]",
       x="",
       y= "Proportion of 
Catchment in Niche")+
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length=unit(.25, "cm"))

# precipitation
p2<- ggplot(results_990)+
  geom_bar(aes(Year,V1/10),stat="identity",fill="cadetblue4")+
  scale_x_continuous(limits=c(989.5,1000.8),breaks=c(990:999))+
  scale_y_continuous(limits=c(0,66),breaks=seq(0,65,by=10))+
  labs(title="[B]",
       x="",
       y= "Precip. 
(cm)")+
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length=unit(.25, "cm"))


# variation in arable area
p3<-ggplot(results_total_ex_non) +
  geom_hline(yintercept=4.2*25,col="cadetblue",size=1.2,alpha=0.4,linetype="dashed")+
  geom_hline(yintercept=121,col="darkolivegreen",size=1.2,alpha=0.4,linetype="dashed")+
  geom_jitter(aes(Year-0.1, mesatop_arable_area),size=3,col="darkolivegreen3",width=0.25,alpha=0.3)+
  stat_summary(aes(Year-0.1,mesatop_arable_area),fun = "mean", geom = "crossbar", size = .7,color="darkolivegreen4", width=0.4) +
  geom_jitter(aes(Year+0.1,drainage_arable_area*25),size=3,col="cadetblue3",width=0.25,alpha=0.3)+
  stat_summary(aes(Year+0.1,drainage_arable_area*25),fun = "mean", geom = "crossbar", size = .7,color="cadetblue4", width=0.4) +
  scale_x_continuous(limits=c(989.5,1000.8),breaks=c(990:999))+
  scale_y_continuous(name="Arable Mesa 
  Top Hectares",limits=c(0,130),breaks=seq(0,125,by=25),sec.axis=sec_axis(~./25,name="Arable 
  Drainage Hectares"))+
  
  
  annotate("text",x= 998.6, y = 113, label = "Maximum Potential Drainage Farming Area", size =6, fontface= "bold", col="cadetblue4")+
  
  annotate("text",x= 991.7, y = 129, label = "Maximum Potential Mesa Top Farming Area", size =6, fontface= "bold", col="darkolivegreen4")+
  
  labs(title="[C]",x="",
       y= "Arable Areas (HA)")+
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.title.y.right = element_text(size=26,color="black",face="bold", vjust = +1.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length=unit(.25, "cm"))


p4<-ggplot(results_total_ex_non)+
  geom_jitter(aes(Year-0.4,drainage_cal_produced/1000000),alpha=0.4,size=3,col="cadetblue3",width=0.2)+
  geom_jitter(aes(Year-0.2,mesatop_cal_produced/1000000),alpha=0.4,size=3,col="darkolivegreen3",width=0.2)+
  geom_jitter(aes(Year,total_cal_produced/1000000),alpha=0.4,size=3,col="mediumpurple2",width=0.2)+
  stat_summary(aes(Year-0.4,drainage_cal_produced/1000000),fun="mean",geom="crossbar", size = 0.7, color="cadetblue4",width=0.3)+
  stat_summary(aes(Year-0.2,mesatop_cal_produced/1000000),fun = "mean", geom = "crossbar", size = .7,color="darkolivegreen4", width=0.3) +
  stat_summary(aes(Year,total_cal_produced/1000000),fun = "mean", geom = "crossbar", size = .7,color="mediumpurple4", width=0.3)+
  coord_cartesian(ylim=c(0,400))+
  scale_x_continuous(limits=c(989.5,1000.8),breaks=c(990:999))+
  labs(title="[D]",
       x="Year (AD)",
       y= "Maize Production
(Millions of Calories)")+
  
  annotate("text",x= 1000.3, y = 100, label = "Simulated 
Drainage
Production", size =6, fontface= "bold", col="cadetblue3")+
  
  
  annotate("text",x= 1000.3, y = 200, label = "Simulated
Mesa Top
Production", size =6, fontface= "bold", col="darkolivegreen4")+
  
  
  annotate("text",x= 1000.3, y = 300, label = "Total
Simulated
Production", size =6, fontface= "bold", col="mediumpurple4")+
  
  theme(title=element_text(size=24,color="black"),
        axis.text=element_text(size=24,color="black"),
        axis.title=element_text(size=26,color="black",face="bold"),
        axis.ticks.length=unit(.25, "cm"))






fig_4<-ggarrange(p1,p2,p3,p4,
                heights = c(.75,.75,1,1.5),
                ncol=1, nrow=4,align="v")

png('./FIGURES/Figure 4.png',height=1500,width=1000)
fig_4
dev.off()





### STEP 8: Figure 5
p5<-ggplot(data=climate)+
  geom_line(aes(Year,In_Niche),col="black",size=0.75)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  coord_cartesian(ylim=c(0,1))+
  labs(title="[A]",
       x= "",
       y = "Proportion of 
Catchment in Niche")+
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length.x=unit(.25, "cm"))



p6<-ggplot(data=results_total)+
  geom_line(aes(Year,mesatop_arable_area,group=Trial),col="darkolivegreen3",alpha=0.04,size=0.2)+
  geom_line(data=results_total_summary,aes(Year,mesatop_arable_area_1),col="darkolivegreen3",alpha=0.95,size=0.95)+
  geom_line(data=results_total_summary,aes(Year,rollmean(mesatop_arable_area_1,10,fill=NA,align=c("center"))),col="darkolivegreen4",size=1.15,alpha=0.85)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  labs(title="[B]",
       x="",
       y = "Arable Mesa Top 
Area (HA)")+
  annotate("text", x = 1255, y = 105, col="darkolivegreen3",size=6.8, fontface =2,hjust=0,alpha=0.65,
           label = "Single Iteration")+
  annotate("text", x = 1255, y = 85, col="darkolivegreen3",size=7, fontface =2,hjust=0,
           label = "Annual Mean 
of All Iteration")+
  annotate("text", x = 1255, y = 60, col="darkolivegreen4",size=7, fontface =2,hjust=0,
           label = "Ten Year Mean 
of All Iteration")+
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length.x=unit(.25, "cm"))


p7<-ggplot(data=results_total)+
  geom_line(aes(Year,drainage_arable_area,group=Trial),col="cadetblue3",alpha=0.04,size=0.3)+
  geom_line(data=results_total_summary,aes(Year,drainage_arable_area_1),col="cadetblue3",alpha=0.85,size=0.75)+
  geom_line(data=results_total_summary,aes(Year,rollmean(drainage_arable_area_1,10,fill=NA,align=c("center"))),col="cadetblue4",size=1.15,alpha=0.85)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  labs(title="[C]",
       x="",
       y = "Arable Mesa Slope 
Area (HA)")+
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length.x=unit(.25, "cm"))


p8<-ggplot(data=results_total)+
  geom_line(aes(Year,total_cal_produced/1000000,group=Trial),col="mediumpurple",alpha=0.04,size=0.1)+
  geom_line(data=results_total_summary,aes(Year,total_cal_produced_1/1000000),col="mediumpurple",alpha=0.85,size=0.75)+
  geom_line(data=results_total_summary,aes(Year,rollmean(total_cal_produced_1/1000000,10,fill=NA,align=c("center"))),col="mediumpurple4",size=1.15,alpha=0.85)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  coord_cartesian(ylim=c(0,400))+
  labs(title="[D]",
       x= "Year (AD)",
       y = "Maize Production 
(Millions of Calories)")+
  theme(axis.text=element_text(size=24,color="black"),
        axis.title=element_text(size=26,color="black",face="bold"),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"))


fig_5<-ggarrange(p5,p6,p7,p8,
                    heights = c(.85,.85,.85,1),
                    ncol=1, nrow=4,align="v")


png('./FIGURES/Figure 5.png',height=1700,width=1400)
fig_5
dev.off()





### STEP 9: Figure 6
p9<-ggplot(data=results_total)+
  geom_line(data=results_total_summary,aes(Year,total_cal_produced_1/1000000),col="mediumpurple",alpha=0.85,size=0.75)+
  geom_line(data=results_total_summary,aes(Year,rollmean(total_cal_produced_1/1000000,10,fill=NA,align=c("center"))),col="mediumpurple4",size=1.15,alpha=0.85)+
  geom_line(data=results_total_summary,aes(Year,pop_cal_needs_1/1000000),col="#FAAA14",size=1.45,alpha=0.95)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limits=c(0,300),breaks=seq(0,300,by=100))+
  labs(title="[A]",
       y="Maize Production & Needs
(Millions of Calories)",
       x = "")+
  
  annotate("text", x = 1252, y = 150, col="mediumpurple",size=6.8, fontface =2,hjust=0,
           label = "Maize Production")+
  annotate("text", x = 1255, y = 35, col="#FAAA14",size=7, fontface =2,hjust=0,
           label = "Maize Needs")+
  
  
  theme(title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        panel.grid.minor.y = element_blank())


p10<-ggplot(data=results_total)+
  geom_line(data=accumulated,aes(Year,accumulated_cal/1000000,col=in_out,group=1),size=1.15,alpha=0.25)+
  geom_line(data=accumulated,aes(Year,rollmean(accumulated_cal/1000000,10,fill=NA,align=c("center")),col=in_out,group=1),size=1.45,alpha=0.85)+
  scale_color_manual(values=c("coral3","darkseagreen4"),name="",labels=c("In Deficit","Out of Deficit"))+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limits=c(-200,400),breaks=seq(-200,400,by=100))+
  
  annotate("text", x = 1255, y = 200, col="darkseagreen4",size=7, fontface =2,hjust=0,
           label = "Surplus Years")+
  annotate("text", x = 1162, y = -50, col="coral3",size=7, fontface =2,hjust=0,
           label = "Deficit Years")+
  
  
  labs(title="[B]",
       y="Maize Accumulated
(Millions of Calories)",
       x = "")+
  theme(legend.position="none",
        title=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        panel.grid.minor.y = element_blank())


p11<-ggplot(data=accumulated)+
  geom_col(data=deficit_years,aes(x=Year, y = 1),fill="salmon3",width=1,alpha=0.75)+
  labs(title = "[C]",
       x = "",
       y = "")+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  theme(axis.text.x=element_text(size=24,color="black"),
        axis.text.y=element_blank(),
        axis.title=element_text(size=26,color="black",face="bold"),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



fig_6<-ggarrange(p9,p10,p11,
                 heights = c(0.65,1,0.5),
                 ncol=1, nrow=3,align="v")


png('./FIGURES/Figure 6.png',height=1300,width=1500)
fig_6
dev.off()


### STEP 10: Figure 7
p12<-ggplot(increase_sim_for_plot,aes(fill=type,x=increase_amount,y = pct_plot))+
  geom_bar(stat="identity")+
  geom_text(aes(label=pct_in),y=1.05,col="grey30",size=6.5,fontface=2)+
  geom_text(aes(label=pct_out),y=-0.05,col="grey30",size=6.5,fontface=2)+
  scale_fill_manual(values=c("coral3","darkseagreen3"),name="",labels=c("Remaining In Deficit","Out of Deficit"))+
  coord_cartesian(ylim=c(-0.1,1.1))+
  scale_x_continuous(limits=c(1,1.65),breaks=seq(1,1.6,by=0.12),labels=c("0%","12%","24%","36%","48%","60%"))+
  labs(title="Proportion of Years In/Out of Deficit",
       y= "",
       x = "Increase in Farming Area & Productivity")+
  theme(legend.text=element_text(size=24, color="black",face="bold"),
        axis.text.x=element_text(size=24,color="black"),
        axis.title.x=element_text(size=26,color="black",face="bold"),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

png('./FIGURES/Figure 7.png',height=700,width=1900)
p12
dev.off() 





### STEP 11: Figure 8
pop <- read.csv("https://raw.githubusercontent.com/sfield2/Acute-Climate-Stress/main/DATA/all_community_occupations.csv")%>%
  .[,1:4]

pop[1,c("fv_pop_av")]<-0
pop[551,c("fv_pop_av")]<-0
pop[552:581,c("fv_pop_av")]<- NA

fv_stress_regime <- read.csv("https://raw.githubusercontent.com/sfield2/Acute-Climate-Stress/main/DATA/FVC_stress_regime.csv")
potential_stress_periods <- subset(fv_stress_regime, gen_av > "1")

precip <- read.csv("./OUTPUT/UPDATES/results_nichesize_precipamt.csv",header=T,fileEncoding = 'UTF-8-BOM')

low_precip <- subset(precip,V2<13.77)

sp <- NULL 
sapply(seq(nrow(potential_stress_periods)), function(x)
{
  ifelse((sum(diff(potential_stress_periods[x:(x+4), "Year"], 1)) == 4 &
            sum(diff(potential_stress_periods[x:(x+4), "Year"], 1) == 1) == 4),
         sp <<- rbind(sp, potential_stress_periods[x:(x+4),]),"")
})

fv_stress_periods <- unique(sp)%>%
  subset(Year <= 1250)




p8<-ggplot(data=accumulated)+
  geom_line(data=pop,aes(x=Year,y=fv_pop_av),col="black",size=1.5)+
  labs(title = "[A]",
       x = "Year (AD)",
       y = "Population")+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  theme(axis.text.y=element_text(size=24,color="black"),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=26,color="black",face="bold"),
        axis.title.x = element_blank(),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

p9<-ggplot()+
  geom_col(data=persistent_deficits,aes(x=Year, y = 1),fill="salmon3",width=2,alpha=0.75)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limit=c(0,1))+
  labs(title="[B]")+
  annotate("text", x = 1250, y = 0.8, col="grey30",size=7, fontface =2,hjust=0,
           label = "Persistent 
Deficit Years")+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

p10<-ggplot()+
  geom_col(data=low_precip,aes(x=Year, y = 1),fill="salmon3",width=2,alpha=0.75)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limit=c(0,1))+
  labs(title="[C]")+
  annotate("text", x = 1250, y = 0.8, col="grey30",size=7, fontface =2,hjust=0,
           label = "Low Precipitation 
Years")+
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


p11<-ggplot()+
  geom_col(data=fv_stress_periods,aes(x=Year, y = 1),fill="salmon3",width=2,alpha=0.75)+
  scale_x_continuous(limits=c(700,1300),breaks=seq(600,1300,by=50))+
  scale_y_continuous(limit=c(0,1))+
  annotate("text", x = 1250, y = 0.8, col="grey30",size=7, fontface =2,hjust=0,
           label = "Acute Stress 
Periods")+
  labs(title="[D]",
       x = "Year (AD)")+
  theme(axis.text.x=element_text(size=24,color="black"),
        axis.title.x=element_text(size=26,color="black",face="bold"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        title=element_text(size=24,color="black"),
        axis.ticks.length.x=unit(.25, "cm"),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())


fig_8<-ggarrange(p8,p9,p10,p11,
                 heights = c(1,0.5,0.5,0.6),
                 ncol=1, nrow=4,align="v")


png('./FIGURES/Figure 8.png',height=1300,width=1600)
fig_8
dev.off()




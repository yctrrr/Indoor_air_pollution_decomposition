setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
dat_dir <- "D:\\shaoyanchuan\\data\\"
library(data.table)
library(dplyr)
library(sp)
library(sf)
library(rgdal)
library(pinyin)
library(terra)
library(lubridate)
library(tidyr)
library(scales)
library(Cairo)
library(colorspace)

source("D:/shaoyanchuan/codebook/function/Sf_visulaization.R")
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('./script/Core.R', encoding = 'UTF8')
county_sf <- readRDS(paste0(dat_dir,"Shp\\county_China\\county_shp.Rds"))%>%
  dplyr::select(County,Province,City)
county_sp <- as_Spatial(county_sf)

##========Process activity pattern and infiltration factor========
##  Activity pattern
prov_info <- read.csv("result\\Table\\Input_data\\Province_information.csv", encoding = "UTF-8")
time_fraction <- fread("data\\Activity_pattern_EHB.csv")%>%
  mutate_at(.vars = vars(4:11), .funs = function(x) x/24)%>% ## calculate fraction
  melt(id.vars = c("Subarea","Condition","Activity"), variable.name = "Gender", value.name = "Time") ## wide to long

##  Infiltration data
inf <- fread("data\\Infiltration_factor_by_province.csv")

hs_infwide <- dplyr::select(inf,Winter.Mean:Winter.P95)%>%
  setnames(gsub("Winter","HS",colnames(.)))

hs_inf <- cbind(dplyr::select(inf,Province), hs_infwide)%>%
  melt("Province",variable.name = "Metric", value.name = "Infl")%>%
  filter(Metric%in%c("HS.Mean","HS.P5","HS.P95"))%>%
  mutate(Condition = "H")

##  Divide by heating and non-heating
nhs_inf_v <- vector(mode = "list", length = nrow(inf))
for (i in 1:nrow(inf)) {
  nhs_inf_sample <- c(rnorm(1000,inf$Summer.Mean[i],inf$Summer.SD[i]),
              rnorm(2000,inf$SA.Mean[i], inf$SA.SD[i]))
  nhs_inf_q <- quantile(nhs_inf_sample,c(.05,.25,.5,.75,.95))
  
  nhs_inft  <- c(mean(nhs_inf_sample),sd(nhs_inf_sample),nhs_inf_q)%>%
    t()%>%
    as.data.frame()%>%
    setnames(c("NHS.Mean","NHS.SD","NHS.P5","NHS.P25","NHS.50","NHS.P75","NHS.P95"))%>%
    mutate(across(where(is.numeric), ~ round(., 2)))
  
  nhs_inf_v[[i]] <- nhs_inft
}

nhs_infwide <- bind_rows(nhs_inf_v)
nhs_inf <- cbind(dplyr::select(inf,Province), nhs_infwide)%>%
  melt("Province", variable.name = "Metric", value.name = "Infl")%>%
  filter(Metric%in%c("NHS.Mean","NHS.P5","NHS.P95"))%>%
  mutate(Condition = "N-H")

## Match heating and non-heating infiltration data
inf_byh <- full_join(hs_inf, nhs_inf)%>%
  left_join(prov_info)

activ_info <- inner_join(inf_byh, time_fraction)%>%
  mutate(sex = Gender%>%str_extract_all("[aA-zZ]+")%>%tolower(),
         agediv = Gender%>%str_extract_all("[0-9]+")%>%as.numeric())%>%
  dplyr::select(-Gender,-whoname)

fwrite(activ_info, "result\\Table\\Input_data\\activity_information.csv")


##==========Integrated air pollution==========
activ_info <- fread("result\\Table\\Input_data\\activity_information.csv")%>%
  filter(Activity == "Outdoor", Metric%in%c("HS.Mean","NHS.Mean"))

IPM25_lst <- vector(mode = "list", length = length(2000:2020))
for (y in 2000:2020) {
  cat("=======Process integrated PM2.5 year: ",y,"=======\n")
  ydlen <- yday(as.Date(paste0(y,"-12-31")))
  yl <- y

  ## calculate indoor and outdoor air pollution weighted by time fraction and infiltration
  IPM25_byt <- rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",y,".tif"))%>%
    as.data.frame(xy = TRUE, na.rm =FALSE)%>%
    dplyr::select(x,y,APM25,hd_count,pindex,countyindex)%>%
    na.omit()%>%mutate(hd_f = hd_count/ydlen)%>%
    inner_join(activ_info)%>%
    mutate(AAP = APM25*Time , HAP_O =  APM25*(1-Time)*Infl,
           IAP_UCl = AAP + HAP_O,
           IAP_UB = 223*(1-Time)+ AAP + HAP_O,
           IAP_UC = 38*(1-Time)+ AAP + HAP_O, 
           IAP_RCl = AAP + HAP_O,
           IAP_RB = 250*(1-Time)+ AAP + HAP_O, 
           IAP_RC = 117*(1-Time)+ AAP + HAP_O) 
 
  ## weight integrated PM2.5 by heating and non-heating time
  IPM25_byh <- IPM25_byt%>%
    group_by(x,y,sex,agediv)%>%
    summarise_at(vars(AAP,HAP_O,IAP_UCl,IAP_UB,IAP_UC,IAP_RCl,IAP_RB,IAP_RC),
                 ~ .x[Condition == "H"]*hd_f[Condition == "H"] + 
                   .x[Condition == "N-H"]*(1-hd_f[Condition == "H"]))%>%
    ungroup()%>%mutate(year = yl)

  ## convert data frame to long format
  IPM25_grid <- 1:6%>%purrr::map_dfr(function(t){
    expo_type <- c("IAP_UCl","IAP_UB","IAP_UC","IAP_UCl","IAP_RB","IAP_RC")
    area_type <- c(rep("Urban",3),rep("Rural",3))
    fuel_type <- rep(c("Clean","Biomass","Coal"),2)
    
    IPM25_byh%>%dplyr::select(x,y,year,sex,agediv,AAP,HAP_O,concentration = !!expo_type[t])%>%
      mutate(HAP_I = concentration - AAP - HAP_O)%>%
      mutate(Area = area_type[t], Fuel = fuel_type[t])
    })
  
  fwrite(IPM25_grid, paste0("result\\Table\\Integrated_PM2.5_exposure\\Integrated_PM2.5_",y,".csv"))
}

## extract grid information
coords <- rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_2000.tif"))%>%
  as.data.frame(xy = TRUE, na.rm =FALSE)%>%
  dplyr::select(x,y,pindex,countyindex)%>%
  na.omit()
fwrite(coords, paste0("result\\Table\\Integrated_PM2.5_exposure\\Grid_information.csv"))


##=======Calculate population-weighted concentration========
year = 2000
PM_by_poplst = vector(mode = "list", length = length(2000:2020))
AAP_by_poplst = vector(mode = "list", length = length(2000:2020))
PM_poptable = expand_grid(concentration = seq(0, 300, by = 0.1)%>%round(1), year = 2000:2020)

fread(paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_",year,".csv"))

## population-weighted concentration
IPM_pop_group <- 2000:2020%>%
  purrr::map_dfr(function(year){
    
    cat("======Year: ",year,"======\n")
    yl = year
    PM_r <- fread(paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_",year,".csv"))%>%
      mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)

    AgeGroup <- fread('./data/IHME_GBD_2019_POP_SYA_2000-2020.csv')%>%
      mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                          labels = 1:4, include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
      mutate(sex = tolower(sex), agegroup = as.character(agegroup))%>%
      filter(year == yl)

    Pop_long <- fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))%>%
      filter(year == yl)%>%mutate(Pop = Pop*Popf, .keep = "unused")%>%
      mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)
    
    list(PM_r, AgeGroup, Pop_long)%>%reduce(left_join)%>%
      mutate(Pop = Pop*AgeStruc, Pop_all = sum(Pop), Popnf = Pop/Pop_all)%>%
      group_by(x,y,year,Area,Fuel)%>%
      summarise(concentration = weighted.mean(concentration, Pop),
                AAP = weighted.mean(AAP, Pop),
                HAP_O = weighted.mean(HAP_O, Pop),
                HAP_I = weighted.mean(HAP_I, Pop),
                Popnf = sum(Popnf), Pop = sum(Pop))%>%ungroup()%>%
      mutate_at(vars(concentration, AAP, HAP_O, HAP_I), round, digits = 1)
  })
fwrite(IPM_pop_group, "result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5_by_group.csv")


IPM_pop <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5_by_group.csv")%>%
  rename(IAP = concentration)%>%group_by(x, y, year)%>%
  summarise(IAP = weighted.mean(IAP, Pop)%>%replace_na(0),
            AAP = weighted.mean(AAP, Pop)%>%replace_na(0),
            HAP_O = weighted.mean(HAP_O, Pop)%>%replace_na(0),
            HAP_I = weighted.mean(HAP_I, Pop)%>%replace_na(0),
            Popnf = sum(Popnf), Pop = sum(Pop))%>%ungroup()%>%
  mutate_at(vars(IAP,AAP,HAP_O,HAP_I), round, digits = 1)

fwrite(IPM_pop, "result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5.csv")

## calculate exposed population fraction
IPM_pop <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5.csv")%>%
  mutate(AAP_ind = AAP + HAP_O, HAP_ind = HAP_I)

IPM_by_pop <- c("IAP","AAP_ind","HAP_ind")%>%
  purrr::map(~ IPM_pop%>%group_by_at(vars("year", .x))%>%
        summarise(Popnf = sum(Popnf), Pop = sum(Pop))%>%ungroup()%>%
        rename(concentration = .x)%>%left_join(PM_poptable,.)%>%
        replace_na(list(Popf = 0, Pop = 0))%>%mutate(Pop = round(Pop))%>%
        rename(!!paste0("Popnf_",.x):=Popnf, !!paste0("Pop_",.x):= Pop)
        )%>%reduce(left_join)
fwrite(IPM_by_pop, "result/Table/Integrated_PM2.5_exposure/PM2.5_group_by_population.csv")


##=======Visualize exposed population at grid scale========
EPop_plot <- function(data ,value = "Pop"){
  grid_plot(data = data, label = NULL, value = value)+
    scale_fill_manual(name = expression(paste("Fuel use population (#/ 100km"^"2",")")),
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.hjust=0.5,title.vjust=1,
                                           byrow = TRUE, ticks= F))+
    theme(panel.grid.major = element_line(colour = alpha("black", 0.4), size = 0.6))+
    geom_sf(data = polygon_sf, fill = "white", alpha = 0.1, size = 0.4)+
    geom_sf(data = Nline_sf, fill = "white", size = 0.6)+
    xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(strip.background = element_rect(color="black", fill="white", size = 0.5))+ ##change the facetwrap background
    theme(axis.text=element_text(size = 30))+
    theme(legend.title = element_text(colour = "black", size = 40))
}


polygon <- readOGR(dsn="data/Shp/Province.shp",
                   layer="Province",use_iconv=TRUE, encoding = "UTF-8")

Nline <- readOGR("data/Shp/Nline.shp",
                 layer="Nline",use_iconv=TRUE, encoding = "UTF-8")
polygon_sf <- st_as_sf(polygon)
Nline_sf <- st_as_sf(Nline)

IPM_spatial_pop <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5.csv")%>%
  mutate(Pop = ifelse(IAP < 35|is.na(IAP), 0, Pop%>%round()))

breaks <- c(-Inf, 10, 50, 100, 500, 1000, 5000, 10000, 20000, +Inf)
labels <- c("0 - 10", "10 - 50", "50 - 100","100 - 500","500 - 1000", 
            "1000 - 5000", "5000 - 10000", "10000 - 20000", "> 20000")
colors <- RColorBrewer::brewer.pal(length(labels),"YlGnBu")

year = 2020
for (year in 2000:2020) {
  yl <- year
  Pop_spatial <- IPM_spatial_pop%>%
    filter(year == yl)%>%
    dplyr::select(x,y,Pop)%>%
    dfproj()%>% mutate(Pop = cut(Pop, breaks, labels))

  Pop_spatial_plot <- EPop_plot(Pop_spatial)
  
  Cairo::CairoPNG(paste0("result/Fig/PM2.5_Exposure/Exposured_population_",year,".png"),
                  width = 2800,height = 1500, res=150)
  print(Pop_spatial_plot)
  dev.off()
}
##=======Visualize PM2.5 concentration at grid scale=======
IPM_plot <- function(data ,value = "IAP", limit){
  grid_plot(data = data, label = NULL, value = value)+
    scale_fill_gradientn(name=expression("Population-weighted PM"[2.5]*"("*"μg/m"^3* ")"),
                         colours = cols(10), limits = limit, na.value="white",
                         guide = guide_colorbar(frame.colour = "black",frame.linewidth=2,
                                                title.position="right",title.hjust=0.5,title.vjust=1,
                                                byrow = TRUE, ticks= F))+ ##set colour
    theme(panel.grid.major = element_line(colour = alpha("black", 0.4), size = 0.6))+
    geom_sf(data = polygon_sf, fill = "white", alpha = 0.1, size = 0.4)+
    geom_sf(data = Nline_sf, fill = "white", size = 0.6)+
    xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(strip.background = element_rect(color="black", fill="white", size = 0.5))+ ##change the facetwrap background
    theme(axis.text=element_text(size = 30))+
    theme(legend.key.width=unit(0.8, "cm"))+
    theme(legend.key.height=unit(3, "cm"))+##legend height
    theme(legend.title = element_text(colour = "black", size = 30))
}

IPM <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5.csv")
IPM%>%
  group_by(year)%>%
  summarise(AAP_ind = weighted.mean((AAP+HAP_O),Pop),
            HAP_ind = weighted.mean(HAP_I,Pop))%>%
  as.data.frame()

cols <-  colorRampPalette(c( "#009966", "#ffde33", "#ff9933",
                             "#cc0033", "#660099", "#7e0023"))
year = 2000
exp_table = tibble(exp = c("AAP_ind","HAP_ind"), limit = list(c(0.01,120),c(0.01,200)))

for (year in 2000:2020) {
  cat("======Year: ",year,"======\n")
  yl <- year

  IPM_spatial <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5.csv")%>%
    filter(year == yl)%>%
    mutate(AAP_ind = AAP+HAP_O,HAP_ind = HAP_I)%>%
    dplyr::select(x,y,IAP,AAP_ind,HAP_ind)%>%
    dfproj(method = "near")
    # mutate(IAP = cut(IAP, breaks, labels))
  
  1:2%>%
    purrr::map(~ {
      IPM_spatial_plot <- IPM_plot(IPM_spatial,exp_table$exp[.x],exp_table$limit[[.x]])
      
      Cairo::CairoPNG(paste0("result/Fig/PM2.5_Exposure/",exp_table$exp[.x],"_Exposure_",year,".png"),
                      width = 2800,height = 1500, res=150)
      print(IPM_spatial_plot)
      dev.off()
    })
}

##=======Integrated exposure change during 2000-2020========
IPM_national <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5.csv")%>%
  group_by(year)%>%
  dplyr::summarise(IAP = weighted.mean(IAP, Pop),
            AAP = weighted.mean(AAP, Pop),
            HAP_O = weighted.mean(HAP_O, Pop),
            HAP_I = weighted.mean(HAP_I, Pop),
            # check = weighted.mean(HAP_O+AAP,Pop),
            Popnf = sum(Popnf), Pop = sum(Pop))%>%ungroup()%>%
  mutate(HAP_ind = HAP_I, AAP_ind = HAP_O + AAP)%>%
  dplyr::select(year,HAP_ind,AAP_ind)%>%
  pivot_longer(HAP_ind:AAP_ind, names_to = "Exposure")%>%
  mutate(Exposure = ifelse(Exposure == "HAP_ind", "Indoor source","Outdoor source"))

IPM_national_bygroup <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5_by_group.csv")%>%
  rename(IAP = concentration)%>%
  group_by(year,Area,Fuel)%>%
  summarise(IAP = weighted.mean(IAP, Pop),
            AAP = weighted.mean(AAP, Pop),
            HAP_O = weighted.mean(HAP_O, Pop),
            HAP_I = weighted.mean(HAP_I, Pop),
            # check = weighted.mean(HAP_O+AAP,Pop),
            Popnf = sum(Popnf), Pop = sum(Pop))%>%ungroup()%>%
  mutate(HAP_ind = HAP_I, AAP_ind = HAP_O + AAP)%>%
  dplyr::select(year,Area,Fuel,HAP_ind,AAP_ind,Pop,Popnf)%>%
  pivot_longer(HAP_ind:AAP_ind, names_to = "Exposure")%>%
  mutate(Exposure = ifelse(Exposure == "HAP_ind", "Indoor source","Outdoor source"))
fwrite(IPM_national_bygroup, "result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5_national.csv")

IPM_time_plot <- function(data){
  ggplot(data,aes(x = year, y = value))+
    geom_bar(aes(fill = Exposure), stat = "identity", 
             position = "stack", color = "white")+
    labs(x= NULL, y = expression(paste("Population-weighted exposure (µg/m"^"3",")")))+
    scale_fill_manual(name = NULL, values = c("#69AFDF","#E0BCC1"))+
    theme(panel.background = element_rect(size=1,fill='white', colour='black'),
          panel.border = element_rect(size = 1,colour = "black", fill=NA),
          legend.background = element_blank(),
          legend.box.background = element_rect(size= 0.5,colour = "black"))+
    theme(axis.text=element_text(colour = "black",size = rel(2.5)), ##axis text size
          axis.title = element_text(colour = "black",size = rel(3)),
          plot.margin=unit(c(5,0,5,0),"mm"))+ ##axis title size
    theme(legend.key.width=unit(1, "cm"),
          legend.key.height=unit(1, "cm"),  ##legend size
          legend.position = "top",
          legend.text = element_text(colour = "black",size = rel(2.5))) ##legend text size
  
}

CairoPNG(paste0("result\\Fig\\PM2.5_Exposure\\Exposure_change.png"),height = 1700,width = 2800, res = 200)
print(IPM_time_visu)
dev.off()


IPM_national_groupplot <- IPM_time_plot(IPM_national_bygroup)+
  facet_grid(Fuel ~ Area,labeller = labeller(groupwrap = label_wrap_gen(10)))+
  theme(strip.background = element_rect(color="black", size = 1))+ ##change the facetwrap background
  theme(strip.text = element_text(size = 32))
CairoPNG(paste0("result\\Fig\\PM2.5_Exposure\\Exposure_change_by_group.png"),height = 4500,width = 5000, res = 200)
print(IPM_national_groupplot)
dev.off()


##=======Cumulative percentage of population exposure during 2000-2020========
library(RColorBrewer)
IPM_by_pop <- fread("result/Table/Integrated_PM2.5_exposure/PM2.5_group_by_population.csv")%>%
  filter(year%in%c(2000,2020), concentration < 300)%>%
  mutate(year = factor(year))%>%group_by(year)%>%
  mutate(Popf_IAP = cumsum(Popf_IAP), Popf_AAP = cumsum(Popf_AAP), Popf_HAP = cumsum(Popf_HAP))%>%
  ungroup()%>%pivot_longer(cols = c("Popf_IAP","Popf_AAP", "Popf_HAP"))%>%
  mutate(label = paste0(year,", ", ifelse(name == "Popf_IAP","IAP",
         ifelse(name == "Popf_AAP","AAP","HAP")))%>%
         factor(levels = c(paste0(c(2000,2020),", IAP"), 
                           paste0(c(2000,2020),", AAP"),
                           paste0(c(2000,2020),", HAP"))))

IPM_by_popf <- IPM_by_pop%>%filter(name%in%c("Popf_IAP","Popf_AAP"))

color = c("#CE161F","#EC684F","#0A4990","#6DADD3")
bgcol = RColorBrewer::brewer.pal(6,"Blues")
rects = data.frame(xstart = c(0,10,15,25,35), xend = c(10,15,25,35,300),col = letters[1:5])

PM_by_popplot = ggplot()+
  geom_rect(data = rects, aes(xmin = xstart,
            xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.2)+
  geom_line(data = IPM_by_popf, aes(x = concentration, y = value, color = label), size = 1.2)+
  scale_y_continuous("Cumulative population fraction (%)",
                     expand =  expansion(mult = 0))+
  scale_x_continuous(expression(paste("PM"[2.5]," exposure (µg/m"^"3",")")),
                     expand = expansion(mult = c(0,0.04)))+
  theme_bw()+ 
  theme(legend.key.width = unit(1, "cm"))+
  theme(axis.text = element_text(colour = "black", size = rel(2.5)))+
  theme(axis.title.x = element_text(colour = "black", size = rel(3),
                                    margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(colour = "black", size = rel(3),
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  theme(legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(colour = "black", size = rel(2.8)), ##legend text size
        legend.title = element_blank(), legend.spacing.y = unit(1, "cm"),
        legend.position = c(0.8, 0.3),
        legend.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))+
  guides(color = guide_legend(byrow = TRUE))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,0,0,0), "cm"), strip.background = element_blank())+
  scale_color_manual(values = color)+
  scale_fill_manual(values = bgcol, guide = "none")

CairoPNG(paste0("result/Fig/PM2.5_Exposure/Cumulative_Concentration.png"),
         width = 1800,height = 1500,res = 180)
PM_by_popplot
dev.off()

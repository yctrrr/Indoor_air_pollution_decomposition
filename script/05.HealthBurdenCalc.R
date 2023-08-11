# Core functions load
library(Cairo)

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source('./script/Core.R', encoding = 'UTF8')
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')

#---- Calculate PM2.5 mortality ----
# C-R function setting 
# this section used to choose C-R function for Health Impact Calculation.
# Supported C-R functions:  
#              'IER', 'NCD+LRI'(Part of GEMM), '5COD'(Part of GEMM), 'MRBRT'

use_CR('MRBRT')

# Data load
read_files(
  GRID = './result/Table/Integrated_PM2.5_exposure/GRID_information.csv',
  Pop = paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'),
  PM_real = paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_2000.csv"),
  PM_cf = './data/PM_Ctrl.csv', # PM_cf works only in counter-fact scenario
  MortRate = './data/GBD_incidence_China_2000-2019.csv',
  AgeGroup = './data/GBD_agestructure_China_2000-2017.csv'
)


AgeGroupt = fread('./data/IHME_GBD_2019_POP_SYA_2000-2020.csv')%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4,include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(sex = tolower(sex), agegroup = as.character(agegroup))

MortRatet = fread('./data/IHME_GBD_2019_Deathrate_mod_2000-2020.csv')%>%
  filter(metric == "Rate")%>%
  dplyr::select(endpoint, sex, agegroup, MortRate, MortRateU, MortRateL, year)%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4, include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(endpoint = tolower(endpoint),sex = tolower(sex), agegroup = as.character(agegroup))
# RR_table$MEAN$PTB_ALL
RR = (RR_table$MEAN)%>%RR_std()%>%
  mutate(agediv = cut(as.numeric(agegroup), breaks = c(0,5,15,65,96),
                      labels = 1:4,
                      include.lowest = TRUE, right = FALSE)%>%as.numeric())

for (y in as.character(2000:2020)) {
  yl = y

  tic(paste0("=====Calculate PM2.5 Mortality in: ",y,"=====\n"))

  Grids = Grid_info
  Popy = Pop %>% filter(year == as.numeric(yl))%>%
    mutate(Pop = Pop*Popf,.keep = "unused")%>%dplyr::select(-year)

  ##  filter corresponding age and sex groups
  ag = AgeGroupt%>%filter(year == as.numeric(yl))%>%
    dplyr::select(sex, agegroup, AgeStruc, agediv)

  mRate = MortRatet%>%filter(year == as.numeric(yl))%>%
    dplyr::select(endpoint, sex, agegroup, MortRate, agediv)

  expo_type <- c("AAP","HAP_O","HAP_I","AAP_ind","IAP")
  PM25_Morty <- 1:length(expo_type)%>%
    purrr::map_dfr(
      function(n){
        PM_r <- fread(paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_",yl,".csv"))%>%
          rename(IAP = concentration)%>%
          mutate(AAP_ind = AAP + HAP_O, HAP_I = ifelse(HAP_I < 0.1, 0, HAP_I))%>%
          dplyr::select(x,y,sex,agediv,Area,Fuel,concentration := !! expo_type[n])%>%
          mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)%>%
          mutate(concentration = matchable(concentration))

        ##  calculate mortality at grid level
        M = Mortality(Grids = Grids, RR = RR, PM_r = PM_r, pop = Popy, 
                      ag = ag, mRate = mRate)%>%
          mutate(expo_type = expo_type[n])
        return(M)
      }
    )

  fwrite(PM25_Morty, paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1_',y,'.csv'))
  toc()
} ## end of year

PM_r <- fread(paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_",y,".csv"))%>%
  rename(IAP = concentration)%>%
  mutate(AAP_ind = AAP + HAP_O, HAP_I = ifelse(HAP_I < 0.1, 0, HAP_I))



## calculation check
IPM_pop <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5_by_group.csv")%>%
  mutate(AAP_indc = AAP+HAP_O,HAP_indc = HAP_I)%>%
  dplyr::select(x:Fuel,year,Popnf,Pop,AAP_indc,HAP_indc)%>%
  filter(year == yl)
IPM_pop

## output tables
2000:2020%>%
  map_dfr(~ fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1_',.x,'.csv'))%>%
            filter(expo_type == "IAP"|expo_type == "AAP_ind"|expo_type == "HAP_I")%>%
            group_by(x,y,expo_type)%>%
            summarise(Mort = sum(Mort))%>%ungroup()%>%
            pivot_wider(x:y, names_from = expo_type, values_from = Mort)%>%
            mutate(AAP_ind = ifelse(AAP_ind + HAP_I == 0,0,IAP*AAP_ind/(AAP_ind + HAP_I)),
                   HAP_I = IAP - AAP_ind)%>%
            pivot_longer(cols = AAP_ind:IAP, names_to = "expo_type", values_to = "Mort")%>%
            mutate(year = .x))%>%
  fwrite(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1.csv'))

fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1.csv'))%>%
  filter(year %in% c(2002,2007,2012,2017))%>%
  pivot_wider(names_from = expo_type, values_from = Mort)%>%
  group_by(year)%>%
  dplyr::summarise_at(vars(AAP_ind, HAP_I, IAP),sum)%>%ungroup()%>%
  mutate(AAP_ratio = AAP_ind/IAP)

##==== Visualization of PM2.5 premature death =====
## Mortality population plot function
Mort_plot <- function(data ,value = "Mort"){
  grid_plot(data = data, label = NULL, value = value)+
    scale_fill_manual(name = expression(paste("PM"[2.5],"-related death (#/ 100km"^"2",")")),
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.vjust = 0.5,
                                           byrow = TRUE, ticks = F))+
    geom_sf(data = polygon_sf, fill = "white", alpha = 0.1, linewidth = 0.4)+
    geom_sf(data = Nline_sf, fill = "white", linewidth = 0.6)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(axis.text=element_text(size = 35),
          plot.margin=unit(c(5,0,5,0),"mm"))+
    theme(panel.grid.major = element_line(colour = alpha("black", 0.4), size = 0.4))+
    theme(legend.title = element_text(colour = "black", size = 40),
          legend.text=element_text(size=30),
          legend.key.width=unit(1, "cm"),
          legend.key.height=unit(1, "cm"),##legend height
          legend.spacing.y = unit(.8, 'cm'),
          legend.box.just = "center",
          legend.title.align = 0.5)
}
library(Cairo)
polygon <- readOGR(dsn="data/Shp/Province.shp",
                   layer="Province",use_iconv=TRUE, encoding = "UTF-8")

Nline <- readOGR("data/Shp/Nline.shp",
                 layer="Nline",use_iconv=TRUE, encoding = "UTF-8")
polygon_sf <- st_as_sf(polygon)
Nline_sf <- st_as_sf(Nline)

breaks <- c(-Inf, 0.1, 1, 5, 10, 25, 50, 75, 100, 150, +Inf)
labels <- c("0 - 0.1", "0.1 - 1", "1 - 5","5 - 10", "10 - 25", "25 - 50", 
            "50 - 75", "75 - 100", "100 - 150", "> 150")
colors <- rev(RColorBrewer::brewer.pal(length(labels),"Spectral"))

for (expo in c("IAP","AAP_ind","HAP_I")) {
  grid_aggr <- c(2000,2020)%>%
    map_dfr(function(t){
      fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1.csv'))%>%
        filter(year == t, expo_type == expo)%>%
        dplyr::select(-expo_type)%>%dfproj()%>%
        mutate(Mort = cut(Mort, breaks, labels),
               year = t)
    })
  
  grid_aggr_plot <- Mort_plot(grid_aggr)+
    facet_wrap("year", nrow = 1, labeller = labeller(groupwrap = label_wrap_gen(20)))+
    theme(strip.background = element_blank())+ ##change the facetwrap background
    theme(strip.text = element_text(size = 40),
          axis.text = element_text(size = 30),
          legend.box.margin=margin(c(120,0,0,0)))
  
  grid_aggr_plot <- Nanhai_zoom(grid_aggr_plot)
  
  CairoPNG(paste0("result/Fig/PM2.5_Death/PM2.5_Death_Decomposed_Year_",expo,".png"),width = 3000,height = 1500,res=150)
  print(grid_aggr_plot)
  dev.off()
}


##========Change in Death by time series======
Mort_national <- fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1.csv'))%>%
  filter(expo_type%in%c("AAP_ind","HAP_I"))%>%
  group_by(year,expo_type)%>%
  summarise(Mort = sum(Mort))%>%ungroup()%>%
  mutate(expo_type = ifelse(expo_type == "AAP_ind","Outdoor source","Indoor source"),
         Mort = Mort/1e3)

Mort_time_visu <- ggplot(Mort_national, aes(x = year, y = Mort))+
  geom_bar(aes(fill = expo_type), stat = "identity", 
           position = "stack", color = "white")+
  labs(x= NULL, y = expression(paste("PM"[2.5],"-related death (thousand)")))+
  scale_fill_manual(name = NULL, values = c("#b3e2cd","#f4cae4"))+
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
CairoPNG(paste0("result\\Fig\\PM2.5_Death\\Death_change.png"),height = 1700,width = 2800, res = 200)
print(Mort_time_visu)
dev.off()

##========Changes in Death at grid scael========
PM25_Mortsel <- fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1.csv'))%>%
  filter(year == 2000|year == 2020)%>%
  pivot_wider(x:expo_type, names_from = year, values_from = Mort)%>%
  mutate(Mortc = `2020` - `2000`,.keep = "unused")%>%
  pivot_wider(x:y, names_from = expo_type, values_from = Mortc)

PM25_Mortsel%>%
  summarise_at(vars(AAP_ind,HAP_I,IAP),sum)
PM25_Mortsel

breaks <- c(-Inf, -100, -50,-10, -5, -1, 1, 5, 10, +Inf)
labels <- c("< -100", "-100 - -50","-50 - -10", "-10 - -5", "-5 - -1","-1 - 1",
            "1 - 5","5 - 10", "> 10")
colors <- colorspace::sequential_hcl(length(labels)-4, palette = "agGrnYl")%>%
  c("white","#FEE08B","#FDAE61","#F46D43")

PM25_Mortsel_mod <- PM25_Mortsel%>%dfproj()%>%
  pivot_longer(AAP_ind:IAP,values_to = "Mortc")%>%
  mutate(Mortc = cut(Mortc,breaks,labels),
         name = ifelse(name == "AAP_ind","Outdoor source",
                ifelse(name == "HAP_I","Indoor source","Total")))

PM25_Mortsel_plot <- Mort_plot(PM25_Mortsel_mod, value = "Mortc")+
  facet_wrap("name", nrow = 1, labeller = labeller(groupwrap = label_wrap_gen(20)))+
  scale_fill_manual(name = expression(paste("Change in deaths (#/ 100km"^"2",")")),
                    values = colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                         title.position="right",title.vjust = 0.5,
                                         byrow = TRUE, ticks = F))+
  theme(strip.background = element_blank())+ ##change the facetwrap background
  theme(strip.text = element_text(size = 40),
        axis.text = element_text(size = 30),
        legend.box.margin=margin(c(120,0,0,0)))
PM25_Mortsel_plot <- Nanhai_zoom(PM25_Mortsel_plot)
CairoPNG(paste0("result/Fig/PM2.5_Death/Integrated_PM2.5_Death_Change_2000-2020.png"),
         width = 4000,height = 1500,res = 150)
print(PM25_Mortsel_plot)
dev.off()

##========Changes in Death by province(Table)========
prov_info <- fread("result/Table/Input_data/province_information.csv", encoding = 'UTF-8')
grid_info <- fread('result/Table/Integrated_PM2.5_exposure/GRID_information.csv')

PM25_Mort_prov <- fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1.csv'))%>%
  filter(year == 2000|year == 2020)%>%
  pivot_wider(x:expo_type, names_from = year, values_from = Mort)%>%
  mutate(Mortc = `2020` - `2000`)%>%
  left_join(grid_info)%>%left_join(prov_info)%>%
  group_by(Province, expo_type)%>%
  summarise(Mort_2000 = sum(`2000`)%>%round(),Mort_2020 = sum(`2020`)%>%round(),
            Mortc = sum(Mortc))%>%ungroup()%>%
  mutate(Mortc = (Mortc/Mort_2000*100)%>%formatC(format = "f",
                                                 digits = 1))%>%
  pivot_wider(Province, names_from = expo_type, 
              values_from = c(Mort_2000, Mort_2020, Mortc),
              names_vary = "slowest")

fwrite(PM25_Mort_prov, 'result/Table/Integrated_PM2.5_death/Mortality_group_by_province.csv')

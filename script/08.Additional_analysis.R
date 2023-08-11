library(Cairo)
library(terra)
library(tictoc)
library(data.table)
library(dplyr)
library(foreach)
library(furrr)

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source('./script/Core.R', encoding = 'UTF8')
source('./script/Core_MonteCarlo.R',encoding = 'UTF8')
source('./script/Decomposition.R', encoding = 'UTF8')
GRID = './result/Table/Integrated_PM2.5_exposure/GRID_information.csv'
use_CR('MRBRT')
read_files(
  GRID = './result/Table/Integrated_PM2.5_exposure/GRID_information.csv',
  Pop = paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'),
  PM_real = paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_2000.csv"),
  PM_cf = './data/PM_Ctrl.csv', # PM_cf works only in counter-fact scenario
  MortRate = './data/GBD_incidence_China_2000-2019.csv',
  AgeGroup = './data/GBD_agestructure_China_2000-2017.csv'
)

AgeGroupt = fread('./data/IHME_GBD_2019_POP_SYA_2000-2020.csv')%>%
  # dplyr::select(sex, agegroup, val, upper,lower,year)%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4, include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(sex = tolower(sex), agegroup = as.character(agegroup))

Pop_long_sum = Pop%>%group_by(year, Area, Fuel)%>%
  summarise(Pop = sum(Pop*Popf))%>%ungroup()%>%
  group_by(year)%>%mutate(Pop_all = sum(Pop))%>%ungroup()%>%
  mutate(Popf = Pop/Pop_all)

##========Calculate population specific mortality in 2010=======
Mortality_by_pop <- 2000:2020%>%
  purrr::map_dfr(~ fread(paste0('./result/Table/Integrated_PM2.5_death/Integrated_PM2.5_Death_g0.1_',.x,'.csv'))%>% 
                   split(., .[, c("Area","Fuel")])%>%
                   purrr::map(function(k) dplyr::select(k, -Area, -Fuel) %>% 
                                mutate(x = round(x,2) %>% str_c, y = round(y,2) %>% str_c)%>%
                                summarise(Mort = sum(Mort, na.rm = TRUE)))%>%
                   unlist()%>%data.frame(Poptype = gsub(".Mort","",names(.)), Mort = .)%>%
                   mutate(year = .x))%>%
  mutate(Area = str_split_fixed(Poptype, "[.]", n = 2)[,1],
         Fuel = str_split_fixed(Poptype, "[.]", n = 2)[,2], Mort = round(Mort))

fwrite(Mortality_by_pop, "result/Table/Integrated_PM2.5_death/Mortality_group_by_population.csv")

##=======Population inequality========
Pop_sum <- fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))%>%
  group_by(year, Area, Fuel)%>%
  summarise(Pop = sum(Pop))%>%ungroup()%>%
  group_by(year)%>%mutate(Pop_all = sum(Pop))%>%ungroup()%>%
  mutate(Popf = Pop/Pop_all)

Mortality_inequ = fread("result/Table/Integrated_PM2.5_death/Mortality_group_by_population.csv")%>%
  group_by(year)%>%mutate(Mort_all = sum(Mort))%>%ungroup()%>%
  mutate(Mortf = Mort/Mort_all)%>%left_join(Pop_sum)%>%
  group_by(year, Area)%>% summarise(Mortf = sum(Mortf), Popf = sum(Popf))%>%ungroup()%>%
  mutate(NGD = (Mortf - Popf)/Popf) ## calculation of disparity

Mortality_inequ_plot <- ggplot(Mortality_inequ, aes(year, NGD, fill = Area))+
  geom_bar(stat = "identity", position ="identity", width = 1, color = "black")+
  theme_bw()+ 
  guides(fill = guide_legend(byrow = TRUE))+
  theme(panel.border = element_blank(),
        axis.text = element_text(colour = "black", size = 35),
        axis.title = element_text(size = 40),
        legend.text = element_text(size = 35), 
        legend.title = element_text(colour = "black", size = 40),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  scale_fill_manual(values = c("#EE6D59","#69AFDF"))+
  labs(x = NULL, y = "Normalized Health Disparity")

Mortality_inequ_plot

CairoPNG(paste0("result/Fig/PM2.5_Death/Integrated_PM2.5_Death_by_time_series.png"),
         width = 2500,height = 1500, res = 120)
print(Mortality_inequ_plot)
dev.off()


##=======Mortality fraction========
Mortality_by_pop <- fread("result/Table/Integrated_PM2.5_death/Mortality_group_by_population.csv")%>%
  mutate(Poptype = gsub("\\."," ",Poptype))

Mortality_by_popf = filter(Mortality_by_pop,
                           year%in%c(2000,2007,2013,2020))%>%
  group_by(year)%>%
  mutate(Mort = (Mort/sum(Mort)*100), pos = cumsum(Mort)-0.5*Mort)%>%ungroup()%>%
  mutate(text = formatC(Mort, format = "f", digits = 1),
         text = ifelse(Mort > 10,
                       paste0(text,"%"),"")
         )

Mortality_by_popplot <- ggplot(Mortality_by_popf, aes(x="", y = Mort, fill = Poptype)) +
  geom_bar(stat = "identity", color = "black", size = 0.3, width = 0.5) +
  coord_polar("y", start = 0, direction = -1)+
  geom_text(aes(x="", y = Mort, label = text),
            position = position_stack(vjust = 0.5),size = 6.5)+
  facet_wrap(vars(year), nrow = 2)+
  theme_void()+labs(fill = "Population type")+
  theme(legend.text = element_text(colour = "black", size = 30))+ ##legend text size
  theme(legend.title = element_text(colour = "black", size = 30))+
  theme(legend.key.width=unit(1, "cm"),
        legend.key.height=unit(1, "cm"),
        legend.spacing.y = unit(.8, 'cm'))+
  theme(strip.text = element_text(size = 40))+
  theme(panel.spacing.x = unit(2, "cm"),
        legend.box.margin=margin(c(50,50,50,50)))+
  guides(fill = guide_legend(byrow = TRUE))+
  scale_fill_brewer(palette="Pastel1")

CairoPNG(paste0("result/Fig/PM2.5_Death/Integrated_PM2.5_Death_by_population_group.png"),
         width = 2400,height = 1500, res=140)
print(Mortality_by_popplot)
dev.off()


##=======**Visualization of age structure, population=======
AgeGroup_toplot <- AgeGroupt%>%
  group_by(year, agegroup)%>%
  summarise(AgeStruc = sum(AgeStruc))%>%
  ungroup()%>%
  mutate(agegroup = ifelse(agegroup < 15, "0 - 14",ifelse(agegroup<65,"15 - 64","65+")))

AgeGroup_plot <- ggplot(AgeGroup_toplot, aes(year, AgeStruc, fill = agegroup))+
  geom_bar(stat = "identity", width = 0.3, position = position_stack())+
  theme_bw()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "white"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  labs(x = "Year", y = "Population Fraction")
ggsave(paste0("result/Fig/PM2.5_Death/Age_structure.png"),AgeGroup_plot)


MortRate_toplot <- MortRatet%>%
  group_by(year, endpoint)%>%
  summarise(MortRate = sum(MortRate))%>%
  ungroup()

Pop_toplot <- Pop%>%
  group_by(year, Area)%>%
  summarise(Pop = sum(Pop*Popf)/1e9)%>%
  ungroup()


Pop_plot <- ggplot(Pop_toplot)+
  geom_smooth(aes(x = year, y = Pop, group = Area, color = Area, fill = Area), size = 1.5)+
  geom_point(aes(x = year, y = Pop, group = Area, color = Area), size = 1.5)+
  theme_bw()+
  scale_y_continuous(limits = c(0.2,1.2), breaks = seq(0.2,1.2,0.2))+
  scale_color_manual(name = "",values = c("#FFC000","#4EB272"))+
  scale_fill_manual(name = "",values = c("#FFC000","#4EB272"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth = 0.5, fill = NA, colour = 'black'),
        panel.border = element_rect(linewidth = 1, colour = "black", fill = NA)
        )+
  labs(x = "", y = "Population (billion)")+
  theme(axis.text = element_text(size = 25),
        axis.title.y = element_text(size = 28,
                                  margin = margin(0,20,0,0)),
        legend.text=element_text(size = 25),
        legend.position = c(.8,.9),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.key.width=unit(2.5, "cm"),
        legend.key.height=unit(1, "cm"),
        legend.spacing.y = unit(.8, 'cm'),
        legend.spacing.x = unit(.5, 'cm'))+
  guides(color = guide_legend(byrow = TRUE), title.hjust = 1)

CairoPNG(paste0("result/Fig/Fuel_use/Population_change.png"), 
      width = 2200,height = 1200,res = 140)
print(Pop_plot)
dev.off()
ggsave("result/Fig/Fuel_use/Population_change.png",Pop_plot,
       width = 12, height = 7)

##=======**Visualization of PM2.5 RR curve=======
RR_mean = (RR_table$MEAN)%>%RR_std()%>%rename(RR_mean = RR)
RR_low = (RR_table$LOW)%>%RR_std()%>%rename(RR_low = RR)
RR_up = (RR_table$UP)%>%RR_std()%>%rename(RR_up = RR)

for (age in c(25,50,75,999)) {
  # RR_mean%>%filter(concentration == 200, agegroup == 75)
  if(age == 999){
    RR_curve <- list(RR_mean,RR_low,RR_up)%>%reduce(left_join)%>%
      mutate(concentration = as.numeric(concentration),endpoint = toupper(endpoint))%>%
      filter(concentration <= 300)%>%
      group_by(concentration,endpoint)%>%
      summarise(RR_mean = mean(RR_mean), RR_low = mean(RR_low),
                RR_up = mean(RR_up))%>%ungroup()
  }else{
    RR_curve <- list(RR_mean,RR_low,RR_up)%>%reduce(left_join)%>%
      mutate(concentration = as.numeric(concentration))%>%
      filter(agegroup == age, concentration <= 300)%>%
      mutate(endpoint = toupper(endpoint))
  }
  
  color = c('#FD6D5A', '#FEB40B', '#6DC354', '#994487', '#518CD8', '#443295')
# RR_curve%>%filter(concentration ==)
  RR_plot <- ggplot(RR_curve, aes(x = concentration, y = RR_mean))+
    geom_line(aes(group = endpoint, color = endpoint),linewidth = .6)+
    facet_wrap(~ endpoint)+
    scale_fill_manual(name = "Endpoint",values = color)+
    scale_color_manual(name = "Endpoint",values = color)+
    # geom_point(size = 0.1)+
    geom_ribbon(aes(ymax = RR_up, ymin = RR_low, fill = endpoint), alpha = 0.5)+
    theme_bw()+
    guides(fill = guide_legend(byrow = TRUE))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank())+
    labs(x = expression(paste("Concentration (µg/m"^"3",")")), y = "Relative risk")

  ggsave(paste0("result/Fig/Supp/Relative_risk_age_",age,".png"), RR_plot)
}

##=======**Visualization of economic growth and temperature change=======
prov_info <- fread("result/Table/Input_data/province_information.csv", encoding = 'UTF-8')%>%
  mutate(Province = whoname)
grid_info <- fread('./result/Table/Integrated_PM2.5_exposure/GRID_information.csv')

prov_cov <- readRDS(paste0("data/panel_data/Chinese_provincial_panel_data.Rds"))%>%
  dplyr::select(Province, Year, all_of(IV))%>%
  left_join(prov_info,.)%>%
  mutate(loc_id = pindex)
population_prov <- fread("result\\Table\\Input_data\\province_population.csv")%>%
  dplyr::select(-Province,-Subarea)

gdp_dev <- read_rds(paste0("D:/shaoyanchuan/data/Economic/CSMAR/CSMAR_province_panel.Rds"))%>%
  filter(Province == "中国",Year >= 2000)%>%
  rename(pcGDP = Gdp0116, pcUE = Urban_expenditure_pc,pcRE = Rural_expenditure_pc)%>%
  dplyr::select(Year,pcGDP,pcUE,pcRE)

temp_inf <- fread("result\\Table\\Covariates\\Temperature_g0.1.csv")%>%
  left_join(grid_info)%>%na.omit()%>%
  group_by(year)%>%
  summarise(hdd = mean(hd), cdd = mean(cd), t2m = mean(t2m))%>%
  ungroup()%>%rename(Year = year)
annual_change <- left_join(gdp_dev, temp_inf)

trend_plot <- function(data, value = "value", group = "name"){
  ggplot(data)+
    geom_point(aes_string(x = "Year", y = value, group = group, color = group), size = 1.5)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(linewidth = 0.5, fill = NA, colour = 'black'),
          panel.border = element_rect(linewidth = 1, colour = "black", fill = NA)
    )+
    theme(axis.text=element_text(size = 25),
          axis.title.y=element_text(size = 28, margin = margin(0,20,0,0)),
          legend.text=element_text(size = 25),
          legend.position = "top",
          legend.background = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.key.width=unit(2.5, "cm"),
          legend.key.height=unit(1, "cm"),
          legend.spacing.y = unit(.8, 'cm'),
          legend.spacing.x = unit(.5, 'cm'))+
    guides(color = guide_legend(byrow = TRUE), title.hjust = 1)
  
}
GDP_plot <- annual_change%>%dplyr::select(Year,pcGDP,pcUE,pcRE)%>%
  pivot_longer(pcGDP:pcRE)%>%
  mutate(value = value/1e3,name = factor(name, levels = c("pcGDP","pcUE","pcRE")))%>%
  trend_plot()+
  geom_smooth(aes_string(x = "Year", y = "value", group = "name", 
                         color = "name", fill = "name"), size = 1.5)+
  scale_color_manual(name = "",values = c("#B6122A","#FAD336","#4EB272"))+
  scale_fill_manual(name = "",values = c("#B6122A","#FAD336","#4EB272"))+
  labs(x = "", y = "Thousand Chinese yuan/person")
ggsave("result/Fig/Supp/GDP_trend_plot.png",GDP_plot,
       width = 12, height = 7)


##=======Abstract graph=======
## Fuel use plot
Areavec = c("Urban","Rural")

Shape_eco <- fread("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")%>%
  mutate(Mort = replace_na(Mort,0))%>%
  filter(Decomp %in% c("FF","ECO"))%>%
  group_by(year, Decomp, sample_id)%>%
  dplyr::summarise(value = sum(Mort,na.rm = TRUE))%>%ungroup()%>%
  Decomp_to_plot()%>%
  filter(!Decomp%in%"dock")%>%
  mutate(lower = round(lower/1000), upper = round(upper/1000))%>%
  mutate(year = as.character(year))

Shape_eco_group <- 1:2%>%
  purrr::map_dfr(~ fread("result\\Table\\Uncertainty\\Mortality_Decompisition_uncertainty.csv")%>%
                   group_by(year, Area, Decomp, sample_id)%>%
                   dplyr::summarise(value = sum(Mort, na.rm = TRUE))%>%ungroup()%>%
                   filter(Area == Areavec[.x])%>%
                   Decomp_to_plot()%>%
                   mutate(name = .x, Area = Areavec[.x])%>%
                   filter(Decomp%in%c("ECO"))%>%
                   mutate(mean = -mean, 
                          dock_upper = -round(lower/1000),
                          dock_lower = -round(upper/1000),
                          year = as.character(year))
                 )

deco_num = 0
deco_yname = c("2000 - 2007","2007 - 2013","2013 - 2020")

color_palette = c('2000' = '#ffeda0', '2007' = '#feb24c',
                  '2013' = '#f03b20')


Shape_facet_plot <- function(data, limits = NULL){
  ggplot(data)+
    geom_bar(aes(x = seq, y = mean, fill = year),
             color = NA, linewidth = 1.5,
             stat = 'identity', position ="identity", width = .8)+
    geom_errorbar(mapping = aes(x = seq, ymin = dock_lower, ymax = dock_upper),
                  linewidth = 2, width = .3, color = "black")+
    labs(y = NULL)+
    scale_fill_manual(values = color_palette, name="")+
    scale_x_discrete(labels = deco_yname)+
    scale_y_continuous(
      name = "", limits = limits,
      expand = c(0.01, 0))+
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA, colour = "black", linewidth= 2.5),
      plot.margin = unit(c(.2, .5, .1, .2), "cm"),
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(color = 'black', size = 50),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 50,hjust = 0.5)
      )
}

## save graphs
limitsvec = list(c(0,70),c(0,220))
1:2%>%
  map(.f = function(x){
    Shape_eco_groupplot <- Shape_eco_group%>%
      filter(Area == Areavec[x])%>%
      Shape_facet_plot(limits = limitsvec[[x]])
    
    CairoPNG(paste0("result/Fig/Decomposition/Decomposition_of_Economic_Growth_",
                    Areavec[x],".png"),
             width = 3000,height = 1500, res = 180)
    print(Shape_eco_groupplot)
    dev.off()
  })


##=======Visualization of population=======
breaks <- c(-Inf, 0.1, 10, 50, 100, 500, 1000, 5000, 10000, 20000, +Inf)
labels <- c("0", "0 - 10", "10 - 50", "50 - 100","100 - 500","500 - 1000", 
            "1000 - 5000", "5000 - 10000", "10000 - 20000", "> 20000")
colors <- rev(colorspace::sequential_hcl(length(labels), palette = "Heat"))

fuel_pop_aggr <- fread("result\\Table\\Population\\Fuel_use_population_g0.1.csv")%>%
  filter(year == yl)%>%mutate(Pop = Pop*Popf)%>%
  group_by(x,y)%>%summarise(Pop = sum(Pop))%>%ungroup()%>%
  dfproj()%>%
  mutate(Pop = cut(Pop, breaks, labels))

grid_plot(data = fuel_pop_aggr, label = NULL, value = "Pop")+
  scale_fill_manual(name = "Population",
                    values = colors,
                    guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                         title.position="right",title.hjust=0.5,title.vjust=1,
                                         byrow = TRUE, ticks= F))+
  geom_sf(data = polygon_sf, fill = "white", alpha = 0.1, size = 0.4)+
  geom_sf(data = Nline_sf, fill = "white", size = 0.6)+
  coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)
  
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(dplyr)
library(pinyin)
library(reshape2)
library(stringr)
library(ggplot2)
library(Cairo)
dat_dir <- "D:\\shaoyanchuan\\data\\"
source("D:/shaoyanchuan/codebook/function/Scatter_plot.R")
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Data_Setup.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Functions_parallel2.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Model_parallel_provincial.R"))


##========Prepare Fuel use related data========
##  Read in fuel use survey data
survey_data <- readRDS(paste0("data/population/Fuel_type_survey_province_long_format.Rds"))
fuel_indx <- c("Clean","Biomass","Coal")
survey_data_lf <- survey_data%>%
  rename(Area = Population_type, year = Year, Wood = Firewood, Crop = Other_energy)%>%
  mutate(Area = ifelse(Area == "all","Overall", ifelse(Area == "urban", "Urban", "Rural")),
         Biomass = Wood + Crop, Clean = Gas + Electricity)%>%
  filter(Area != "Overall")%>%
  dplyr::select(Province, year, Area, all_of(fuel_indx))%>%
  pivot_longer(cols = fuel_indx, names_to = "Fuel")%>%
  mutate(Fuel = factor(Fuel, levels = fuel_indx))

fwrite(survey_data_lf, "result/Table/Input_data/survey_data_lf.csv")

## Prepare input data format for Bayesian model
input_data <- survey_data%>%
  dplyr::rename(whoname = Province, year = Year, area = Population_type, 
                gas = Gas, electricity = Electricity, 
                coal = Coal, wood = Firewood, other = Other_energy)%>%
  mutate(area = ifelse(area == "urban",1, ifelse(area == "rural",2,3)))%>%
  dplyr::select(whoname:other)%>%
  group_by(whoname, year)%>%mutate(id = cur_group_id())%>%ungroup()%>%
  mutate(biomass = wood + other, kerosene = 0, charcoal = 0, cropwaste = 0, 
         dung = 0, total_others = 0, total_biomass = biomass + kerosene,
         solid_2 = total_biomass + coal,clean = gas + electricity,
         total_fuels = gas+electricity+coal+biomass+kerosene)

fwrite(input_data,"result/Table/Input_data/input_survey_data.csv")

## extract unique provincial information
mypy <- pydic(method = 'toneless', dic = c("pinyin2")) ## convert Chinese character into English

prov_info <- distinct(input_data, whoname)%>%
  group_by(whoname)%>%mutate(pindex = cur_group_id())%>%ungroup()%>%
  mutate(Province = ifelse(whoname == "陕西","Shaanxi",
                    ifelse(whoname == "内蒙古","Inner Mongolia",
                    ifelse(whoname == "西藏","Tibet",
                    str_to_title(gsub("_","",py(whoname, dic = mypy)))))))%>%
  mutate(Subarea = ifelse(whoname%in%c("北京","天津","河北","山西","内蒙古","河南"),"N",
                   ifelse(whoname%in%c("上海","江苏","浙江","福建","山东","安徽","江西"),"E",
                   ifelse(whoname%in%c("湖北","湖南","广东","广西","海南"),"S",
                   ifelse(whoname%in%c("陕西","甘肃","青海","宁夏","新疆"),"NW",
                   ifelse(whoname%in%c("吉林","黑龙江","辽宁"),"NE","SW")))))
  )

prov_index <- unique(prov_info$pindex)
fwrite(prov_info,"result/Table/Input_data/province_information.csv")

## county population data derived from satellite
population_county <- readRDS("data/panel_data/Chinese_county_panel_data.Rds")%>%
  dplyr::select(County, Province, Year, Pop_all, Pop_U, Pop_O)%>%
  rename(total = Pop_all, rural = Pop_O, urban = Pop_U)%>%
  mutate(urban = ifelse(urban != 0 ,urban , 1),
         rural = ifelse(rural!= 0, rural, 1),
         total = ifelse(total!= 0, total, 1),
         urban_proportion = urban/total)%>% ## in case of zero value
  rename(whoname = Province)%>%left_join(prov_info)

fwrite(population_county, "result\\Table\\Input_data\\county_population.csv")

population_prov <- readRDS("data/panel_data/Chinese_provincial_panel_data.Rds")%>%
  dplyr::select(Province, Year, Pop_all, Pop_U, Pop_O)%>%
  rename(total = Pop_all, rural = Pop_O, urban = Pop_U)%>%
  mutate(urban = ifelse(urban != 0 ,urban , 1),
         rural = ifelse(rural!= 0, rural, 1),
         total = ifelse(total!= 0, total, 1),
         urban_proportion = urban/total)%>% ## in case of zero value
  rename(whoname = Province)%>%left_join(prov_info)%>%na.omit()

fwrite(population_prov, "result\\Table\\Input_data\\province_population.csv")


##==========Survey data visualization=========
mypy <- pydic(method = 'toneless', dic = c("pinyin2")) ## convert Chinese character into English
## convert survey data to long-format for plot
survey_data_to_plot <- survey_data_lf%>%
  mutate(Area = factor(Area, levels = c("Urban","Rural")))%>%
  rename(whoname = Province)%>%
  left_join(prov_info)
  
## visualize these data
survey_visu <- ggplot(survey_data_to_plot, aes(Province, value, fill = Fuel))+
  facet_grid(year ~ Area)+
  geom_bar(stat = "identity", width = 0.3, position = position_stack())+
  theme_bw()+
  scale_fill_brewer(name = NULL, palette = "Spectral")+
  scale_y_continuous(limits = c(0, 1))+
  guides(fill = guide_legend(byrow = TRUE))+
  labs(x = NULL, y = NULL)+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  theme(legend.position="top",
        legend.key.width=unit(1.5, "cm"),
        legend.key.height=unit(1.5, "cm"),
        legend.text=element_text(colour = "black", size = 45),
        legend.box = "horizontal",
        legend.spacing.y = unit(.8, 'cm'),
        legend.background = element_blank(),
        legend.justification = c(0, 1),
        legend.box.margin=margin(0,0,15,5,unit = "mm"))+ ##legend text size
  theme(strip.text = element_text(size = 40,  color = "black",margin = margin(.1, -.1, .1, -.1, "cm")),
        axis.text.x = element_text(colour = "black", size = 30,margin = unit(c(3, 0, 0, 0), "mm"), angle = 90), 
        axis.text.y = element_text(colour = "black", size = 30,margin = unit(c(0, 3, 0, 0), "mm")))+
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        panel.spacing.x=unit(15, "mm") , panel.spacing.y=unit(15, "mm"))## spacing for different facets

CairoPNG("result\\Fig\\Fuel_use\\Survey_visu_by_province.png",height = 1500,width = 2400)
survey_visu
dev.off()


##==========Validation of population survey data and landscan-based urban and rural population=========
survey_data_wide <- survey_data%>%
  dplyr::select(Province,Year,Population_type,Pop)%>%
  na.omit()%>%pivot_wider(names_from = Population_type, values_from = Pop)%>%
  mutate(total = rural + urban)
##  read in provincial and county-level landscan for overall, urban and rural population
population_prov <-  readRDS(paste0("data/panel_data/Chinese_provincial_panel_data.Rds"))%>%
  dplyr::select(Province, Year, Pop_all, Pop_U, Pop_O)%>%
  rename(total_grid = Pop_all, rural_grid = Pop_O, urban_grid = Pop_U)%>%
  left_join(., survey_data_wide)%>%na.omit()


## population derived from Chinese county yearbook and satellite
population_county <- readRDS(paste0("data/panel_data/Chinese_county_panel_data.Rds"))%>%
  dplyr::select(County, Province, Year, Pop_all, Pop_U, Pop_O,
                Population_SYB, Rural_population)%>%
  rename(total_grid = Pop_all, rural_grid = Pop_O, urban_grid = Pop_U,
         total = Population_SYB, rural = Rural_population)%>%
  group_by(County,Province,Year)%>%
  filter(!any(row_number() > 1))%>% ungroup()%>%
  mutate(rural = ifelse(rural> 300, rural/1e4, rural), ## in case wrong statistical units
         urban_grid = total_grid - rural_grid, urban = total - rural)%>%
  filter(urban >= 0)

  
area_indx <- c("urban","rural","total")
for (i in 1:3) {
  province_to_valid <- population_prov%>%
    rename(y = paste0(area_indx[i]),ypred = paste0(area_indx[i],"_grid"))%>%
    mutate_at(vars(y,ypred), ~ round(.x/1e6, digits = 0))

  county_to_valid = population_county%>%
    rename(y = paste0(area_indx[i]),ypred = paste0(area_indx[i],"_grid"))%>%
    dplyr::select(County, Province, Year, ypred, y)%>%na.omit()%>%
    mutate_at(vars(y,ypred), ~ round(.x/1e4, digits = 0))%>%
    mutate(Year = factor(Year))
  
  cor_value1 <- cor(province_to_valid$y,province_to_valid$ypred)%>%
    formatC(digits = 2,width = 2,format = "f")
  MAE1 = mean(abs(province_to_valid$y - province_to_valid$ypred))%>%
    formatC(digits = 2,width = 2,format = "f")

  cor_value2 = cor(county_to_valid$y,county_to_valid$ypred)%>%
    formatC(digits = 2,width = 2,format = "f")
  
  tag1 <- paste0("(",letters[i],")","<br><br>R = ",cor_value1, "<br><br>MAE = ",MAE1)
  tag2 <- paste0("(",letters[i],")", "<br><br>R = ",cor_value2)
  
  population_plot1 <- ggplot(province_to_valid, aes(y, ypred))+
    geom_abline(slope = 1, intercept = 0, colour="black", size = .8)+
    geom_point(size = 5,  alpha = 0.5,shape = 20,color = "#3182bd")+
    coord_cartesian(xlim = c(0, 1.1e2), ylim = c(0, 1.1e2)) +
    theme(panel.background = element_rect(linewidth = 0.5, fill = NA, colour = 'black'),
          panel.border = element_rect(linewidth = 1.5, colour = "black", fill = NA))+
    labs(x = expression(atop("Observed Population (million)")),
         y = expression("Predicted Population (million)"))+
    theme(strip.background = element_blank())+ ##change the facetwrap background
    tag_facets(tag_pool = tag1, tag_suffix = "",
               position = list(x = 0.004, y = .838, hjust = 0))+
    theme(tagger.panel.tag.text = element_text(color = "black", size = 35),
          tagger.panel.tag.background = element_rect(color = "black", size = 2))+
    theme(axis.text = element_text(size = 28), 
          axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"), size = 40),
          axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"), size = 40))
          
  
  CairoPNG(paste0("result/Fig/Fuel_use/Population_validation_",area_indx[i],".png"),width = 2800,height = 2800,res = 300)
  print(population_plot1)
  dev.off()
  
  # xlim = ylim = c(0, 300)
  # 
  # population_plot2 <- cv.plot(county_to_valid%>%full_join(data.frame(y= max(xlim)+100, ypred = max(ylim)+100)),
  #                             xlim = xlim, ylim = ylim,
  #                             cex = 10, alpha = 0.3, binwidth = c(1.2, 1.2),
  #                             ylimrange = 0.075)+
  #   stat_binhex(bins = 250)+
  #   scale_fill_gradientn(name = expression("Count"), trans = "log10",
  #                        colours = c("#114bf9","#0ffff4","#60f818",
  #                                    "#f2f607","#f2d01e","#f20505"),
  #                        guide = guide_colorbar(frame.colour = "black", frame.linewidth=2,
  #                                               title.position="top", title.hjust=0.5,
  #                                               title.vjust = 1, ticks= F))+
  #   labs(x = expression(atop("Observed Population (10,000)")),
  #        y = expression("Predicted Population (10,000)"))+
  #   theme(legend.position = "right", legend.key.width = unit(.6, "cm"), 
  #         legend.key.height=unit(2.8, "cm"), 
  #         legend.text = element_text(colour = "black", size = rel(2)), ##legend text size
  #         legend.title = element_text(colour = "black", size = rel(2.5)),
  #         legend.margin = margin(r= 10))+
  #   theme(plot.subtitle = element_text(colour = "black", size = rel(2.8),hjust = 0.5))+
  #   theme(axis.text = element_text(size=28), axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm"), size = 35),
  #         axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm"), size = 40))
  # 
  # 
  # CairoPNG(paste0("result/Fig/Fuel_use/Population_validation_county_",area_indx[i],".png"),width = 3200,height = 2800,res = 300)
  # print(population_plot2)
  # dev.off()
}

##  compare the sum of population by year
population_prov%>%
  na.omit()%>%
  group_by(Year)%>%summarise_at(vars(total_grid:total),sum)

population_county%>%
  na.omit()%>%
  group_by(Year)%>%summarise_at(vars(total_grid:urban),sum)

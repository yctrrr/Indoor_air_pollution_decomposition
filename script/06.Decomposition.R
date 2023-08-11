#====================================================#
# retrieve driving forces to the PM2.5 health burden #
#====================================================#

# Core Module Load ----
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
library(data.table)
library(dplyr)
library(tictoc)
library(foreach)
library(terra)
library(furrr)
library(Cairo)
dat_dir <- "D:\\shaoyanchuan\\data\\"

source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Data_Setup.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Model_parallel_provincial.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Functions_parallel2.R"))
source('D:/shaoyanchuan/codebook/function/Grid_visualization.R')
source('D:/shaoyanchuan/codebook/function/Sf_visualization.R')
source('./script/Core.R', encoding = 'UTF8')
source('./script/Decomposition.R', encoding = 'UTF8')

##========Decomposition of fuel use population (grid-level)========
deco_years = c(2000,2007,2013,2020)
METE = c("hd","cd","t2m")
ECO = c("GDP_pc_from_CSMAR","Urban_expenditure_pc","Rural_expenditure_pc")

## county covartaites and related information
county_cov <- fread("result\\Table\\Input_data\\covariates_by_county.csv", encoding = "UTF-8")
county_info <- fread("result\\Table\\Input_data\\county_information.csv", encoding = "UTF-8")%>%
  mutate(index = countyindex)

## fuel use fraction model
model_samples <- readRDS("result\\Rdata\\Fuel_pop\\Fit_all_samples.Rds")

## Decomposition for each period
deco_list <- vector(mode = "list", length = length(deco_years)-1)
names(deco_list) <- deco_years[1:length(deco_years)-1]

Grids = fread('./result/Table/Integrated_PM2.5_exposure/GRID_information.csv')%>%
  mutate(index = countyindex)

Popfy = fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))%>%
  filter(year%in%deco_years)%>%
  left_join(Grids,.)

tic("Fuel use decomposition")
Popf_decomp = Pop_Decomposition(deco_years, county_cov, county_info, Popfy)
toc()

## match on grid
fwrite(Popf_decomp, "result\\Table\\Decomposition\\Fuel_use_population_decomposition_detail.csv")

Popf_decomp <- fread("result\\Table\\Decomposition\\Fuel_use_population_decomposition_detail.csv")
Popf_fuel1 <- Popf_decomp%>%
  dplyr::select(x,y,Area,Fuel,index,year,Start,ECO,MET)%>%
  rename(Popf = Start)

Popf_fuel2 <- Popf_decomp%>%filter(year == 2013)%>%
  dplyr::select(x,y,Area,Fuel,index,year,End)%>%
  rename(Popf = End)%>%mutate(year = 2020, ECO = 0, MET = 0)
Popf_fuel <- full_join(Popf_fuel1, Popf_fuel2)

fwrite(Popf_fuel, "result\\Table\\Decomposition\\Fuel_use_population_decomposition.csv")

##  Nonlinear relationship 
Pop_nonline <- fread(paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'))

Pop_decomp_aggre1 <- fread("result\\Table\\Decomposition\\Fuel_use_population_decomposition_detail.csv")%>%
  left_join(Pop_nonline)%>%group_by(Area,Fuel,year)%>%
  summarise(MET_0 = sum(MET_0*Pop, na.rm = TRUE),
            MET_1 = sum(MET_1*Pop, na.rm = TRUE),
            ECO_0 = sum(ECO_0*Pop, na.rm = TRUE),
            ECO_1 = sum(ECO_1*Pop, na.rm = TRUE),
            MET = sum(MET*Pop, na.rm = TRUE),
            ECO = sum(ECO*Pop, na.rm = TRUE),
            Pop = sum(Pop, na.rm = TRUE))%>%ungroup()%>%arrange(year)

Pop_decomp_aggre <- Pop_decomp_aggre1%>%
  mutate(MET_0 = MET_0/Pop*100, MET_1 = MET_1/Pop*100,
         ECO_0 = ECO_0/Pop*100, ECO_1 = ECO_1/Pop*100,
         MET = MET/Pop*100, ECO = ECO/Pop*100)%>%
  mutate_at(vars(MET_0:ECO), .funs = function(x){x%>%formatC(format = "f", digits = 2)})
Pop_decomp_aggre
fwrite(Pop_decomp_aggre, "result\\Table\\Decomposition\\Fuel_use_population_decomposition_nonlinear.csv")

  
##========Decomposition of PM2.5 related death (grid-level)========
deco_years = c(2000,2007,2013,2020)

## Data load
use_CR('MRBRT')
read_files(
  GRID = './result/Table/Integrated_PM2.5_exposure/GRID_information.csv',
  Pop = paste0('./result/Table/Population/Fuel_use_population_g0.1.csv'),
  PM_real = paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_2000.csv"),
  PM_cf = './data/PM_Ctrl.csv', # PM_cf works only in counter-fact scenario
  MortRate = './data/GBD_incidence_China_2000-2019.csv',
  AgeGroup = './data/GBD_agestructure_China_2000-2017.csv'
)

PM_r <- paste0("./result/Table/Integrated_PM2.5_exposure/Integrated_PM2.5_",deco_years,".csv")%>%
  purrr::map_dfr(~ .x %>%fread%>% 
                   mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)%>%
                   mutate(concentration = matchable(concentration)))%>%
  dplyr::select(-AAP,-HAP_O,-HAP_I)

Pop_grid = Pop%>%filter(year%in%deco_years)%>%dplyr::select(-Popf)
Popf_grid = fread("result\\Table\\Decomposition\\Fuel_use_population_decomposition.csv")%>%
  dplyr::select(-index)%>%
  mutate_at(c("x","y"),.funs = function(num, dgt = 2) num %>% round(dgt) %>% str_c)
Popf_grid%>%filter(year == 2007)

AgeGroupt = fread('./data/IHME_GBD_2019_POP_SYA_2000-2020.csv')%>%
  filter(year%in%deco_years)%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4,include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(sex = tolower(sex), agegroup = as.character(agegroup))

MortRatet = fread('./data/IHME_GBD_2019_Deathrate_mod_2000-2020.csv')%>%
  filter(metric == "Rate", year%in%deco_years)%>%
  dplyr::select(endpoint, sex, agegroup, MortRate, MortRateU, MortRateL, year)%>%
  mutate(agediv = cut(agegroup, breaks = c(0,5,15,65,96),
                      labels = 1:4, include.lowest = TRUE, right = FALSE)%>%as.numeric())%>%
  mutate(endpoint = tolower(endpoint),sex = tolower(sex), agegroup = as.character(agegroup))

RR = (RR_table$MEAN)%>%RR_std()%>%
  mutate(agediv = cut(as.numeric(agegroup), breaks = c(0,5,15,65,96),
                      labels = 1:4,
                      include.lowest = TRUE, right = FALSE)%>%as.numeric())

## Decomposition
for (i in 1:(length(deco_years)-1)) {
  start.y <- deco_years[i]
  end.y <- deco_years[i+1]
  
  cat("=======Decompose PM2.5 Mortality from ",start.y," to ",end.y, "========\n")
  Grids = Grid_info%>%dplyr::select(x,y)
  
  ag = AgeGroupt%>%
    dplyr::select(year, sex, agegroup, AgeStruc, agediv)
  
  mRate = MortRatet%>%
    dplyr::select(year, endpoint, sex, agegroup, MortRate, agediv)


  ## parallel process for all variable combination
  cl <- makePSOCKcluster(10)
  registerDoParallel(cl)
  tic(paste0("Decomposition"))
  
  PM25_Decomptmp <- Decomposition(start.y = start.y, end.y = end.y,
                                  RR = RR, Grids = Grids, Popf = Popf_grid, 
                                  Pop = Pop_grid, AgeGroup = ag,
                                  PM_real = PM_r, mRate = mRate)%>%
    mutate(year = deco_years[i])
  toc()
  stopCluster(cl)

  fwrite(PM25_Decomptmp, paste0("./result/Table/Decomposition/PM25-related_death_decomposition_P",i,".csv"))
}


##====Overall decomposition====
## three periods
PM25_Shape_P = 1:3%>%
  purrr::map_dfr(~ fread(paste0("./result/Table/Decomposition/PM25-related_death_decomposition_P",.x,".csv"))%>%
                   mutate(x = as.numeric(x), y = as.numeric(y), period = .x))

fread(paste0("./result/Table/Decomposition/PM25-related_death_decomposition_P",1,".csv"))%>%
  filter(Decomp == "EXP")%>%
  group_by(Area, Fuel)%>%
  summarise(Mort = sum(Mort))
## whole period
PM25_Shape_4 <- PM25_Shape_P%>%
  group_by(x,y,Decomp,Area,Fuel)%>%
  summarise(Mort = sum(Mort))%>%ungroup()%>%
  mutate(period = 4)
fwrite(PM25_Shape_4,"./result/Table/Decomposition/PM25-related_death_decomposition_P4.csv")
PM25_Shape_by_pop <- full_join(PM25_Shape_P, PM25_Shape_4)

## overall shape value
PM25_Shape = PM25_Shape_by_pop%>%
  group_by(x,y,period,Decomp)%>%
  summarise(Mort = sum(Mort))%>%ungroup()

fwrite(PM25_Shape,"./result/Table/Decomposition/Overall_decomposition.csv")


##=======Visualization of gridded decomposition results=======
library(sf)
polygon <- readOGR(dsn="data/Shp/Province.shp",
                   layer="Province",use_iconv=TRUE, encoding = "UTF-8")

Nline <- readOGR("data/Shp/Nline.shp",
                 layer="Nline",use_iconv=TRUE, encoding = "UTF-8")

polygon_sf <- st_as_sf(polygon)
Nline_sf <- st_as_sf(Nline)


breaks <- c(-Inf, -100, -10, -5, -1, 1, 5, 10, 100, +Inf)
labels <- c("< -100", "-100 - -10", "-10 - -5", "-5 - -1","-1 - 1", "1 - 5","5 - 10", "10 - 100","> 100")
colors <- rev(RColorBrewer::brewer.pal(length(labels),"RdBu"))
name_order = c("MET","ECO","FF","PC","PA","ORF","EXP")
Decomp_vname <- c("MET","ECO","PC","PA","ORF","EXP")

PM25_Shape <- fread("./result/Table/Decomposition/Overall_decomposition.csv")
PM25_Shape_mod = 1:4%>%
  purrr::map_dfr(function(i){
    PM25_Shape%>%filter(period == i)%>%
      dplyr::select(x,y,Decomp,Mort)%>%
      pivot_wider(names_from = "Decomp", values_from = "Mort")%>%
      dfproj()%>%dplyr::select(x,y,Decomp_vname)%>%
      pivot_longer(Decomp_vname)%>%
      mutate(period = i)
  })%>%mutate(value = cut(value, breaks = breaks, labels = labels),
              Decomp = factor(name, levels = name_order))%>%
  mutate_name()


Spatial_Decompplot <- function(data){
  data%>%
    grid_facet_plot(value = "value", facet = "name")+
    facet_wrap("name", nrow = 2, labeller = labeller(groupwrap = label_wrap_gen(20)))+
    geom_sf(data = polygon_sf, fill = "white", alpha = 0.1, linewidth = 0.4)+
    geom_sf(data = Nline_sf, fill = "white", linewidth = 0.6)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    scale_fill_manual(name = expression(paste("Changes  in  Death (#/ 100km"^"2",")")),
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.hjust=0.5,title.vjust=1,
                                           byrow = TRUE, ticks= F))+
    theme(strip.text.x = element_text(size = 32),
          axis.text = element_text(size = 25))+
    theme(panel.grid.major = element_line(colour = alpha("black", 0.4), size = 0.3),
          strip.background = element_blank())+
    theme(legend.title = element_text(colour = "black", size = 40),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(1, "cm"),##legend height
          legend.spacing.y = unit(.8, 'cm'))
}

## visualize the decomposition of different periods
1:4%>%purrr::map(~{
  Spatial_Decomp <- PM25_Shape_mod %>%
    filter(period == .x, Decomp%in%Decomp_vname)%>%
    Spatial_Decompplot
  Spatial_Decomp <- Nanhai_zoom(Spatial_Decomp)
  CairoPNG(paste0("./result/Fig/Decomposition/PM25-related_death_decomposition_P",.x,".png"),
           width = 4000,height = 2200, res = 180)
  print(Spatial_Decomp)
  dev.off()
})

##=======Visualization of gridded decomposition results by areas and fuels=======
fname <- c("Urban","Biomass")
breaks <- c(-Inf, -5, -1, -0.1, -0.005, 0, 0.005, 0.1, 1, 5, +Inf)
colors <- rev(RColorBrewer::brewer.pal(length(breaks)-1,"RdBu"))

PM25_Shape_modgroup <- fread(paste0("./result/Table/Decomposition/PM25-related_death_decomposition_P",1,".csv"))%>%
  filter(Area == fname[1], Fuel == fname[2])%>%
  dplyr::select(x,y,Decomp,Mort)%>%
  pivot_wider(names_from = "Decomp", values_from = "Mort")%>%
  dfproj()%>%dplyr::select(x,y,Decomp_vname)%>%
  pivot_longer(Decomp_vname)%>%
  mutate(value = cut(value, breaks = breaks),
          name = factor(name, levels = name_order))

Spatial_Decomp_group <- Spatial_Decompplot(PM25_Shape_modgroup)


CairoPNG(paste0("./result/Fig/Decomposition/PM25-related_death_decomposition_",
                fname[1],"_",fname[2],".png"),
         width = 4000,height = 2500, res = 180)
print(Spatial_Decomp_group)
dev.off()

##========Changes in decomposed PM2.5 death by province(Table)========
prov_info <- fread("result/Table/Input_data/province_information.csv", encoding = 'UTF-8')
grid_info <- fread('result/Table/Integrated_PM2.5_exposure/GRID_information.csv')
PM25_Shapeall <- fread("./result/Table/Decomposition/PM25-related_death_decomposition_P4.csv")


PM25_Shape_prov <- PM25_Shapeall%>%
  pivot_wider(names_from = "Decomp", values_from = "Mort")%>%
  left_join(grid_info)%>%left_join(prov_info)%>%
  dplyr::select(Start,PC,ECO,MET,PA,ORF,EXP,End,Province,Area,Fuel)%>%
  group_by(Province,Area,Fuel)%>%
  summarise_all(function(x){sum(x)%>%round()})%>%ungroup()%>%
  mutate(Mortc = ((End-Start)/Start*100)%>%formatC(format = "f",digits = 1))

fwrite(PM25_Shape_prov,"result/Table/Decomposition/Decomposition_group_by_province.csv")
PM25_Shape_prov

## explore the rationale of low exposure contribution
IPM_national <- fread("result/Table/Integrated_PM2.5_exposure/Population-weighted_PM2.5_national.csv")%>%
  group_by(year,Area,Fuel)%>%
  summarise(value = sum(value),Pop = mean(Pop), Popnf = mean(Popnf))%>%
  ungroup()

IPM_national%>%
  mutate(value = matchable(value))%>%
  filter(Area == "Urban", Fuel == "Coal")%>%
  mutate(value = ifelse(year == 2000, value[year == 2020], value))%>%
  rename(concentration = value)%>%
  list(RR,mRate,ag)%>%
  reduce(left_join)%>%na.omit()%>%
  dplyr::select(-concentration) %>% dplyr::rename(RR_cf = RR)%>%
  mutate(Mort = Pop * AgeStruc * MortRate * (RR_cf - 1)/ RR_cf / 1e5,.keep = 'unused')%>%
  group_by(year,Area,Fuel)%>%
  summarise(Mort = sum(Mort))%>%
  ungroup()

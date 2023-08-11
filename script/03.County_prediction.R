setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
dat_dir <- "D:\\shaoyanchuan\\data\\"
library(data.table)
library(dplyr)
library(sp)
library(sf)
library(abind)
library(doParallel)
library(nimble)
library(furrr)
library(fasterize)

source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Data_Setup.R"))
# source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Functions_parallel2.R"))
source(paste0(dat_dir,"Population\\Fuel_Use_WHO\\Code\\HAP_Model_parallel_provincial.R"))


#========Prepare county-level covariates========
## read in county boundary
county_sf <- readRDS(paste0(dat_dir,"Shp\\county_China\\county_shp.Rds"))%>%
  dplyr::select(County,Province,City)

county_sp <- as_Spatial(county_sf)

IV = c("GDP_pc_from_CSMAR","hd","cd","t2m",
       "Urban_expenditure_pc","Rural_expenditure_pc")

## Downscale provincial data to county
province_ds <- readRDS(paste0(dat_dir,"Economic/Chinese_panel/Annual/Chinese_provincial_panel_data.Rds"))%>%
  ungroup()%>%dplyr::select(Province, Year, Urban_expenditure_pc, Rural_expenditure_pc)

prov_info <- fread("result\\Table\\Input_data\\Province_information.csv", encoding = 'UTF-8')

county_panel <- readRDS("data/panel_data/Chinese_county_panel_data.Rds")%>%
  left_join(province_ds)%>%
  dplyr::select(all_of(IV), Year, Province, City, County)%>%
  distinct(County, City, Province, Year, .keep_all = TRUE)%>%
  na.omit()

## extract county information
county_info <- county_panel%>%
  distinct(County, City, Province)%>%
  left_join(dplyr::select(prov_info, Province, pindex, Subarea))%>%
  mutate(countyindex = 1:nrow(.))

fwrite(county_info,"result\\Table\\Input_data\\county_information.csv")


county_cov <- county_panel%>%
  left_join(county_info)%>%
  mutate(loc_id = countyindex)%>%
  na.omit()

fwrite(county_cov,"result\\Table\\Input_data\\covariates_by_county.csv")


population_county <- fread("result\\Table\\Input_data\\county_population.csv", encoding = 'UTF-8')

##==========County Prediction for fuel use fraction==========
library(data.table)
plot_years = 2000:2020
# province_samples
fuel_names=c('Solid','Biomass','Wood', 'Clean','Coal','Crop')
area_names=c('Urban','Rural')

model_samples <- readRDS("result\\Rdata\\Fuel_pop\\Fit_all_samples.Rds")

cl <- makePSOCKcluster(20)
registerDoParallel(cl)
county_samples=abind(foreach(i = 1:nrow(county_info))%dopar%{  ## Notice: sequence of county_info is same as the value of countyindex
  library(nimble)
  library(data.table)
  library(dplyr)
  return(fuel_use_output(plot_id = county_info$countyindex[i],
                         plot_index = county_info$pindex[i], cov_data = county_cov,
                         plot_years = plot_years, all_samples = model_samples, 
                         type='survey', output='samples',sampling_bias = TRUE))
  },along=5)
stopCluster(cl)


# Compute 2.5%, 50% and 97.5% quantiles.
county_quantiles <- apply(county_samples, c(2,3,4,5),quantile,c(0.025,0.5,0.975),na.rm = TRUE)

dimnames(county_quantiles) <- list(quantile=c('lower','median','upper'),
                                   year=as.character(plot_years),
                                   Area=area_names, Fuel=fuel_names,
                                   countyindex = county_info$countyindex)

county_quantiles_df <- county_quantiles%>%
  melt(value.name = "pred")%>%
  spread(quantile,pred)%>%
  left_join(county_info)

saveRDS(county_quantiles_df,"result\\Rdata\\Fuel_pop\\Fuel_use_county_prediction.Rds")


##==========Data matching for gridded fuel use population and ambient PM2.5==========
library(terra)
year = 2020
dat_dir <- "D:\\shaoyanchuan\\data\\"
fuel_table <- expand.grid(Area = c("Urban","Rural"),Fuel = c("Clean","Biomass","Coal"))%>%
  mutate(name = paste0(Area,"_",Fuel))

county_info <- read.csv("result\\Table\\Input_data\\County_information.csv", encoding = "UTF-8")
county_sf_join <- county_sf%>%left_join(county_info)

county_quantiles_sf <- read_rds("result\\Rdata\\Fuel_pop\\Fuel_use_county_prediction.Rds")%>%
  filter(Fuel%in%c("Clean","Biomass","Coal"))%>%
  dplyr::select(year, Area, Fuel, Province, City, County, pindex, countyindex, median)%>%
  pivot_wider(names_from = c("Area","Fuel"), values_from = "median")%>%
  mutate(Urban_Clean = 1-Urban_Biomass-Urban_Coal, Rural_Clean = 1-Rural_Biomass-Rural_Coal)%>%
  left_join(county_sf, .)

for (year in 2000:2020) {
  cat("======Prepare raster data for matching in year: ",year,"=======\n")
  y = year
  ##  working with TAP 0.1 grid
  PM25_file <- paste0(dat_dir,"Air_pollution\\TAP\\PM\\",year,"\\",year,"_Annual_Average_PM25.csv")
  if(file.exists(PM25_file)){
    PM25 <- fread(PM25_file)
    coordinates(PM25) <- ~ X_Lon +Y_Lat
    gridded(PM25) <- TRUE
    proj4string(PM25) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
    PM25_raw <- rast(PM25, c("PM2.5"))
    PM25_R0.1 <- PM25_raw$PM2.5
    names(PM25_R0.1) <- "APM25"
  }
  
  ## prepare urban population data
  Urban_Pop_file <- paste0(dat_dir, "Population\\RU_Pop\\ESA-CCL\\Urban_Pop_0.00833g_",year,".tif")
  Urban_Pop <- rast(Urban_Pop_file)
  Urban_area <- ifel(Urban_Pop[[2]] > 0, 1, 0)
  names(Urban_area) <- "U_area"
  Urban_area$Pop_all <- Urban_Pop[[1]]
  Urban_area$Pop_U <- ifel(Urban_area$U_area == 1,Urban_area$Pop_all, 0)
  Urban_area$Pop_O <- ifel(Urban_area$U_area == 0,Urban_area$Pop_all, 0)
  
  ##  resample population data
  Urban_area0.1 <- resample(Urban_area, PM25_R0.1, method = "sum")

  Urban_area0.1$U_popf <- Urban_area0.1$Pop_U/Urban_area0.1$Pop_all
  
  Area_info_1 <- fasterize(county_sf_join, raster(PM25_R0.1), field = "pindex",fun = "max")%>%
    rast%>%ifel(. > 0, ., NA)
  names(Area_info_1) = "pindex"

  Area_info_2 <- fasterize(county_sf_join, raster(PM25_R0.1), field = "countyindex",fun = "max")%>%
    rast%>%ifel(. > 0, ., NA)
  names(Area_info_2) = "countyindex"
  
  ## match heating degree days
  hdd_file <- paste0(dat_dir,"ERA5\\annual\\",year,"\\era5-heating-degree-day-g0.25-",year,".Rds")
  if(file.exists(hdd_file)){
    hdd_df <- readRDS(hdd_file)%>%
      dplyr::select(-t2m)
    
    coordinates(hdd_df) <- ~ x + y
    gridded(hdd_df) <- TRUE
    proj4string(hdd_df) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
    hdd_R0.1 <- rast(hdd_df)[[1:4]]%>%
      resample(PM25_R0.1, method = "average")
  }

  ## match meteorological variables
  era5_file <- paste0(dat_dir,"ERA5\\annual\\",year,"\\era5-reanalysis-single-an-g0.25-",year,".grd")
  if(file.exists(era5_file)){
    era5_R0.1 <- rast(era5_file)[[c("u10","v10","t2m","d2m")]]%>%
      resample(PM25_R0.1, method = "average")
  }
  
  ##  assign county fuel use fraction to raster
  fuel_fraction_raster <- fuel_table$name%>%
    purrr::map(function(x){
      fuel_rast = county_quantiles_sf%>%
        filter(year == y)%>%
        fasterize(raster(PM25_R0.1), field = x)%>%rast
      names(fuel_rast) = x
      fuel_rast})%>%purrr::reduce(c)

  ##  Match raster on 0.1*0.1 grid
  match_rast0.1 <- c(Urban_area0.1, PM25_R0.1, hdd_R0.1, era5_R0.1, fuel_fraction_raster, Area_info_1, Area_info_2)
  writeRaster(match_rast0.1, paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",year,".tif"),overwrite = TRUE)
}

##  extract fuel use population table
for (year in 2000:2020) {
  # match_rast1km <- rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.0083_",y,".tif"))
  cat("========Process raster data Year:",year,"========\n")
  
  fuel_pop_df <- rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",year,".tif"))%>%
    as.data.frame(xy = TRUE, na.rm = FALSE)%>%
    dplyr::select(x,y,Pop_U,Pop_O,Urban_Clean:Rural_Coal)%>%na.omit()%>%
    pivot_longer(Urban_Clean:Rural_Coal, names_to = c("Area","Fuel"),
                 names_sep = "_", values_to = "Popf")%>%
    mutate(Pop = ifelse(Area == "Urban", Pop_U, Pop_O), year = year)%>%
    dplyr::select(-Pop_U,-Pop_O)

  MET_df <- rast(paste0("result\\Tif\\Match_data\\Matching_covariate_raster_g0.1_",year,".tif"))%>%
    as.data.frame(xy = TRUE, na.rm = FALSE)%>%
    dplyr::select(x,y,hd,cd,hd_count,cd_count,t2m)%>%na.omit()%>%
    mutate(year = year)
  
  fwrite(fuel_pop_df, paste0("result\\Table\\Population\\Fuel_use_population_g0.1_",year,".csv"))
  fwrite(MET_df, paste0("result\\Table\\Covariates\\Temperature_g0.1_",year,".csv"))
}
## Population
2000:2020%>%
  map_dfr(~ fread(paste0("result\\Table\\Population\\Fuel_use_population_g0.1_",.x,".csv")))%>%
  fwrite(paste0("result\\Table\\Population\\Fuel_use_population_g0.1.csv"))
## Temperature
2000:2020%>%
  map_dfr(~ fread(paste0("result\\Table\\Covariates\\Temperature_g0.1_",.x,".csv")))%>%
  fwrite(paste0("result\\Table\\Covariates\\Temperature_g0.1.csv"))

##=======Visualization of fuel use population========
## fuel use plot function
Fueluse_plot <- function(data ,value = "Pop", colors, title = NULL){
  if(is.null(title)){
    title = expression(paste("Fuel use population (#/ 100km"^"2",")"))
  }
  
  grid_plot(data = data, label = NULL, value = value)+
    facet_grid(Area ~ Fuel, labeller = labeller(groupwrap = label_wrap_gen(10)))+
    scale_fill_manual(name = title,
                      values = colors,
                      guide = guide_legend(frame.colour = "black",frame.linewidth=2,
                                           title.position="right",title.hjust=0.5,title.vjust=1,
                                           byrow = TRUE, ticks= F))+
    geom_sf(data = polygon_sf, fill = "white", alpha = 0.1, linewidth = 0.4)+
    geom_sf(data = Nline_sf, fill = "white", linewidth = 0.6)+
    coord_sf()+xlim(-2879760, 2300893)+ylim(1800000,5940000)+
    theme(strip.background = element_blank())+ ##change the facetwrap background
    theme(strip.text.x  = element_text(size = 32),
          strip.text.y  = element_text(size = 32),
          axis.text = element_text(size = 25))+
    theme(panel.grid.major = element_line(colour = alpha("black", 0.4), size = 0.3))+
    theme(legend.title = element_text(colour = "black", size = 30))
}


##=======Change in fuel use population========
fuel_pop_aggrmulti <- fread(paste0('result\\Table\\Population\\Fuel_use_population_g0.1.csv'))%>%
  filter(year%in%c(2000,2020))%>%mutate(Pop = Pop*Popf)
breaks <- c(-Inf, -10000, -5000, -1000, -500, -100, 100, 500, 1000, 5000, 10000, +Inf)

labels <- c("< -10000", "-10000 - -5000", "-5000 - -1000", "-1000 - -500","-500 - -100",
            "-100 - 100", "100 - 500","500 - 1000", "1000 - 5000", "5000 - 10000","> 10000")
colors <- colorspace::diverging_hcl(length(labels), palette = "Blue-Red")

pop_aggrmulti <- 1:nrow(fuel_table)%>%
  map_dfr(~ fuel_pop_aggrmulti%>%
            filter(Area == fuel_table$Area[.x],Fuel == fuel_table$Fuel[.x])%>%
            dplyr::select(x,y,year,Pop)%>%
            pivot_wider(names_from = year, values_from = Pop)%>%
            mutate(Popc = `2020` - `2000`,.keep = "unused")%>%dfproj%>%
            mutate(Area = fuel_table$Area[.x],Fuel = fuel_table$Fuel[.x]))%>%
  mutate(Popc = cut(Popc, breaks, labels))

pop_aggrmulti_plot <- Fueluse_plot(pop_aggrmulti, value = "Popc", colors = colors)
pop_aggrmulti_plot <- Nanhai_zoom(pop_aggrmulti_plot)

Cairo::CairoPNG(paste0("result/Fig/Fuel_use/Fuel_use_population_Change_2000-2020.png"),
                width = 4000,height = 2500, res = 180)
print(pop_aggrmulti_plot)
dev.off()

## Fuel use fraction change
breaks <- c(-1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1)
labels <- c("-1 - -0.5", "-0.5 - -0.2", "-0.2 - -0.1", "-0.1 - 0",
            "0 - 0.1","0.1 - 0.2", "0.2 - 0.5", "0.5 - 1")
colors <- rev(brewer.pal(n = length(labels),"RdBu"))

popf_aggrmulti <- 1:nrow(fuel_table)%>%
  map_dfr(~ fuel_pop_aggrmulti%>%
            filter(Area == fuel_table$Area[.x],Fuel == fuel_table$Fuel[.x])%>%
            dplyr::select(x,y,year,Popf)%>%
            pivot_wider(names_from = year, values_from = Popf)%>%
            mutate(Popfc = `2020` - `2000`,.keep = "unused")%>%dfproj(method = "near")%>%
            mutate(Area = fuel_table$Area[.x],Fuel = fuel_table$Fuel[.x]))%>%
  mutate(Popfc = cut(Popfc, breaks, labels))

popf_aggrmulti_plot <- Fueluse_plot(popf_aggrmulti, value = "Popfc", 
                                    colors = colors, title = "Fuel use fraction change")
popf_aggrmulti_plot <- Nanhai_zoom(popf_aggrmulti_plot)
Cairo::CairoPNG(paste0("result/Fig/Fuel_use/Fuel_use_fraction_Change_2000-2020.png"),
                width = 4000,height = 2500, res = 180)
print(popf_aggrmulti_plot)
dev.off()

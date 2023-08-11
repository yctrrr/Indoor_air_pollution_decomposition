# PM25 Health Impact Calc Core
# garbage code adapted from Yifan LIU 2021/05/15
# Modify by Yanchuan Shao 2022/04/09

library(tidyverse)
library(writexl)
library(readxl)
library(tictoc)
library(data.table)

tell_CR <- function(CR = .CR_fun) return(
  if (CR %>% str_detect('IER')) str_c(CR)
  else if (CR %in% c('5COD', 'NCD+LRI'))  str_c('GEMM', CR, sep = '_')
  else if (CR == 'MRBRT') str_c(CR)
)

use_CR <- function(CR_fun) {
  assign(".CR_fun", CR_fun, envir = globalenv())
  if (str_detect(CR_fun, "IER|NCD\\+LRI|5COD|MRBRT")) {
    # cat(str_glue("C-R function \"{tell_CR(CR_fun)}\" is set as the default methodology\n"))
  } else {
    warning(str_glue("an exogenous C-R function \"{CR_fun}\" was specified, \\
                     please provide a corrosponding `RR_table` after `read_file()`"))
  }
}

matchable <- function(num, dgt = 1) ifelse(num <=300, num %>% round(dgt) %>% str_c, round(num)%>% str_c)

# dataload module

read_files <- function(
    GRID = './Data/GRID_information_sample.xlsx',
    Pop = './Data/GridPop_sample.xlsx',
    PM_real = './Data/GridPM25_sample.xlsx',
    PM_cf = './Data/PM_Ctrl.csv', # PM_cf works only in counter-fact scenario
    MortRate = './Data/GBD_incidence_China_2000-2019.csv',
    AgeGroup = './Data/GBD_agestructure_China_2000-2017.csv') {
  
  fuse_read <- function(filename) filename %>% {
    if (str_detect(., 'csv$')) fread(.)
    else if (str_detect(., 'xlsx$')) read_xlsx(.)
  }

  assign(
    'Grid_info', 
    envir = globalenv(), 
    fuse_read(GRID) %>% 
      mutate(across(where(is.numeric) & x:y, matchable, dgt = 2))
  )
  
  assign(
    "Pop", envir = globalenv(),
    fuse_read(Pop) %>% 
      mutate(across(where(is.numeric) & x:y, matchable, dgt = 2))
  )
  
  assign(
    'PM_real', envir = globalenv(), 
    fuse_read(PM_real) %>% mutate(
      across(where(is.numeric) & x:y, matchable, dgt = 2),
      across(where(is.numeric) & matches('^\\d{4}'), matchable, dgt = 1)
    )
  )
  
  # Specify UNREAL PM2.5 data, used for only counter-fact scenario.
  
  assign(
    'PM_cf', 
    envir = globalenv(),
    if (file.exists(PM_cf)) {
      fuse_read(PM_cf) %>% mutate(
        across(where(is.numeric) & x:y, matchable, dgt = 2),
        across(where(is.numeric) & matches('^\\d{4}'), matchable, dgt = 1)
      )
    } else NULL
  )
  
  assign(
    'MortRate', 
    envir = globalenv(),
    fuse_read(MortRate) %>% pivot_longer(
      cols = c(-year,-endpoint), names_to = 'agegroup',values_to = 'MortRate'
    ) %>% pivot_wider(
      names_from = 'year', values_from = 'MortRate'
    )
  )
  
  assign(
    "AgeGroup", 
    envir = globalenv(),
    fuse_read(AgeGroup) %>% pivot_longer(
      cols = -`year`, names_to = 'agegroup', values_to = 'AgeStruc'
    ) %>% pivot_wider(
      names_from = 'year', values_from = 'AgeStruc'
    ) %>% mutate(across(-agegroup, prop.table))
  )

  CR_file <- if (!exists('.CR_fun', envir = globalenv())) {
    warning("Please specify a C-R function with `use_CR()`")
  } else if (.CR_fun == 'MRBRT') {
    './Data/RR_index/MRBRT2019_Lookup_Table_LYF220601.xlsx'
  } else if (.CR_fun %in% c('NCD+LRI', '5COD')) {
    './Data/RR_index/GEMM_Lookup_Table_Build_220601.xlsx'
  } else if (.CR_fun %in% c('IER', 'IER2017')) {
    './Data/RR_index/IER2017_Lookup_Table_Build_220601.xlsx'
  } else if (.CR_fun == 'IER2015') {
    './Data/RR_index/IER2015_Lookup_Table_Build_220601.xlsx'
  } else if (.CR_fun == 'IER2013') {
    './Data/RR_index/IER2013_Lookup_Table_Build_220601.xlsx'
  } else if (.CR_fun == 'IER2010') {
    './Data/RR_index/IER2010_Lookup_Table_Build_220601.xlsx'
  } else {
    NA_character_
  }

  assign(
    "RR_table", envir = globalenv(),
    if (is.na(CR_file)) NA_character_
    else expand_grid(excel_sheets(CR_file), CR_file) %>% deframe %>% imap(
      ~ read_excel(.x, sheet = .y) %>% 
        mutate(across(where(is.numeric) & concentration, matchable, dgt = 1)))
  )
}

mortrate_std <- function(x, year) x %>%
  mutate(endpoint = tolower(endpoint)) %>%
  dplyr::select(endpoint, agegroup, MortRate = !!year)

RR_std <- function(RR_index) {
  
  CR <-  if (!exists('.CR_fun', envir = globalenv())) {
    stop("Please specify a C-R function with `use_CR()`")
  } else get(".CR_fun", envir = globalenv())
  
  RR_tbl <- RR_index %>% pivot_longer(
    cols = -concentration,
    values_to = "RR",
    names_to = c("endpoint", "agegroup"),
    names_sep = '_'
  ) %>% mutate(endpoint = tolower(endpoint))
  
  RR_reshape <- if (CR == '5COD') {
    expand_grid(
      concentration = RR_index %>% pull(concentration),
      endpoint = c('copd', 'ihd', 'lc', 'lri', 'stroke'),
      agegroup = c('ALL', seq(25, 95, 5) %>% matchable(0))
    ) %>% left_join(RR_tbl) %>% 
      group_by(concentration, endpoint) %>% fill(RR) %>% ungroup %>% filter(agegroup != 'ALL')
    
  } else if (CR == 'NCD+LRI') {
    expand_grid(
      concentration = RR_index %>% pull(concentration),
      endpoint = c('ncd+lri'),
      agegroup = c('ALL', seq(25, 95, 5) %>% matchable(0))
    ) %>% left_join(RR_tbl) %>% 
      group_by(concentration, endpoint) %>% fill(RR) %>% ungroup %>% filter(agegroup != 'ALL')
    
  } else if (CR %>% str_detect('IER')) {
    expand_grid(
      concentration = RR_index %>% pull(concentration),
      endpoint = c('copd', 'ihd', 'lc', 'stroke', 'lri'),
      agegroup = c('ALL', seq(0, 95, 5) %>% matchable(0))
    ) %>% left_join(RR_tbl) %>% 
      group_by(concentration, endpoint) %>% fill(RR) %>% filter(agegroup != 'ALL') %>% 
      filter(
        (endpoint  %>% str_detect('copd|ihd|lc|stroke') & as.integer(agegroup) >= 25) | 
          endpoint == 'lri') %>% ungroup
  } else if (CR == 'MRBRT') {
    expand_grid(
      concentration = RR_index %>% pull(concentration),
      endpoint = c('copd', 'dm', 'ihd', 'lc', 'lri', 'stroke'),
      agegroup = c('ALL', seq(0, 95, 5) %>% matchable(0))
    ) %>% left_join(RR_tbl) %>% 
      group_by(concentration, endpoint) %>% fill(RR) %>% filter(agegroup != 'ALL') %>% 
      filter(
        (endpoint %>% str_detect('copd|dm|ihd|lc|stroke') & as.integer(agegroup) >= 25) | 
          endpoint == 'lri') %>% ungroup
  }
  
  return(RR_reshape)
}

#' Calculate gridded PM2.5 attributed mortality
#'
#' @param Grids a vector of grid coords
#' @param PM_r a 2-column `data.frame` stores real PM2.5 concentration of each grid
#' @param PM_c a 2-column `data.frame` stores virtual PM2.5 concentration of each grid, equals PM_r at default
#' @param ag proportions of 20 age-groups inside the population structure
#' @param mRate the mortality rates of each endpoints and each age group
#' @param pop a 2-column dataframe stores population volume of each grid
#' @param RR the lookup-table of Concentration-Response functions of PM2.5 exposure
#' @param CR a character string instructs the name of the C-R function
#'
#' @return a table of death estimates for each endpoint & age-groups(columns) for every grids(rows)
#'
#' @examples
Mortality <- function(Grids, PM_r, PM_c = NULL, ag, mRate, pop, RR) {

  if (is.null(PM_c)) PM_c <- PM_r

  RR_tbl <- RR

  result <- list(Grids, PM_c, pop, mRate, ag) %>% 
    reduce(left_join)%>%suppressMessages()
  

  result2 <- result%>%list(RR_tbl) %>%
    reduce(left_join)%>%na.omit()%>%
    dplyr::select(-concentration) %>% dplyr::rename(RR_cf = RR)%>%
    suppressMessages()
  # result2%>%nrow
  # result%>%nrow

  result3 <- result2%>%
    # mutate(Mort = Pop * AgeStruc * MortRate * (RR_cf - 1) / PWRR_real / 1e5,.keep = 'unused')%>%
    mutate(Mort = Pop * AgeStruc * MortRate * (RR_cf - 1)/ RR_cf / 1e5,.keep = 'unused')%>%
    group_by(x,y,Area,Fuel)%>%
    dplyr::summarise(Mort = sum(Mort))%>%
    ungroup()%>%
    suppressMessages()
  return(result3)
}

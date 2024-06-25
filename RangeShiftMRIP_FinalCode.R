##
## Date Created: 2024-06-12
##
## Copyright (c) Rachel Roday, 2024 and Rieligh Hudock, 2024
## Email: rroday@udel.edu and hudockr@udel.edu
##
---------------------------
  ##
  ## Notes: Range Shift Code 
  ##         Collaborators : Ed Hale, Aaron Carlisle, Rileigh Hudock, Dan Millea,
  ##                         Taylor Hoffman, Willa Lane, Caitlin Wilson
  ##       Note to self:  RangeDF4 is the same as  MRIP
  ## ---------------------------

######################## Working Directory and Packages ########################
# Set working directory
setwd("C:/Users/RER/Documents/Masters UD/Range Shift")

# Load libraries 

library(ggplot2)
library(tidyverse)
library(dplyr)
library(sjlabelled)
library(sjmisc)
library(geofacet)
library(boot)
library(zoo)
library(roll)
library(ggiraph)
library(readr)
library(lubridate)
library(broom)
library(FSA)
library(tseries)
library(forecast)
library(Kendall)
library(trend)
library(stats)
library(scales)
library(viridis)
library(RColorBrewer)

######################## Loading Data  #########################################

# Sea surface temperature from raster cells in ERDDAP
SST <- read_csv("State SST.csv")  %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate(Year = year(time), 
         Month = month(date(time), label = T),
         Month_Yr = paste(Year, Month, sep = " "))

# MRIP downloaded data
RangeDF <- read.csv("MAST630 Project - Updated Data.csv", header = T)

#State coordinates
Coords <- read.csv("State_Coast_Coords.csv", header = T)

# NEFSC Spring Trawl Catch Data
catch.data <- read.csv("22561_UNION_FSCS_SVCAT.csv", header = T)

# NEFSC Spring Trawl Station Data
station.data <- read.csv("22561_UNION_FSCS_SVSTA.csv", header = T)

#MRIP <- read.csv("~/Desktop/Range Shift Paper/Bigelow Trawl Survey/22561_NEFSCSpringFisheriesIndependentBottomTrawlData/MRIP.DF.csv")
#catch22 <- read.csv("catch22.csv")


# Species Labellers
Spp = c("Atlantic Mackerel" = "Atlantic mackerel", "Black Sea Bass" = "Black sea bass", 
        "Bluefish" = "Bluefish", "Cobia" = "Cobia", "Florida Pompano" = "Florida pompano",
        "Summer Flounder" = "Summer flounder","Winter Flounder" = "Winter flounder", 
        "Yellowfin Tuna" = "Yellowfin tuna")

# Create formatter function
formatter <- function(...){
  function(x) format(round(x, 1), ...)
}

######################## Tidying  Data  ########################################

## Tidying MRIP dataframe
colnames(RangeDF) = c("Year", "State", "Species", "Catch", "PSE", "TE",
                      "TEPSE", "CPUE")

RangeDF[RangeDF == "." & !is.na(RangeDF)] <- NA
RangeDF[RangeDF == "" & !is.na(RangeDF)] <- NA
RangeDF <- type.convert(RangeDF, as.is = TRUE)
RangeDF$CPUE <- as.numeric(RangeDF$CPUE)

RangeDF1 <- RangeDF %>%
  mutate(State = recode(State, "DELAWARE" = "Delaware",
                        "GEORGIA" = "Georgia", 
                        "MARYLAND" = "Maryland",
                        "MASSACHUSETTS" = "Massachusetts",
                        "Massachussets" = "Massachusetts",
                        "NEW JERSEY" = "New Jersey", 
                        "NORTH CAROLINA" = "North Carolina",
                        "SOUTH CAROLINA" = "South Carolina",
                        "VIRGINIA" = "Virginia",
                        "Connecticut " = "Connecticut"),
         Species = recode(Species, "ATLANTIC MACKEREL" = "Atlantic Mackerel",
                          "BLACK SEA BASS" = "Black Sea Bass",
                          "BLUEFISH" = "Bluefish",
                          "COBIA" = "Cobia",
                          "FLORIDA POMPANO" = "Florida Pompano",
                          "SUMMER FLOUNDER" = "Summer Flounder",
                          "WINTER FLOUNDER" = "Winter Flounder",
                          "YELLOWFIN TUNA" = "Yellowfin Tuna")) %>%
  drop_na(State, Species) %>%
  group_by(State, Species) %>%
  arrange(State, Species, Year) %>%
  mutate(Bin1 = cut(Year, breaks=c(1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015,2020,2025)),
         Bin2 = cut(Year, breaks = c(1980, 1990, 2000, 2010, 2020, 2030)),
         # I dont like the rolling averages and I will only use binned averages from here on out
         RollMean5yr = zoo::rollmean(CPUE, k = 5, fill = NA),
         RollMean10yr = zoo::rollmean(CPUE, k = 10, fill = NA),
         RollSD5yr = sqrt(sum((CPUE-RollMean5yr)^2)/5),
         ROLL = roll::roll_mean(CPUE, width = 5, min_obs = 4),
         roll2 = rollapplyr(1:42, 5, mean, fill = NA)) %>%
  ungroup() %>%
  group_by(State, Species, Bin1) %>%
  mutate(Avg5 = mean(CPUE),
         SD5 = sd(CPUE)) %>%
  ungroup() %>%
  group_by(State, Species, Bin2) %>%
  mutate(Avg10 = mean(CPUE),
         SD10 = sd(CPUE)) 

# Join state coordiantes to MRIP data frame, now called RangeDF2
RangeDF2 <- left_join(RangeDF1, Coords, by = "State") %>%
  mutate(State = factor(State, levels = c("Maine", "New Hampshire", "Massachusetts", "Rhode Island", "Connecticut", "New York", "New Jersey", "Delaware",
                                          "Maryland", "Virginia", "North Carolina", "South Carolina", "Georgia")))

# Subset MRIP data to exclude southern states
RangeDF.subset <- RangeDF2 %>%
  filter(State != "South Carolina",
         State != "Georgia")

## Tidying NEFSC data

#remove unimportant columns
catch.data <- catch.data %>% select(-TOW, -STATION, -STATUS_CODE, -ID, -CATCHSEX, 
                                    -CATCH_COMMENT, -X, -X.1, -X.2, -X.3)

# remove unimportant columns from station data
station.data <- station.data %>% select(-TOW, -STATION, -STATUS_CODE, -ID, -GEARCOND, -SHG, 
                                        -GEAR_COMMENT, -ACQUISITION_CODE, -ACQUISITION_COMMENT, 
                                        -TOGA, -SVVESSEL, -CRUNUM, -BEGIN_GMT_TOWDATE, -END_GMT_TOWDATE, 
                                        -EST_JULIAN_DAY, -GMT_YEAR, -GMT_MONTH, -GMT_DAY, -GMT_JULIAN_DAY, 
                                        -GMT_TIME, -BEGEKVLOG, -ENDEKVLOG, -BEGLAT, -BEGLON, 
                                        -ENDLAT, -ENDLON, -LORS1, -LORS2, -LORE1, -LORE2, -CABLE, 
                                        -PITCH, -HEADING, -COURSE, -RPM, -DOPDISTB, -DOPDISTW, 
                                        -DESSPEED, -GEARID, -DOORID, -CLOUD, -WINDDIR, -SWELLDIR, 
                                        -SWELLHGT, -TRASHAMT, -TRASHSHL, -TRASHBIO, -TRASHSUB, -XBT, 
                                        -FULD, -NO_DETAIL_SVSPP, -BOTSPEED, -WATCH_CHIEF_COMMENTS, 
                                        -STATION_COMMENTS, -HABITAT_COMMENTS) %>%
  filter(!is.na(DECDEG_BEGLAT))

# # removing unnecessary columns from MRIP data frame
# MRIP <- MRIP %>% select(-X) %>% mutate(Data.Type = "MRIP") %>% rename(range = range.lat)

# filtering catch data and redoing scientific names
filtered.catch <- catch.data %>%
  separate(SCIENTIFIC_NAME, into = c("genus","species","common", "common2", "common3")) %>%
  mutate(sci = paste(genus, species), 
         com = paste(common, common2, common3)) %>%
  filter(sci %in% c("Centropristis striata", "Pomatomus saltatrix", "Scomber scombrus", "Rachycentron canadum", "Paralichthys dentatus", "Pseudopleuronectes americanus"))

unfiltered.catch <- catch.data %>%
  separate(SCIENTIFIC_NAME, into = c("genus","species","common", "common2", "common3")) %>%
  mutate(sci = paste(genus, species), 
         com = paste(common, common2, common3)) 

##Organized stations and catch, join data
catch <- filtered.catch %>%
  group_by(CRUISE6, CRUISE, STRATUM, com) %>%
  summarise(avg.catch = mean(EXPCATCHNUM),
            exp.catch = EXPCATCHNUM) %>%
  filter(!any(is.na(exp.catch)))

stations <- station.data %>%
  group_by(CRUISE6, CRUISE, STRATUM) %>%
  summarise(surf.temp = mean(SURFTEMP), 
            bkt.temp = mean(BKTTEMP),
            bot.temp = mean(BOTTEMP),
            year = EST_YEAR,
            month = EST_MONTH,
            day = EST_DAY,
            lat = mean(DECDEG_BEGLAT),
            long = mean(DECDEG_BEGLON),
            range.lat = DECDEG_BEGLAT) %>%
  mutate(qt = surf.temp - bkt.temp) %>%
  distinct(CRUISE6, CRUISE, STRATUM, .keep_all = T)


# joining catch data and station data
NEFSC <- left_join(catch, stations) %>%
  mutate(xm = lat * avg.catch) %>%
  summarise(xbar = (sum(xm))/sum(avg.catch),
            com = com,
            year = year,
            month = month,
            day = day,
            range.lat = range.lat,
            exp.catch = exp.catch,
            avg.catch = avg.catch) %>%
  mutate(date = paste(month,day, year, sep = "-"),
         date2 = mdy(date))

##Filtering Bigelow trawl locations by state 
NEFSC2 <- NEFSC %>% 
  mutate(state = case_when((range.lat > 30.1 & range.lat < 32.039091) ~ "GEORGIA",
                           (range.lat > 32.039091 & range.lat < 33.850616) ~ "SOUTH CAROLINA",
                           (range.lat > 33.850616 & range.lat < 36.545642) ~ "NORTH CAROLINA",
                           (range.lat > 36.545642 & range.lat < 38.014752) ~ "VIRGINIA",
                           (range.lat > 38.014752 & range.lat < 38.450555) ~ "MARYLAND",
                           (range.lat > 38.450555 & range.lat < 38.872217) ~ "DELAWARE",
                           (range.lat > 38.872217 & range.lat < 40.476778) ~ "NEW JERSEY",
                           (range.lat > 40.476778 & range.lat < 41.005902) ~ "NEW YORK",
                           (range.lat > 41.005902 & range.lat < 41.371745) ~ "CONNECTICUT",
                           (range.lat > 41.371745 & range.lat < 41.493239) ~ "RHODE ISLAND",
                           (range.lat > 41.493239 & range.lat < 42.871019) ~ "MASSACHUSETTS",
                           (range.lat > 42.871019 & range.lat < 43.04283) ~ "NEW HAMPSHIRE",
                           (range.lat > 43.04283 & range.lat < 44.823452) ~ "MAINE",
                           TRUE ~ "NA"),  # Use TRUE as the condition for the default case
         center_lat = case_when((state == "GEORGIA") ~ 31.399899,
                                (state == "SOUTH CAROLINA") ~ 32.892435,
                                (state == "NORTH CAROLINA") ~ 34.762265,
                                (state == "VIRGINIA") ~ 37.36721,
                                (state == "MARYLAND") ~ 38.230127,
                                (state == "DELAWARE") ~ 38.633806,
                                (state == "NEW JERSEY") ~ 39.821385,
                                (state == "NEW YORK") ~ 40.704224,
                                (state == "CONNECTICUT") ~ 41.259409,
                                (state == "RHODE ISLAND") ~ 41.390066,
                                (state == "MASSACHUSETTS") ~ 41.971489,
                                (state == "NEW HAMPSHIRE") ~ 42.96563,
                                (state == "MAINE") ~ 44.013106,
                                TRUE ~ NA_real_)) %>% # Use TRUE as the condition for the default case and NA_real_ as the default value
  mutate(Data.Type = "Bigelow") %>%
  mutate(com2 = case_when((com == "Atlantic mackerel ") ~ "Atlantic Mackerel",
                          (com == "black sea bass") ~ "Black Sea Bass",
                          (com == "summer flounder ") ~ "Summer Flounder",
                          (com == "winter flounder ") ~ "Winter Flounder",
                          (com == "cobia  NA") ~ "Cobia",
                          (com == "bluefish  NA") ~ "Bluefish",
                          TRUE ~ "NA")) %>%
  filter(year >= 1981, year <= 2022)


# creating a new df with the range to graph later
RANGE <- NEFSC2 %>% 
  group_by(year, com2) %>%
  summarise(range = max(center_lat)- min(center_lat),
            min = min(center_lat),
            max = max(center_lat)) %>%
  filter(!is.na(year)) 


# NEFSC Center of Distribution
##Total CPUE/Weighted latitude
DIST <- NEFSC2 %>% group_by(com2,year) %>%
  mutate(range = max(center_lat)- min(center_lat),
         min = min(center_lat),
         max = max(center_lat),
         total_CPUE = sum(avg.catch),
         xm = center_lat * total_CPUE) %>%
  filter(!is.na(year)) %>%
  summarise(xbar = (sum(xm))/(sum(total_CPUE)),
            com2 = com2,
            year = year,
            month = month,
            day = day,
            range.lat = range.lat,
            exp.catch = exp.catch,
            avg.catch = avg.catch,
            range = range,
            center_lat = center_lat) %>%
  mutate(date = paste(month, day, year, sep = "/"),
         date2 = mdy(date)) %>%
  mutate(Data.Type = "Bigelow") %>%
  distinct(year, com2, xbar, range, Data.Type)

DIST2 <- left_join(RANGE, DIST, by = c("year", "com2", "range")) %>%
  filter(!com2 == "Cobia")
#

######################## Visualizing Raw MRIP Data  ############################

# Total catch by state
ggplot(RangeDF2, aes(x = Species, y = Catch)) +
  geom_bar(stat ="identity") +
  facet_wrap(~State) +
  labs(x = "", 
       y= "Catch", 
       title ="",
       caption = "Source: MRIP| Viz: @rachel_roday") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Total CPUE over timespan for each specfies in each state
ggplot(RangeDF2, aes(x = State, y = CPUE)) +
  geom_bar(stat ="identity") +
  facet_wrap(~Species, ncol = 4) +
  labs(x = "", 
       y= "Total CPUE 1980 - 2022", 
       title ="",
       caption = "Source: MRIP| Viz: @rachel_roday") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# CPUE over time by state
ggplot(RangeDF2, aes(x=Year, y= RollMean10yr))+
  geom_line(aes(color = Species))+
  facet_wrap(~State, ncol =2,  scales = "free", strip.position = "right")+
  theme_bw() 

# CPUE over time by Species
ggplot(RangeDF2, aes(x=Year, y= RollMean10yr))+
  geom_line(aes(color = State))+
  facet_wrap(~Species, scales = "free")+
  theme_bw()

# CPUE over state by time 
ggplot(RangeDF2, aes(x=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                "Connecticut", "Rhode Island", "Massachusetts", 
                                                "New Hampshire", "Maine")), y = CPUE, 
                     group = factor(Year),color = factor(Year)))+
  geom_line()+
  facet_wrap(~Species, scales = "free")+
  labs(x = "State") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Avg CPUE over state by Binned 5 years 
ggplot(RangeDF2, aes(x=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                "Connecticut", "Rhode Island", "Massachusetts", 
                                                "New Hampshire", "Maine")), y = Avg5, 
                     group = Bin1,color = Bin1))+
  geom_line()+
  facet_wrap(~Species, scales = "free")+
  labs(x = "State") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Avg CPUE over state by Binned 10 years 
ggplot(RangeDF2, aes(x=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                "Connecticut", "Rhode Island", "Massachusetts", 
                                                "New Hampshire", "Maine")), y = Avg10, 
                     group = Bin2,color = Bin2))+
  geom_line()+
  facet_wrap(~Species, scales = "free")+
  labs(x = "State") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Rolling 10 year average CPUE across states 
ggplot(RangeDF2, aes(x=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                "Connecticut", "Rhode Island", "Massachusetts", 
                                                "New Hampshire", "Maine")), y = RollMean10yr, 
                     group = factor(Year),color = factor(Year)))+
  geom_line()+
  facet_wrap(~Species, scales = "free")+
  labs(x = "State") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Avergae 5 years CPUE BUBBLE PLOT by year across 5 years bins
p <- ggplot(RangeDF2, aes(y=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                     "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                     "Connecticut", "Rhode Island", "Massachusetts", 
                                                     "New Hampshire", "Maine")), 
                          x= Bin1)) +
  geom_point_interactive(aes(size = Avg5, color = SD5, tooltip = Avg5, data_id = Avg5)) +
  #  geom_raster(aes(fill = Avg10))+
  facet_wrap(~Species) +
  scale_size(breaks = seq(0,40,2)) +
  scale_color_continuous(type = "viridis")+
  labs(size = "Average CPUE \n(5 yrs bins)", y = "", x = "Binned Years") +
  theme_bw()

girafe(ggobj = p)




# Raster plot showinf rollings 10 year average 
ggplot(RangeDF2, aes(y=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                "Connecticut", "Rhode Island", "Massachusetts", 
                                                "New Hampshire", "Maine")), 
                     x= Year)) +
  # geom_point(aes(size = RollMean10yr, color = RollMean10yr)) +
  geom_tile(aes(fill = cut(RollMean5yr, c(0, 5, 10, 15, 20, 25, 30,40))))+
  facet_wrap(~Species) +
  scale_size(breaks = seq(0,40,2)) +
  # scale_fill_continuous(type = "viridis")+
  labs(size = "Average CPUE \n(5 yrs bins)", y = "", x = "Year", fill = "Rolling 5 year CPUE mean") +
  theme_bw() +
  scale_fill_brewer(type = "seq", palette = "RdYlGn") +
  #scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(panel.grid = element_blank())






######################## Figures   #############################################
########################     Figure 2  -  MRIP Range Graphs   ##################

##Min and max latitudes over time 
RangeDF4 <- RangeDF2 %>%
  filter(!is.na(CPUE)) %>%
  group_by(Year, Species) %>%
  mutate(xm = Lat * CPUE) %>%
  summarise(xbar = (sum(xm))/sum(CPUE),
            max=max(Lat),
            min = min(Lat),
            range = max(Lat)- min(Lat))

RangeDF4 %>%
  ggplot(aes(x = Year, y = max)) +
  geom_segment(aes(y = min, yend = max, xend = Year)) +
  geom_point(aes(y = min)) +
  geom_point(aes(yend = max))+
  facet_wrap(~Species, nrow = 2, labeller = as_labeller(Spp)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Latitude (°N)", x = "Year") 

# ggsave("CPUE.Range.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
#        path = "C:/Users/RER/Documents/Masters UD/Range Shift")

##Plots for latitudinal expansion and restrictions
RangeDF4 %>%
  ggplot(aes(x = Year, y = range)) +
  geom_line(size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Species, nrow = 2, labeller = as_labeller(Spp)) +
  labs(y = "Range (° Latitude)") +
  stat_smooth(method = "lm", se = F, linetype = "dashed", color = "black", size = .5)

ggsave("Figure2.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")


########################     Figure 2  -  Statistics  ##########################
# Monotonic Trend Test (Mann Kendall) for each Species range
AM2 <- as.ts(RangeDF4 %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(AM2)$range)     # not normal
summary(MannKendall(AM2))                  # NOT significant trend
MannKendall(AM2[,2])

BSB2 <- as.ts(RangeDF4 %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(BSB2)$range)     # not normal
MannKendall(BSB2[,2])                           # significant trend

BF2 <- as.ts(RangeDF4 %>% filter(Species == "Bluefish") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(BF2)$range)     # not normal
MannKendall(BF2[,2])                           # NOT significant trend

C2 <- as.ts(RangeDF4 %>% filter(Species == "Cobia") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(C2)$range)     # not normal
MannKendall(C2[,2])                           # significant trend

FP2 <- as.ts(RangeDF4 %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(FP2)$range)     # not normal
MannKendall(FP2[,2])                           # significant trend

SF2 <- as.ts(RangeDF4 %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(SF2)$range)     # not normal
MannKendall(SF2[,2])                           # NOT significant trend

WF2 <- as.ts(RangeDF4 %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(WF2)$range)     # not normal
MannKendall(WF2[,2])                           # NOT significant trend

YFT2 <- as.ts(RangeDF4 %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(YFT2)$range)     # not normal
MannKendall(YFT2[,2])                           # significant trend

# Quantifying degress of change from 1981-2022
# mx +b
AM2.lm <- RangeDF4 %>% filter(Species == "Atlantic Mackerel")
summary(lm(range ~ Year, AM2.lm))
(( -0.01287 * 2022) + 32.40940) - (( -0.01287 * 1981) + 32.40940)

BSB2.lm <- RangeDF4 %>% filter(Species == "Black Sea Bass")
summary(lm(range ~ Year, BSB2.lm))
(( 0.021759 * 2022) - 32.454136) - (( 0.021759 * 1981) - 32.454136)

BF2.lm <- RangeDF4 %>% filter(Species == "Bluefish")
summary(lm(range ~ Year, BF2.lm))
(( -0.010655 * 2022) + 33.767842) - (( -0.010655 * 1981) + 33.767842)

C2.lm <- RangeDF4 %>% filter(Species == "Cobia")
summary(lm(range ~ Year, C2.lm))
(( 0.06857 * 2022)  -131.16272) - (( 0.06857 * 1981)  -131.16272)

FP2.lm <- RangeDF4 %>% filter(Species == "Florida Pompano")
summary(lm(range ~ Year, FP2))
(( 0.08723 * 2022) - 169.88938 ) - (( 0.08723 * 1981) -169.88938)

SF2.lm <- RangeDF4 %>% filter(Species == "Summer Flounder") 
summary(lm(range ~ Year, SF2.lm))
(( 0.02127 * 2022) -32.52282) - (( 0.02127 * 1981) -32.52282)

WF2.lm <- RangeDF4 %>% filter(Species == "Winter Flounder") 
summary(lm(range ~ Year, WF2.lm))
(( -0.009789 * 2022) + 23.846853) - (( -0.009789 * 1981) + 23.846853)

YFT2.lm <- RangeDF4 %>% filter(Species == "Yellowfin Tuna")
summary(lm(range ~ Year, YFT2.lm))
(( -0.02842 * 2022) + 64.64939) - (( -0.02842 * 1981) + 64.64939)


# Monotonic Trend Test (Mann Kendall) for each Species min and max
AM.Max <- as.ts(RangeDF4 %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, max) %>% arrange(Year))
summary(MannKendall(AM.Max))                  # significant trend
MannKendall(AM.Max[,2])
AM.Min <- as.ts(RangeDF4 %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, min) %>% arrange(Year))
summary(MannKendall(AM.Min[,2]))                  # significant trend

BSB.Max <- as.ts(RangeDF4 %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(BSB.Max[,2])                           # significant trend
BSB.Min <- as.ts(RangeDF4 %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(BSB.Min[,2])                           # significant trend

BF.Max <- as.ts(RangeDF4 %>% filter(Species == "Bluefish") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(BF.Max[,2])                           # significant trend
BF.Min <- as.ts(RangeDF4 %>% filter(Species == "Bluefish") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(BF.Min[,2])                           # significant trend

C.Max <- as.ts(RangeDF4 %>% filter(Species == "Cobia") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(C.Max[,2])                           # significant trend
C.Min <- as.ts(RangeDF4 %>% filter(Species == "Cobia") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(C.Min[,2])                           # significant trend

FP.Max <- as.ts(RangeDF4 %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(FP.Max[,2])                           # significant trend
FP.Min <- as.ts(RangeDF4 %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(FP.Min[,2]) 

SF.Max <- as.ts(RangeDF4 %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(SF.Max[,2])                           # significant trend
SF.Min <- as.ts(RangeDF4 %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(SF.Min[,2])                           # significant trend

WF.Max <- as.ts(RangeDF4 %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(WF.Max[,2])                           # significant trend
WF.Min <- as.ts(RangeDF4 %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(WF.Min[,2])                           # significant trend

YFT.Max <- as.ts(RangeDF4 %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(YFT.Max[,2])                           # significant trend
YFT.Min <- as.ts(RangeDF4 %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(YFT.Min[,2])                           # significant trend

########################     Figure 2  -  SUBSET + stats    ####################

##Min and max latitudes over time 
RangeDF4.SS <- RangeDF.subset %>%
  filter(!is.na(CPUE)) %>%
  group_by(Year, Species) %>%
  mutate(xm = Lat * CPUE) %>%
  summarise(xbar = (sum(xm))/sum(CPUE),
            max=max(Lat),
            min = min(Lat),
            range = max(Lat)- min(Lat))

RangeDF4.SS %>%
  ggplot(aes(x = Year, y = max)) +
  geom_segment(aes(y = min, yend = max, xend = Year)) +
  geom_point(aes(y = min)) +
  geom_point(aes(yend = max))+
  facet_wrap(~Species, nrow = 2, labeller = as_labeller(Spp)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Latitude (°N)", x = "Year") 

# ggsave("CPUE.Range.Subset.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
#        path = "C:/Users/RER/Documents/Masters UD/Range Shift")

##Plots for latitudinal expansion and restrictions
RangeDF4.SS %>%
  ggplot(aes(x = Year, y = range)) +
  geom_line(size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  facet_wrap(~Species, nrow = 2, labeller = as_labeller(Spp)) +
  labs(y = "Range (° Latitude)") +
  stat_smooth(method = "lm", se = F, linetype = "dashed", color = "black", size = .5)

# ggsave("CPUE.Range.Shift.Subset.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
#        path = "C:/Users/RER/Documents/Masters UD/Range Shift")



# Monotonic Trend Test (Mann Kendall) for each Species range
AM2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(AM2.SS)$range)     # not normal
summary(MannKendall(AM2.SS))                  # NOT significant trend
MannKendall(AM2.SS[,2])

BSB2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(BSB2.SS)$range)     # not normal
MannKendall(BSB2.SS[,2])                           # significant trend

BF2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Bluefish") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(BF2.SS)$range)     # not normal
MannKendall(BF2.SS[,2])                           # NOT significant trend

C2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Cobia") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(C2.SS)$range)     # not normal
MannKendall(C2.SS[,2])                           # NOT significant trend

FP2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(FP2.SS)$range)     # not normal
MannKendall(FP2.SS[,2])                           #significant trend

SF2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(SF2.SS)$range)     # not normal
MannKendall(SF2.SS[,2])                           # NOT significant trend

WF2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(WF2.SS)$range)     # not normal
summary(MannKendall(WF2.SS[,2]))                           # NOT significant trend

YFT2.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, range) %>% arrange(Year))
shapiro.test(as.data.frame(YFT2.SS)$range)     # not normal
MannKendall(YFT2.SS[,2])                           #NOT  significant trend

# Quantifying degress of change from 1981-2022
# mx +b

AM2.lm.SS <- RangeDF4.SS %>% filter(Species == "Atlantic Mackerel")
BSB2.lm.SS <- RangeDF4.SS %>% filter(Species == "Black Sea Bass")
BF2.lm.SS <- RangeDF4.SS %>% filter(Species == "Bluefish")
C2.lm.SS <- RangeDF4.SS %>% filter(Species == "Cobia")
FP2.lm.SS <- RangeDF4.SS %>% filter(Species == "Florida Pompano")
SF2.lm.SS <- RangeDF4.SS %>% filter(Species == "Summer Flounder") 
WF2.lm.SS <- RangeDF4.SS %>% filter(Species == "Winter Flounder") 
YFT2.lm.SS <- RangeDF4.SS %>% filter(Species == "Yellowfin Tuna")

summary(lm(range ~ Year, AM2.lm.SS))
(( -0.01287 * 2022) + 32.40940) - (( -0.01287 * 1981) + 32.40940)

summary(lm(range ~ Year, BSB2.lm.SS))
(( 0.021759 * 2022) - 35.816502) - (( 0.021759 * 1981) - 35.816502)

summary(lm(range ~ Year, BF2.lm.SS))
(( -0.010655 * 2022) + 30.405476) - (( -0.010655 * 1981) + 30.405476)

summary(lm(range ~ Year, C2.lm.SS))
(( 0.02780 * 2022)  -52.43003) - (( 0.02780 * 1981)  -52.43003)

summary(lm(range ~ Year, FP2.lm.SS))
(( 0.07066 * 2022) - 139.55594 ) - (( 0.07066 * 1981) -139.55594)

summary(lm(range ~ Year, SF2.lm.SS))
(( -0.022161 * 2022) + 52.026123  ) - (( -0.022161 * 1981) + 52.026123  )

summary(lm(range ~ Year, WF2.lm.SS))
(( -0.009789 * 2022) + 23.846853) - (( -0.009789 * 1981) + 23.846853)

summary(lm(range ~ Year, YFT2.lm.SS))
(( -0.02787     * 2022)  -49.17820) - (( 0.02787     * 1981) -49.17820)


# Monotonic Trend Test (Mann Kendall) for each Species min and max
AM.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, max) %>% arrange(Year))
summary(MannKendall(AM.Max.SS))                  # significant trend
MannKendall(AM.Max.SS[,2])
AM.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, min) %>% arrange(Year))
summary(MannKendall(AM.Min.SS[,2]))                  # NOT significant trend

BSB.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(BSB.Max.SS[,2])                           # significant trend
BSB.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(BSB.Min.SS[,2])                           #NOT  significant trend

BF.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Bluefish") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(BF.Max.SS[,2])                           # NOT significant trend
BF.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Bluefish") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(BF.Min.SS[,2])                           # NOT significant trend

C.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Cobia") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(C.Max.SS[,2])                           # NOT significant trend
C.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Cobia") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(C.Min.SS[,2])                           # NOT significant trend

FP.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(FP.Max.SS[,2])                           # significant trend
FP.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(FP.Min.SS[,2])                          # NOT sig

SF.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(SF.Max.SS[,2])                           # NOT significant trend
SF.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(SF.Min.SS[,2])                           # NOT significant trend

WF.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(WF.Max.SS[,2])                           # NOT significant trend
WF.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(WF.Min.SS[,2])                           # NOT significant trend

YFT.Max.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, max) %>% arrange(Year))
MannKendall(YFT.Max.SS[,2])                           # NOT significant trend
YFT.Min.SS <- as.ts(RangeDF4.SS %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, min) %>% arrange(Year))
MannKendall(YFT.Min.SS[,2])                           # significant trend



########################     Figure 3  -  MRIP Center of Distribution  #########
# Weighted average of CPUE
# Center of Gravity (latitude biomass weighted average) of CPUE per species
RangeDF3 <- RangeDF2 %>%
  group_by(Year, Species) %>%
  filter(!is.na(CPUE)) %>%
  mutate(xm = Lat * CPUE) %>%
  summarise(xbar = (sum(xm))/sum(CPUE))

# Figure 4
ggplot(RangeDF3, aes(x=Year, y = xbar)) +
  geom_line(size = 1) +
  facet_wrap(~Species, scales ="free_y", ncol = 4, labeller = as_labeller(Spp)) +
  labs(y = "Latitude (°N)", x = "Year") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_smooth(method = "lm", se = F, linetype = "dashed", color = "black", size = .5) +
  scale_y_continuous(labels = formatter(nsmall = 1))


# ggsave("CPUE.COG.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
#        path = "C:/Users/RER/Documents/Masters UD/Range Shift")

ggplot(RangeDF4, aes(x=Year, y = xbar), color ="black") +
  geom_ribbon(aes(x = Year, ymin = min, ymax = xbar), fill = "lightgrey") +
  geom_line(size = 2) +
  geom_line(aes(x= Year, y = min), color ="grey") +
  geom_line(aes(x= Year, y = max), color ="grey") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Latitude (°N)", x = "Year") +
  geom_ribbon(aes(x = Year, ymin = xbar, ymax = max), fill = "lightgrey") +
  stat_smooth(method = "lm", color = "white", linetype = "dashed", level = 100) +
  facet_wrap(~Species, nrow = 2, labeller = as_labeller(Spp)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.5), breaks = c(32.5, 35, 37.5, 40, 42.5, 45)) +
  scale_x_continuous(limits = c(1980,2022), expand = expansion(mult = c(0.05, .01))) #.1 first


ggsave("Figure3.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")

########################     Figure 3  -  Statistics  ##########################

# Monotonic Trend Test (Mann Kendall) for each Species
AM <- as.ts(RangeDF3 %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(AM)$xbar)      #normal
MannKendall(AM[,2])                           # YES significant trend

BSB <- as.ts(RangeDF3 %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(BSB)$xbar)      #normal
MannKendall(BSB[,2])                           # YES significant trend

BF <- as.ts(RangeDF3 %>% filter(Species == "Bluefish") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(BF)$xbar)      #normal
MannKendall(BF[,2])                           # significant trend

C <- as.ts(RangeDF3 %>% filter(Species == "Cobia") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(C)$xbar)      # non normal
MannKendall(C[,2])                           #NOT  significant trend

FP <- as.ts(RangeDF3 %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(FP)$xbar)      # normal
MannKendall(FP[,2])                           # significant trend

SF <- as.ts(RangeDF3 %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(SF)$xbar)      #  normal
MannKendall(SF[,2])                           # NOT significant trend

WF <- as.ts(RangeDF3 %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(WF)$xbar)      #  normal
MannKendall(WF[,2])                           #NOT significant trend

YFT <- as.ts(RangeDF3 %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(YFT)$xbar)      #  normal
MannKendall(YFT[,2])                           #NOT significant trend

# Quantifying degress of change from 1981-2022
# mx +b
AM.lm <- RangeDF3 %>% filter(Species == "Atlantic Mackerel")
BSB.lm <- RangeDF3 %>% filter(Species == "Black Sea Bass")
BF.lm <- RangeDF3 %>% filter(Species == "Bluefish")
C.lm <- RangeDF3 %>% filter(Species == "Cobia")
FP.lm <- RangeDF3 %>% filter(Species == "Florida Pompano")
SF.lm <- RangeDF3 %>% filter(Species == "Summer Flounder") 
WF.lm <- RangeDF3 %>% filter(Species == "Winter Flounder") 
YFT.lm <- RangeDF3 %>% filter(Species == "Yellowfin Tuna")

summary(lm(xbar ~ Year, AM.lm))
(( 0.06191 * 2022)  -82.93557) - (( 0.06191 * 1981)  -82.93557)

summary(lm(xbar ~ Year, BSB.lm))
(( 0.042019 * 2022) - 46.674399) - (( 0.042019 * 1981) - 46.674399)

summary(lm(xbar ~ Year, BF.lm))
(( -0.051753 * 2022) + 141.256797) - (( -0.051753 * 1981) + 141.256797)

summary(lm(xbar ~ Year, C.lm))
(( -0.007376 * 2022) + 49.310395) - (( -0.007376 * 1981)  + 49.310395)

summary(lm(xbar ~ Year, FP.lm))
(( 0.022731 * 2022) - 11.424745 ) - (( 0.022731 * 1981) -11.424745)

summary(lm(xbar ~ Year, SF.lm))
(( -0.005047 * 2022) + 49.352037) - (( -0.005047 * 1981) + 49.352037)

summary(lm(xbar ~ Year, WF.lm))
(( 0.003819 * 2022) + 33.859275) - (( 0.003819 * 1981) + 33.859275)

summary(lm(xbar ~ Year, YFT.lm))
(( 0.009225 * 2022) + 19.993513) - (( 0.009225 * 1981) + 19.993513)

########################     Figure 3  -  SUBSET + Stats  ######################


# Weighted average of CPUE
# Center of Gravity (latitude biomass weighted average) of CPUE per species
RangeDF3.SS <- RangeDF.subset %>%
  group_by(Year, Species) %>%
  filter(!is.na(CPUE)) %>%
  mutate(xm = Lat * CPUE) %>%
  summarise(xbar = (sum(xm))/sum(CPUE))

# Figure 2
ggplot(RangeDF3.SS, aes(x=Year, y = xbar)) +
  geom_line(size = 1) +
  facet_wrap(~Species, scales ="free_y", ncol = 4, labeller = as_labeller(Spp)) +
  labs(y = "Latitude (°N)", x = "Year") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_smooth(method = "lm", se = F, linetype = "dashed", color = "black", size = .5) +
  scale_y_continuous(labels = formatter(nsmall = 1))


# ggsave("CPUE.COG.Subset.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
#        path = "C:/Users/RER/Documents/Masters UD/Range Shift")

ggplot(RangeDF4.SS, aes(x=Year, y = xbar), color ="black") +
  geom_ribbon(aes(x = Year, ymin = min, ymax = xbar), fill = "lightgrey") +
  geom_line(size = 2) +
  geom_line(aes(x= Year, y = min), color ="grey") +
  geom_line(aes(x= Year, y = max), color ="grey") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_ribbon(aes(x = Year, ymin = xbar, ymax = max), fill = "lightgrey") +
  stat_smooth(method = "lm", color = "white", linetype = "dashed", level = 100) +
  facet_wrap(~Species, nrow = 2, labeller = as_labeller(Spp)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.5)) +
  scale_x_continuous(limits = c(1980,2022), expand = expansion(mult = c(.2, 0))) 


ggsave("CPUE.Range.Dist.MRIP.Subset.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift")

# Monotonic Trend Test (Mann Kendall) for each Species
AM.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Atlantic Mackerel") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(AM.SS)$xbar)      #normal
MannKendall(AM.SS[,2])                           # YES significant trend

BSB.SS <- as.ts(RangeDF3 %>% filter(Species == "Black Sea Bass") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(BSB.SS)$xbar)      #normal
MannKendall(BSB.SS[,2])                           # YES significant trend

BF.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Bluefish") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(BF.SS)$xbar)      #normal
MannKendall(BF.SS[,2])                           # significant trend

C.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Cobia") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(C.SS)$xbar)      # non normal
MannKendall(C.SS[,2])                           #NOT  significant trend

FP.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Florida Pompano") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(FP.SS)$xbar)      # NOT normal
MannKendall(FP.SS[,2])                           # significant trend

SF.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Summer Flounder") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(SF.SS)$xbar)      # NOT  normal
MannKendall(SF.SS[,2])                           # YES significant trend

WF.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Winter Flounder") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(WF.SS)$xbar)      #  normal
MannKendall(WF.SS[,2])                           #NOT significant trend

YFT.SS <- as.ts(RangeDF3.SS %>% filter(Species == "Yellowfin Tuna") %>% dplyr::select(Year, xbar) %>% arrange(Year))
shapiro.test(as.data.frame(YFT.SS)$xbar)      #  normal
MannKendall(YFT.SS[,2])                           #NOT significant trend

# Quantifying degress of change from 1981-2022
# mx +b
AM.lm.SS <- RangeDF3.SS %>% filter(Species == "Atlantic Mackerel")
BSB.lm.SS <- RangeDF3.SS %>% filter(Species == "Black Sea Bass")
BF.lm.SS <- RangeDF3.SS %>% filter(Species == "Bluefish")
C.lm.SS <- RangeDF3.SS %>% filter(Species == "Cobia")
FP.lm.SS <- RangeDF3.SS %>% filter(Species == "Florida Pompano")
SF.lm.SS <- RangeDF3.SS %>% filter(Species == "Summer Flounder") 
WF.lm.SS <- RangeDF3.SS %>% filter(Species == "Winter Flounder") 
YFT.lm.SS <- RangeDF3.SS %>% filter(Species == "Yellowfin Tuna")

summary(lm(xbar ~ Year, AM.lm.SS))
(( 0.06191 * 2022)  -82.93557) - (( 0.06191 * 1981)  -82.93557)

summary(lm(xbar ~ Year, BSB.lm.SS))
(( 0.033610 * 2022) - 28.427298) - (( 0.033610 * 1981) - 28.427298)

summary(lm(xbar ~ Year, BF.lm.SS))
(( -0.037620 * 2022) + 114.363571) - (( -0.037620 * 1981) + 114.363571)

summary(lm(xbar ~ Year, C.lm.SS))
(( 0.01525 * 2022) + 6.14286   ) - (( 0.01525 * 1981)  + 6.14286   )

summary(lm(xbar ~ Year, FP.lm.SS))
(( 0.035147 * 2022) - 34.829934 ) - (( 0.035147 * 1981) -34.829934)

summary(lm(xbar ~ Year, SF.lm.SS))
(( 0.007634 * 2022) + 24.385790) - (( 0.007634 * 1981) + 24.385790)

summary(lm(xbar ~ Year, WF.lm.SS))
(( 0.003819 * 2022) + 33.859275) - (( 0.003819 * 1981) + 33.859275)

summary(lm(xbar ~ Year, YFT.lm.SS))
(( -0.010701 * 2022) + 60.187710) - (( -0.010701 * 1981) + 60.187710)


########################     Figure 4  -  NEFSC Range graphs ###################

##Plots for latitudinal expansion and restrictions using center lat
RANGE %>%
  filter(!com2 == "Cobia") %>%
  ggplot(aes(x = year, y = range)) +
  geom_line(size = 1) +
  facet_wrap(~com2, labeller = as_labeller(Spp)) +
  labs(y = "Range (°Latitude)", x = "Year") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_smooth(method = "lm", se = F, linetype = "dashed", color = "black", size = .5) 

ggsave("Figure4.png", plot = last_plot(), width = 8,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")

##Min and max latitudes over time using center lat
RANGE %>%
  filter(!com2 == "Cobia") %>%
  ggplot(aes(x = year, y = max)) +
  geom_segment(aes(y = min, yend = max, xend = year)) +
  geom_point(aes(y = min)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_point(aes(yend = max)) +
  facet_wrap(~com2, labeller = as_labeller(Spp)) +
  labs(y = "Latitude (°N)", x = "Year") 

# ggsave("Bigelow.MinMax.png", plot = last_plot(), width = 8,height = 5, units = "in", dpi = 500,
#        path = "C:/Users/RER/Documents/Masters UD/Range Shift")

########################     Figure 4  -  Statistics #########
#Normality for BTS Range
AM_species <- RANGE[RANGE$com == 'Atlantic mackerel ', ]
shapiro.test(AM_species$range) ##non normal

BSB_species <- RANGE[RANGE$com == 'black sea bass', ]
shapiro.test(BSB_species$range) ##non normal

BF_species <- RANGE[RANGE$com == 'bluefish  NA', ]
shapiro.test(BF_species$range) ##non normal

SF_species <- RANGE[RANGE$com == 'summer flounder ', ]
shapiro.test(SF_species$range) ##non normal

WF_species <- RANGE[RANGE$com == 'winter flounder ', ]
shapiro.test(WF_species$range) ##non normal

#Normality for BTS Distribution
AM_species.d <- DIST[DIST$com == 'Atlantic mackerel ', ]
shapiro.test(AM_species.d$xbar) ##normal

BSB_species.d <- DIST[DIST$com == 'black sea bass', ]
shapiro.test(BSB_species.d$xbar) ##normal

BF_species.d <- DIST[DIST$com == 'bluefish  NA', ]
shapiro.test(BF_species.d$xbar) ##non normal

SF_species.d <- DIST[DIST$com == 'summer flounder ', ]
shapiro.test(SF_species.d$xbar) ##normal

WF_species.d <- DIST[DIST$com == 'winter flounder ', ]
shapiro.test(WF_species.d$xbar) ##non normal


##Trend test Range 
#Atlantic Mackerel
AM.Range <- as.ts(RANGE %>% filter(com == "Atlantic mackerel ") %>% dplyr::select(year, range) %>% arrange(year))
MannKendall(AM.Range[,2])  ## not significant

##Black Sea Bass
BSB.Range <- as.ts(RANGE %>% filter(com == "black sea bass") %>% dplyr::select(year, range) %>% arrange(year))
MannKendall(BSB.Range[,2]) ## not significant

##Bluefish
BF.Range <- as.ts(RANGE %>% filter(com == "bluefish  NA") %>% dplyr::select(year, range) %>% arrange(year))
MannKendall(BF.Range[,2]) ## not significant

##Cobia
CB.Range <- as.ts(RANGE %>% filter(com == "cobia  NA") %>% dplyr::select(year, range) %>% arrange(year))
MannKendall(CB.Range[,2]) ## not significant

##Summer Flounder
SF.Range <- as.ts(RANGE %>% filter(com == "summer flounder ") %>% dplyr::select(year, range) %>% arrange(year))
MannKendall(SF.Range) ## not significant

##Winter Flounder
WF.Range <- as.ts(RANGE %>% filter(com == "winter flounder ") %>% dplyr::select(year, range) %>% arrange(year))
MannKendall(WF.Range[,2]) ##significant

#### Mann Kendall on Max and Min #
#Atlantic Mackerel
AM.Max <- as.ts(RANGE %>% filter(com == "Atlantic mackerel ") %>% dplyr::select(year, max) %>% arrange(year))
MannKendall(AM.Max[,2]) #significant

AM.Min <- as.ts(RANGE %>% filter(com == "Atlantic mackerel ") %>% dplyr::select(year, min) %>% arrange(year))
MannKendall(AM.Min[,2]) # significant

#Black Sea Bass
BSB.Max <- as.ts(RANGE %>% filter(com == "black sea bass") %>% dplyr::select(year, max) %>% arrange(year))
MannKendall(BSB.Max[,2]) #significant

BSB.Min <- as.ts(RANGE %>% filter(com == "black sea bass") %>% dplyr::select(year, min) %>% arrange(year))
MannKendall(BSB.Min[,2]) # significant

#Bluefish
BF.Max <- as.ts(RANGE %>% filter(com == "bluefish  NA") %>% dplyr::select(year, max) %>% arrange(year))
MannKendall(BF.Max[,2]) #significant

BF.Min <- as.ts(RANGE %>% filter(com == "bluefish  NA") %>% dplyr::select(year, min) %>% arrange(year))
MannKendall(BF.Min[,2]) # significant

#Cobia
CB.Max <- as.ts(RANGE %>% filter(com == "cobia  NA") %>% dplyr::select(year, max) %>% arrange(year))
MannKendall(CB.Max[,2]) #significant

CB.Min <- as.ts(RANGE %>% filter(com == "cobia  NA") %>% dplyr::select(year, min) %>% arrange(year))
MannKendall(CB.Min[,2]) # significant

#Summer Flounder
SF.Max <- as.ts(RANGE %>% filter(com == "summer flounder ") %>% dplyr::select(year, max) %>% arrange(year))
MannKendall(SF.Max[,2]) # not significant

SF.Min <- as.ts(RANGE %>% filter(com == "summer flounder ") %>% dplyr::select(year, min) %>% arrange(year))
MannKendall(SF.Min[,2]) # significant

#Winter Flounder
WF.Max <- as.ts(RANGE %>% filter(com == "winter flounder ") %>% dplyr::select(year, max) %>% arrange(year))
MannKendall(WF.Max[,2]) # not significant

WF.Min <- as.ts(RANGE %>% filter(com == "winter flounder ") %>% dplyr::select(year, min) %>% arrange(year))
MannKendall(WF.Min[,2]) # significant

########################     Figure 5  -  NEFSC Center of Distribution  ########

DIST2 %>%
  filter(!com2 == "Cobia") %>%
  ggplot(aes(x=year, y = xbar), color ="black") +
  geom_ribbon(aes(x = year, ymin = min, ymax = xbar), fill = "lightgrey") +
  geom_line(size = 2) +
  geom_line(aes(x= year, y = min), color ="grey") +
  geom_line(aes(x= year, y = max), color ="grey") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Latitude (°N)", x = "Year") +
  geom_ribbon(aes(x = year, ymin = xbar, ymax = max), fill = "lightgrey") +
  stat_smooth(method = "lm", color = "white", linetype = "dashed", level = 100) +
  facet_wrap(~com2, nrow = 2, labeller = as_labeller(Spp))


  ggplot(aes(x= year, y = xbar)) +
  geom_line(size = 1) +
  facet_wrap(~com, scales ="free_y", ncol = 3, labeller = as_labeller(Spp)) +
  labs(y = "Latitude (°N)", x = "Year") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_smooth(method = "lm", se = F, linetype = "dashed", color = "black", size = .5) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

ggsave("Figure5.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")

########################     Figure 5  -  Statistics ################

AM.NEFSC <- DIST %>%
  filter(com2 == "Atlantic Mackerel",
         Data.Type == "Bigelow")
summary(lm(xbar ~ year, AM.NEFSC))
# ((0.05795*2022) - 76.21625) - ((0.05795*1981) - 76.21625)
# 2.37595 degree change

BSB.NEFSC <- DIST %>%
  filter(com2 == "Black Sea Bass",
         Data.Type == "Bigelow")
summary(lm(xbar ~ year, BSB.NEFSC))
# ((0.065344*2022) - 92.692144) - ((0.065344*1981) - 92.692144)
# 2.679104 degree change

BF.NEFSC <- DIST %>%
  filter(com2 == "Bluefish",
         Data.Type == "Bigelow")
summary(lm(xbar ~ year, BF.NEFSC))
# ((0.050792*2022) - 65.904254) - ((0.050792*1981) - 65.904254)
# 2.082472 degree change

CB.NEFSC <- DIST %>%
  filter(com2 == "Cobia",
         Data.Type == "Bigelow")
summary(lm(xbar ~ year, CB.NEFSC))
# ((0.021209*2022) - 7.865950) - ((0.021209*1981) - 7.865950)
# 0.869569 degree change

SF.NEFSC <- DIST %>%
  filter(com2 == "Summer Flounder",
         Data.Type == "Bigelow")
summary(lm(xbar ~ year, SF.NEFSC))
# ((0.035198*2022) - 32.107763) - ((0.035198*1981) - 32.107763)
# 1.443118 degree change

WF.NEFSC <- DIST %>%
  filter(com2 == "Winter Flounder",
         Data.Type == "Bigelow")
summary(lm(xbar ~ year, WF.NEFSC))
# ((0.008702*2022) + 23.872356) - ((0.008702*1981) + 23.872356)
# 0.356782 degree change 

#Atlantic Mackerel
AM.Dist <- as.ts(DIST %>% filter(com2 == "Atlantic Mackerel") %>% dplyr::select(year, xbar) %>% arrange(year))
MannKendall(AM.Dist[,3])  ##significant 

##Black Sea Bass
BSB.Dist <- as.ts(DIST %>% filter(com2 == "Black Sea Bass") %>% dplyr::select(year, xbar) %>% arrange(year))
MannKendall(BSB.Dist[,3]) ##significant

##Bluefish
BF.Dist <- as.ts(DIST %>% filter(com2 == "Bluefish") %>% dplyr::select(year, xbar) %>% arrange(year))
MannKendall(BF.Dist[,3]) ##significant

##Cobia
CB.Dist <- as.ts(DIST %>% filter(com2 == "Cobia") %>% dplyr::select(year, xbar) %>% arrange(year))
MannKendall(CB.Dist[,3]) ##significant

##Summer Flounder
SF.Dist <- as.ts(DIST %>% filter(com2 == "Summer Flounder") %>% dplyr::select(year, xbar) %>% arrange(year))
MannKendall(SF.Dist[,3]) ##significant

##Winter Flounder
WF.Dist <- as.ts(DIST %>% filter(com2 == "Winter Flounder") %>% dplyr::select(year, xbar) %>% arrange(year))
MannKendall(WF.Dist[,3]) ##significant

########################     Figure 6  -  SST Anomalies Time Series  ###########

# # 1980-2022 temperature anomaly average along east coast
# ggplot(sst.trend2, aes(x = Year, y = anom2), color = "red") +
#   geom_line(size =1) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   labs(y = "Temeprature Anomoly (°C)", x = "Year") 


SST.anom <- SST %>%
  ungroup() %>%
  mutate(baseline.avg = mean(sst, na.rm = T)) %>%
  group_by(Year) %>%
  reframe(yearly.avg = mean(sst, na.rm = T),
          anom = yearly.avg - baseline.avg) %>%
  distinct() %>%
  na.omit()

# 1980-2022 temperature anomaly average along east coast
ggplot(SST.anom, aes(x = Year, y = anom)) +
  geom_line(size =1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(y = "Temeprature Anomoly (°C)", x = "Year") +
  geom_line(y = 0, linetype = "dashed")

SST.anom.2 <- as.ts(SST.anom  %>% dplyr::select(Year, anom) %>% distinct() %>% arrange(Year))
MannKendall(SST.anom.2)

ggsave("Figure6.png", plot = last_plot(), width = 5,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")
########################  Supplemental       ####################################
########################     Figure S1 -  SST and CPUE ccf  ###########

Cobia.SST <- full_join(SST.anom, Cobia) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

C.ccf <- ccf(as.ts(Cobia.SST$anom), as.ts(Cobia.SST$xbar), plot = F)  # NOT SIGNIF
plot(C.ccf[C.ccf$lag[1:15],])
ggsave("CobiaCCF.png", plot = last_plot(), width = 3,height = 3, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift")


BF.SST <- full_join(SST.anom, Bluefish) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

BF.ccf <- ccf(as.ts(BF.SST$anom), as.ts(BF.SST$xbar), plot = F)       # SIGNIF (0-15 years)
plot(BF.ccf[BF.ccf$lag[1:15],])

BSB.SST <- full_join(SST.anom, BSB) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

BSB.ccf <- ccf(as.ts(BSB.SST$anom), as.ts(BSB.SST$xbar), plot = F)     # SIGNIF (0-10 years)
plot(BSB.ccf[BSB.ccf$lag[1:15],])


AM.SST <- full_join(SST.anom, AM) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

AM.ccf <- ccf(as.ts(AM.SST$anom), as.ts(AM.SST$xbar), plot = F)      # SIGNIF (-1 - 12 years)
plot(AM.ccf[AM.ccf$lag[1:15],])


WF.SST <- full_join(SST.anom, WF) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

WF.ccf <- ccf(as.ts(WF.SST$anom), as.ts(WF.SST$xbar), plot = F)     # NOT SIGNIF
plot(WF.ccf[WF.ccf$lag[1:15],])



FLP.SST <- full_join(SST.anom, FLP) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

FLP.ccf <- ccf(as.ts(FLP.SST$anom), as.ts(FLP.SST$xbar), plot = F)     # SIGNIF (0 year)
plot(FLP.ccf[FLP.ccf$lag[1:15],])


SF.SST <- full_join(SST.anom, SF) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

SF.ccf <- ccf(as.ts(SF.SST$anom), as.ts(SF.SST$xbar), plot = F)      # SIGNIF (3 year)
plot(SF.ccf[SF.ccf$lag[1:15],])



YFT.SST <- full_join(SST.anom, YFT) %>%
  mutate(lag = lag(Year, n = 1),
         timediff = Year - lag,
         ABORT = timediff > 1) %>%
  # filter(ABORT == T) %>%
  dplyr::select(anom, xbar)

YFT.ccf <- ccf(as.ts(YFT.SST$anom), as.ts(YFT.SST$xbar), plot = F)       # SIGNIF (0-1 year)
plot(YFT.ccf[YFT.ccf$lag[1:15],])

########################     Figure S2 -  PSE Values      ##################################
RangeDF2 %>%
  ungroup() %>%
  group_by(State) %>%
  summarize(minPSE = min(PSE, na.rm = T),
            maxPSE = max(PSE, na.rm = T),
            avgPSE = mean(PSE, na.rm = T))

RangeDF2 %>%
  ungroup() %>%
  group_by(State) %>%
  ggplot(aes(x = Year, y = PSE)) +
  geom_point()+
  facet_wrap(~State)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(y = "Percent Standard Error")


ggsave("PSE_State.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")


RangeDF2 %>%
  ungroup() %>%
  group_by(Species) %>%
  summarize(minPSE = min(PSE, na.rm = T),
            maxPSE = max(PSE, na.rm = T),
            avgPSE = mean(PSE, na.rm = T))

RangeDF2 %>%
  ungroup() %>%
  group_by(Species) %>%
  ggplot(aes(x = Year, y = PSE)) +
  geom_point()+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(y = "Percent Standard Error")


ggsave("PSE_Species.png", plot = last_plot(), width = 10,height = 5, units = "in", dpi = 500,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")

########################     Figure S3 -  Raw Data      ########################


# Create a custom transformation function for the size scale
custom_trans <- scales::trans_new(
  name = "custom",
  transform = function(x) ifelse(x > 200, 200 + (x - 200) / 3, x),  # Flatten the scale after 100
  inverse = function(x) ifelse(x > 200, 200 + (x - 200) * 3, x)  # Inverse transformation
)

# # Create a custom color gradient function
# custom_color_gradient <- function(value) {
#   colors <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(200)
#   color_fun <- colorRamp(colors)
#   return(color_fun((value - min(value)) / (max(value) - min(value))))
# }

# Cap color values at 200
RangeDF2$color_value <- pmin(RangeDF2$CPUE, 150)

ggplot(RangeDF2, aes(y=factor(State, levels = c("Georgia","South Carolina", "North Carolina",
                                                "Virginia", "Maryland", "Delaware", "New Jersey","New York",
                                                "Connecticut", "Rhode Island", "Massachusetts", 
                                                "New Hampshire", "Maine")), 
                     x= Year)) +
  geom_point(aes(size = CPUE, color = color_value))+
  facet_wrap(~Species) +
  labs(size = "CPUE", y = "", x = "Year", fill = "", color = "CPUE (capped at 150)") +
  theme_bw() +
  #scale_color_brewer(type = "seq", palette = "RdYlGn") +
  #scale_x_discrete(expand=c(0,0)) + 
 # scale_y_discrete(expand=c(0,0)) +
 # theme(panel.grid = element_blank())+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())+
  scale_size_continuous(trans = custom_trans, 
                        range = c(1, 10),  # Adjust the range as needed
                        breaks = c(0, 25, 50, 75, 100, 150),
                        labels = c("0", "25", "50", "75", "100", ">150")) +
  scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), 
                        values = scales::rescale(c(0, 25, 50, 75, 100, 150)))   # Adjust color values to cap at 200
  # scale_color_distiller(palette = "Spectral") +
  # scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral")), 
  #                       values = scales::rescale(c(0, 25, 75, 100, 400))) # Apply Spectral color scale
  # scale_size(breaks = seq(0,400,50)) +
#  scale_size(breaks = c(0,25, 50,75, 100,400)) +
  #scale_color_gradient(low = "#5A5A5B", high = "red")
#  scale_color_viridis_c()
  #scale_color_continuous(type = "viridis")
  #scale_colour_gradientn(
    colours = c('darkgreen', 'forestgreen', 'darkseagreen3', 'darkseagreen2',
                           'indianred1', 'indianred2', 'indianred3', 'darkred')


ggsave("SupplFigRawData.png", plot = last_plot(), width = 10,height = 8, units = "in", dpi = 750,
       path = "C:/Users/RER/Documents/Masters UD/Range Shift/FinalFigs")

########################     Percentage of Positive Tows #######################


#### Winter Flounder
WF.tow <- unfiltered.catch %>%
  mutate(presence = case_when(com == ("winter flounder ") ~ 1, 
                              .default = 0)) %>%
  group_by(CRUISE, STRATUM) %>%
  summarise(sum = sum(presence)) %>%
  mutate(pos.tow = case_when((sum > 0) ~ 1,
                             .default =  0)) 
WF.tow %>%
  ungroup %>%
  summarise(length = length(pos.tow),
            sum.tow = sum(pos.tow),
            prop.tow = sum.tow/length * 100)
##45% of all tows

#### Atlantic Mackerel
AM.tow <- unfiltered.catch %>%
  mutate(presence = case_when(com == ("Atlantic mackerel ") ~ 1, 
                              .default = 0)) %>%
  group_by(CRUISE, STRATUM) %>%
  summarise(sum = sum(presence)) %>%
  mutate(pos.tow = case_when((sum > 0) ~ 1,
                             .default =  0)) 
AM.tow %>%
  ungroup %>%
  summarise(length = length(pos.tow),
            sum.tow = sum(pos.tow),
            prop.tow = sum.tow/length * 100)
####31% of all tows

#### Bluefish
BF.tow <- unfiltered.catch %>%
  mutate(presence = case_when(com == ("bluefish  NA") ~ 1, 
                              .default = 0)) %>%
  group_by(CRUISE, STRATUM) %>%
  summarise(sum = sum(presence)) %>%
  mutate(pos.tow = case_when((sum > 0) ~ 1,
                             .default =  0)) 
BF.tow %>%
  ungroup %>%
  summarise(length = length(pos.tow),
            sum.tow = sum(pos.tow),
            prop.tow = sum.tow/length * 100)
####5% of all tows

#### Black Sea Bass
BSB.tow <- unfiltered.catch %>%
  mutate(presence = case_when(com == ("black sea bass") ~ 1, 
                              .default = 0)) %>%
  group_by(CRUISE, STRATUM) %>%
  summarise(sum = sum(presence)) %>%
  mutate(pos.tow = case_when((sum > 0) ~ 1,
                             .default =  0)) 
BSB.tow %>%
  ungroup %>%
  summarise(length = length(pos.tow),
            sum.tow = sum(pos.tow),
            prop.tow = sum.tow/length * 100)
####18% of all tows

#### Summer Flounder
SF.tow <- unfiltered.catch %>%
  mutate(presence = case_when(com == ("summer flounder ") ~ 1, 
                              .default = 0)) %>%
  group_by(CRUISE, STRATUM) %>%
  summarise(sum = sum(presence)) %>%
  mutate(pos.tow = case_when((sum > 0) ~ 1,
                             .default =  0)) 
SF.tow %>%
  ungroup %>%
  summarise(length = length(pos.tow),
            sum.tow = sum(pos.tow),
            prop.tow = sum.tow/length * 100)
####37% of all tows

#### Cobia

C.tow <- unfiltered.catch %>%
  mutate(presence = case_when(com == ("cobia  NA") ~ 1, 
                              .default = 0)) %>%
  group_by(CRUISE, STRATUM) %>%
  summarise(sum = sum(presence)) %>%
  mutate(pos.tow = case_when((sum > 0) ~ 1,
                             .default =  0)) 
C.tow %>%
  ungroup %>%
  summarise(length = length(pos.tow),
            sum.tow = sum(pos.tow),
            prop.tow = sum.tow/length * 100)
####0.5% of all tows
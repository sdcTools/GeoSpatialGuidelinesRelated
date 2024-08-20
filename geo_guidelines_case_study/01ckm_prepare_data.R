##### R-Code for Case Study: Population grids with the Cell Key Method ####### #
#
# 01 - Prepare artificial micro data for the CKM case study
#
# contact: Martin Möhler, martin.moehler@destatis.de
# 2024-08-12
############################################################################## #

library(data.table)
library(dplyr)

# customize to local setting
data_folder <- paste0(getwd(), "/geo_guidelines_case_study")


# ----- get base data from web -----

# Zensus 2011 demographic data by 1ha grid cell
web_dem <- paste0("https://www.zensus2022.de/static/DE/gitterzellen/",
                  "csv_Demographie_100_Meter-Gitter.zip")
# download
download.file(web_dem, paste0(data_folder, "/csv_Demographie_100_Meter-Gitter.zip"),
              mode = "wb") # ~ 315 MB
# unzip
unzip(paste0(data_folder, "/csv_Demographie_100_Meter-Gitter.zip"), exdir = data_folder)


# ----- read in data -----

dem <- fread(file.path(data_folder, "Bevoelkerung100M.csv"))

# get coordinates of grid centroids
dem$x_mp_100m <- as.integer(substr(dem$Gitter_ID_100m_neu, 24, 30)) + 50
dem$y_mp_100m <- as.integer(substr(dem$Gitter_ID_100m_neu, 16, 22)) + 50

# subset to LAEA100kmN31E41
dem <- dem[dem$x_mp_100m > 4100000 & dem$x_mp_100m < 4200000 &
             dem$y_mp_100m > 3100000 & dem$y_mp_100m < 3200000, ]


# ----- create artificial micro data -----

# data for population total 
pop <- dem[dem$Merkmal == " INSGESAMT", ]

# aggregate to micro data
pop <- pop[pop$Anzahl != 0, 
           .(1:Anzahl, x_mp_100m, y_mp_100m), by = .(Gitter_ID_100m)] %>%
  dplyr::select(-V1)

# counts of population aged 50+
dem50p <- dem[dem$Merkmal == "ALTER_KURZ" &
                dem$Auspraegung_Code >= 4, c("Gitter_ID_100m", "Anzahl")] %>%
  group_by(Gitter_ID_100m) %>%
  summarise(weight = sum(Anzahl))

q50p <- sum(dem50p$weight) / sum(pop$Anzahl) # overall share of people aged 50+

############################################################################## #
# Note: Artificial micro data for people aged 50+ is created by sampling from the 
# overall population; the spatial distribution of people aged 50+ is considered
# by means of different inclusion probabilities for grid cells.
############################################################################## #

# create inclusion probabilities
pop <- merge(pop, dem50p, by = "Gitter_ID_100m", all.x = TRUE)
pop$weight <- pop$weight / sum(pop$weight, na.rm = TRUE)
pop$weight[is.na(pop$weight)] <- 0

# sample
set.seed(20240422)
s <- sample(1:nrow(pop), 
            size = round(q50p * nrow(pop)), 
            replace = TRUE, 
            prob = pop$weight)

pop <- pop[s, ] %>% dplyr::select(-weight)
pop$person_id <- 1:nrow(pop)

# simulate exact locations
pop$x <- pop$x_mp_100m + runif(nrow(pop), -50, 50)
pop$y <- pop$y_mp_100m + runif(nrow(pop), -50, 50)


# ----- save data -----

save(pop, file = paste0(data_folder, "/pop_data.RData"), compress = TRUE)

rm(list = ls())


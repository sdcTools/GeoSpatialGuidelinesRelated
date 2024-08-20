##### R-Code for Case Study: Population grids with the Cell Key Method ####### #
#
# 02 - Apply CKM protection and analyse results
#
# contact: Martin Möhler, martin.moehler@destatis.de
# 2024-08-12
############################################################################## #

library(ptable)
library(raster)
library(data.table)
library(dplyr)
library(spdep)
library(sdcSpatial)
library(ggplot2)
library(SpatialKWD)

# customize to local setting
data_folder <- paste0(getwd(), "/geo_guidelines_case_study")

source(file.path(data_folder, "functions_ckm_raster.R"))
source(file.path(data_folder, "functions_infoloss_raster.R"))


# load data from step 01
load(file.path(data_folder, "pop_data.RData"))


# color palette
cbpal <- c("#0072B2", "#E69F00", "#56B4E9", "#D55E00")


# ----- set up data -----

# add record keys to micro data
set.seed(20240307)
pop$rk <- runif(nrow(pop))

# aggregate micro data to raster
r_dat <- sdc_raster_ckm(pop[, c("x", "y")], 
                        r = raster(resolution = 100, 
                                   xmn = min(pop$x_mp_100m) - 50,
                                   xmx = max(pop$x_mp_100m) + 50,
                                   ymn = min(pop$y_mp_100m) - 50,
                                   ymx = max(pop$y_mp_100m) + 50), 
                        variable = 1, min_count = 3,
                        rkey = pop$rk)

# sensitive cells according to sdcSpatial
plot(r_dat, "count", sensitive = TRUE)

# specify local focus regions
fa1 <- extent(r_dat$value, 425, 475, 25, 75)
fa2 <- extent(r_dat$value, 210, 260, 670, 720)
fa3 <- extent(r_dat$value, 925, 975, 475, 525)

# initial inspection

# share of cells at risk
sdcSpatial::sensitivity_score(r_dat)
risk_scores(crop(r_dat$value$count, fa1), type = "area", k = 3)
risk_scores(crop(r_dat$value$count, fa2), type = "area", k = 3)
risk_scores(crop(r_dat$value$count, fa3), type = "area", k = 3)

# share of pop. at risk
risk_scores(r_dat$value$count, type = "pop", k = 3) 
risk_scores(crop(r_dat$value$count, fa1), type = "pop", k = 3)
risk_scores(crop(r_dat$value$count, fa2), type = "pop", k = 3)
risk_scores(crop(r_dat$value$count, fa3), type = "pop", k = 3)

# map view
r_df <- as.data.frame(r_dat$value$count, xy = TRUE) %>% na.omit()
r_df$sens <- r_df$count > 0 & r_df$count < 3

ggplot(r_df, aes(x, y)) +
  geom_raster(aes(fill = count)) +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal() +
  ggtitle("LAEA 100kmN31E41") +
  xlab("E") +
  ylab("N") +
  geom_rect(aes(xmin = fa1@xmin, xmax = fa1@xmax, ymin = fa1@ymin, ymax = fa1@ymax), 
            fill = NA, color = "blue") +
  geom_rect(aes(xmin = fa2@xmin, xmax = fa2@xmax, ymin = fa2@ymin, ymax = fa2@ymax), 
            fill = NA, color = "blue") +
  geom_rect(aes(xmin = fa3@xmin, xmax = fa3@xmax, ymin = fa3@ymin, ymax = fa3@ymax), 
            fill = NA, color = "blue") +
  geom_raster(data = r_df[r_df$sens, ], fill = "red") +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_blank(), 
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))


# ----- apply CKM -----

# make p-tables (for 2 variants)
pt1 <- create_cnt_ptable(D = 4, V = 1.8, js = 1) # weaker protection
pt2 <- create_cnt_ptable(D = 6, V = 2.5, js = 2) # medium-level protection

pt1@pTable$variant <- "I"
pt2@pTable$variant <- "II"

# visualize ptables together
ptabs <- pt1@pTable[pt1@pTable$i > 0, ]
for(i in max(pt1@pTable$i):max(pt2@pTable$i)) {
  pt <- pt1@pTable[pt1@pTable$i == max(pt1@pTable$i), ]
  pt$i <- i
  pt$j <- pt$i + pt$v
  ptabs <- rbind(ptabs, pt)
}
ptabs <- rbind(ptabs, pt2@pTable[pt2@pTable$i > 0, ])
ptabs$i <- factor(ptabs$i, levels = 1:max(ptabs$i), 
                  labels = paste("i:", c(1:(max(ptabs$i) - 1), 
                                         paste0(max(ptabs$i), "+"))))

ptabs$v[ptabs$variant == "I"]  <- ptabs$v[ptabs$variant == "I"] - 0.12
ptabs$v[ptabs$variant == "II"] <- ptabs$v[ptabs$variant == "II"] + 0.12

ggplot(ptabs, 
       aes(x = v, y = p, group = variant)) +
  geom_segment(aes(xend = v, yend = 0, color = variant), alpha = 0.7, lwd = 1.0) +
  scale_fill_manual(values = c(cbpal[1], cbpal[2])) +
  scale_color_manual(values = c(cbpal[3], cbpal[4])) +
  facet_wrap(~i) +
  theme_bw() +
  xlab("noise value added") +
  ylab("probability") +
  theme(strip.background = element_rect(fill = cbpal[3]),
        legend.position = "bottom")


# perturbation step
r_ckm1 <- protect_ckm(r_dat, ptab = pt1)
r_ckm2 <- protect_ckm(r_dat, ptab = pt2)


# ----- inspect results -----

## share of cells < k after CKM

# variant I
risk_scores(r_ckm1$value$count, type = "area", k = 3)
risk_scores(crop(r_ckm1$value$count, fa1), type = "area", k = 3)
risk_scores(crop(r_ckm1$value$count, fa2), type = "area", k = 3)
risk_scores(crop(r_ckm1$value$count, fa3), type = "area", k = 3)
# variant II
risk_scores(r_ckm2$value$count, type = "area", k = 3)
risk_scores(crop(r_ckm2$value$count, fa1), type = "area", k = 3)
risk_scores(crop(r_ckm2$value$count, fa2), type = "area", k = 3)
risk_scores(crop(r_ckm2$value$count, fa3), type = "area", k = 3)

# map view
r_ckm_fa <- rbind(as.data.frame(crop(r_dat$value$count, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "original", fa = "focus area 1"),
                  as.data.frame(crop(r_dat$value$count, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "original", fa = "focus area 2"),
                  as.data.frame(crop(r_dat$value$count, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "original", fa = "focus area 3"),
                  as.data.frame(crop(r_ckm1$value$count, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 1"),
                  as.data.frame(crop(r_ckm1$value$count, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 2"),
                  as.data.frame(crop(r_ckm1$value$count, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 3"),
                  as.data.frame(crop(r_ckm2$value$count, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 1"),
                  as.data.frame(crop(r_ckm2$value$count, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 2"),
                  as.data.frame(crop(r_ckm2$value$count, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 3")) %>%
  na.omit()
r_ckm_fa$sens <- r_ckm_fa$count > 0 & r_ckm_fa$count < 3

ggplot(r_ckm_fa, aes(x, y)) +
  geom_tile(aes(fill = count)) +
  geom_tile(data = r_ckm_fa[r_ckm_fa$sens, ], fill = "red") +
  facet_grid(variant ~ fa) +
  scale_fill_viridis_c(direction = -1) +
  guides(fill = guide_colorbar(barheight = 15, barwidth = 0.5)) +
  scale_x_continuous(name = NULL, expand = c(0, 0)) +
  scale_y_continuous(name = NULL, expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))


# ----- info loss 1: measures of distributional distance -----

## MSE

# full map
distr_distance(r_ckm1$value$count, r_dat$value$count, "mse", rescale = FALSE, include_empty = FALSE)
distr_distance(r_ckm2$value$count, r_dat$value$count, "mse", rescale = FALSE, include_empty = FALSE)
# focus areas
distr_distance(crop(r_ckm1$value$count, fa1), crop(r_dat$value$count, fa1), 
               "mse", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa1), crop(r_dat$value$count, fa1), 
               "mse", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm1$value$count, fa2), crop(r_dat$value$count, fa2), 
               "mse", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa2), crop(r_dat$value$count, fa2), 
               "mse", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm1$value$count, fa3), crop(r_dat$value$count, fa3), 
               "mse", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa3), crop(r_dat$value$count, fa3), 
               "mse", rescale = FALSE, include_empty = FALSE)

## MAE

# full map
distr_distance(r_ckm1$value$count, r_dat$value$count, "mae", rescale = FALSE, include_empty = FALSE)
distr_distance(r_ckm2$value$count, r_dat$value$count, "mae", rescale = FALSE, include_empty = FALSE)
# focus areas
distr_distance(crop(r_ckm1$value$count, fa1), crop(r_dat$value$count, fa1), 
               "mae", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa1), crop(r_dat$value$count, fa1), 
               "mae", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm1$value$count, fa2), crop(r_dat$value$count, fa2), 
               "mae", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa2), crop(r_dat$value$count, fa2), 
               "mae", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm1$value$count, fa3), crop(r_dat$value$count, fa3), 
               "mae", rescale = FALSE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa3), crop(r_dat$value$count, fa3), 
               "mae", rescale = FALSE, include_empty = FALSE)

## HD

# full map
distr_distance(r_ckm1$value$count, r_dat$value$count, "hd", rescale = TRUE, include_empty = FALSE)
distr_distance(r_ckm2$value$count, r_dat$value$count, "hd", rescale = TRUE, include_empty = FALSE)
# focus areas
distr_distance(crop(r_ckm1$value$count, fa1), crop(r_dat$value$count, fa1), 
               "hd", rescale = TRUE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa1), crop(r_dat$value$count, fa1), 
               "hd", rescale = TRUE, include_empty = FALSE)
distr_distance(crop(r_ckm1$value$count, fa2), crop(r_dat$value$count, fa2), 
               "hd", rescale = TRUE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa2), crop(r_dat$value$count, fa2), 
               "hd", rescale = TRUE, include_empty = FALSE)
distr_distance(crop(r_ckm1$value$count, fa3), crop(r_dat$value$count, fa3), 
               "hd", rescale = TRUE, include_empty = FALSE)
distr_distance(crop(r_ckm2$value$count, fa3), crop(r_dat$value$count, fa3), 
               "hd", rescale = TRUE, include_empty = FALSE)

## KWD
# NOTE: This part can take some time to run!

# full map
get_KWD(r_ckm1$value$count, r_dat$value$count,
        type = "floor", method = "approx", L = 3)
get_KWD(r_ckm2$value$count, r_dat$value$count,
        type = "floor", method = "approx", L = 3)
# focus areas
get_KWD(crop(r_ckm1$value$count, fa1), crop(r_dat$value$count, fa1), 
        type = "floor", method = "exact")
get_KWD(crop(r_ckm2$value$count, fa1), crop(r_dat$value$count, fa1), 
        type = "floor", method = "exact")
get_KWD(crop(r_ckm1$value$count, fa2), crop(r_dat$value$count, fa2), 
        type = "floor", method = "exact")
get_KWD(crop(r_ckm2$value$count, fa2), crop(r_dat$value$count, fa2), 
        type = "floor", method = "exact")
get_KWD(crop(r_ckm1$value$count, fa3), crop(r_dat$value$count, fa3), 
        type = "floor", method = "exact")
get_KWD(crop(r_ckm2$value$count, fa3), crop(r_dat$value$count, fa3), 
        type = "floor", method = "exact")


# ----- info loss 2: measures of spatial association -----

## VMR

# full map
vmr(r_dat$value$count,  "vmr", include_empty = TRUE)
vmr(r_ckm1$value$count, "vmr", include_empty = TRUE)
vmr(r_ckm2$value$count, "vmr", include_empty = TRUE)
# focus areas
vmr(crop(r_dat$value$count,  fa1), "vmr", include_empty = TRUE)
vmr(crop(r_ckm1$value$count, fa1), "vmr", include_empty = TRUE)
vmr(crop(r_ckm2$value$count, fa1), "vmr", include_empty = TRUE)
vmr(crop(r_dat$value$count,  fa2), "vmr", include_empty = TRUE)
vmr(crop(r_ckm1$value$count, fa2), "vmr", include_empty = TRUE)
vmr(crop(r_ckm2$value$count, fa2), "vmr", include_empty = TRUE)
vmr(crop(r_dat$value$count,  fa3), "vmr", include_empty = TRUE)
vmr(crop(r_ckm1$value$count, fa3), "vmr", include_empty = TRUE)
vmr(crop(r_ckm2$value$count, fa3), "vmr", include_empty = TRUE)

## MMR

# full map
vmr(r_dat$value$count,  "mmr1", include_empty = FALSE)
vmr(r_ckm1$value$count, "mmr1", include_empty = FALSE)
vmr(r_ckm2$value$count, "mmr1", include_empty = FALSE)
# focus areas
vmr(crop(r_dat$value$count,  fa1), "mmr1", include_empty = FALSE)
vmr(crop(r_ckm1$value$count, fa1), "mmr1", include_empty = FALSE)
vmr(crop(r_ckm2$value$count, fa1), "mmr1", include_empty = FALSE)
vmr(crop(r_dat$value$count,  fa2), "mmr1", include_empty = FALSE)
vmr(crop(r_ckm1$value$count, fa2), "mmr1", include_empty = FALSE)
vmr(crop(r_ckm2$value$count, fa2), "mmr1", include_empty = FALSE)
vmr(crop(r_dat$value$count,  fa3), "mmr1", include_empty = FALSE)
vmr(crop(r_ckm1$value$count, fa3), "mmr1", include_empty = FALSE)
vmr(crop(r_ckm2$value$count, fa3), "mmr1", include_empty = FALSE)


## Moran's I

# full map
(mI_dat  <- get_moranI(r_dat$value$count,  type = "queen", zero.policy = TRUE))
(mI_ckm1 <- get_moranI(r_ckm1$value$count, type = "queen", zero.policy = TRUE))
(mI_ckm2 <- get_moranI(r_ckm2$value$count, type = "queen", zero.policy = TRUE))
# focus areas
(mI_dat_fa1  <- get_moranI(crop(r_dat$value$count,  fa1), type = "queen", zero.policy = TRUE))
(mI_ckm1_fa1 <- get_moranI(crop(r_ckm1$value$count, fa1), type = "queen", zero.policy = TRUE))
(mI_ckm2_fa1 <- get_moranI(crop(r_ckm2$value$count, fa1), type = "queen", zero.policy = TRUE))
(mI_dat_fa2  <- get_moranI(crop(r_dat$value$count,  fa2), type = "queen", zero.policy = TRUE))
(mI_ckm1_fa2 <- get_moranI(crop(r_ckm1$value$count, fa2), type = "queen", zero.policy = TRUE))
(mI_ckm2_fa2 <- get_moranI(crop(r_ckm2$value$count, fa2), type = "queen", zero.policy = TRUE))
(mI_dat_fa3  <- get_moranI(crop(r_dat$value$count,  fa3), type = "queen", zero.policy = TRUE))
(mI_ckm1_fa3 <- get_moranI(crop(r_ckm1$value$count, fa3), type = "queen", zero.policy = TRUE))
(mI_ckm2_fa3 <- get_moranI(crop(r_ckm2$value$count, fa3), type = "queen", zero.policy = TRUE))

## ILM

# full map
(mI_ckm1$estimate[1] - mI_dat$estimate[1]) / 2
(mI_ckm2$estimate[1] - mI_dat$estimate[1]) / 2
# focus areas
(mI_ckm1_fa1$estimate[1] - mI_dat_fa1$estimate[1]) / 2
(mI_ckm2_fa1$estimate[1] - mI_dat_fa1$estimate[1]) / 2
(mI_ckm1_fa2$estimate[1] - mI_dat_fa2$estimate[1]) / 2
(mI_ckm2_fa2$estimate[1] - mI_dat_fa2$estimate[1]) / 2
(mI_ckm1_fa3$estimate[1] - mI_dat_fa3$estimate[1]) / 2
(mI_ckm2_fa3$estimate[1] - mI_dat_fa3$estimate[1]) / 2


# ----- info loss 3: local information loss -----

## cell-level IL measures

# AD
r_ckm1$value$ad <- il_local(r_ckm1$value$count, r_dat$value$count, type = "ad")
r_ckm2$value$ad <- il_local(r_ckm2$value$count, r_dat$value$count, type = "ad")
# SD
r_ckm1$value$sd <- il_local(r_ckm1$value$count, r_dat$value$count, type = "sd")
r_ckm2$value$sd <- il_local(r_ckm2$value$count, r_dat$value$count, type = "sd")
# SDSR
r_ckm1$value$sdsr <- il_local(r_ckm1$value$count, r_dat$value$count, type = "sdsr")
r_ckm2$value$sdsr <- il_local(r_ckm2$value$count, r_dat$value$count, type = "sdsr")

# map view
r_ilcell <- rbind(as.data.frame(crop(r_ckm1$value$ad, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 1"),
                  as.data.frame(crop(r_ckm1$value$ad, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 2"),
                  as.data.frame(crop(r_ckm1$value$ad, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 3"),
                  as.data.frame(crop(r_ckm2$value$ad, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 1"),
                  as.data.frame(crop(r_ckm2$value$ad, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 2"),
                  as.data.frame(crop(r_ckm2$value$ad, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 3"))

ggplot(r_ilcell, aes(x, y, fill = ad)) +
  geom_tile() +
  facet_grid(variant ~ fa) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "AE   ") +
  guides(fill = guide_colorbar(barheight = 15, barwidth = 0.5)) +
  scale_x_continuous(name = NULL, expand = c(0, 0)) +
  scale_y_continuous(name = NULL, expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))


## Moving window approach

# MAE
r_ckm1$value$mw_mae <- mw_loss(r_ckm1$value$ad, wsize = 5, fun = mean)
r_ckm2$value$mw_mae <- mw_loss(r_ckm2$value$ad, wsize = 5, fun = mean)
# MSE
r_ckm1$value$mw_mse <- mw_loss(r_ckm1$value$sd, wsize = 5, fun = mean)
r_ckm2$value$mw_mse <- mw_loss(r_ckm2$value$sd, wsize = 5, fun = mean)
# HD
r_ckm1$value$mw_hd <- mw_loss(r_ckm1$value$ad, wsize = 5, 
                        fun = function(x){(1/sqrt(2))*sqrt(sum(x))})
r_ckm2$value$mw_hd <- mw_loss(r_ckm2$value$ad, wsize = 5, 
                        fun = function(x){(1/sqrt(2))*sqrt(sum(x))})

# map view
r_il_mw <- rbind(as.data.frame(crop(r_ckm1$value$mw_mae, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 1"),
                 as.data.frame(crop(r_ckm1$value$mw_mae, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 2"),
                 as.data.frame(crop(r_ckm1$value$mw_mae, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 3"),
                 as.data.frame(crop(r_ckm2$value$mw_mae, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 1"),
                 as.data.frame(crop(r_ckm2$value$mw_mae, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 2"),
                 as.data.frame(crop(r_ckm2$value$mw_mae, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 3"))

ggplot(r_il_mw, aes(x, y, fill = mw_mae)) +
  geom_tile() +
  facet_grid(variant ~ fa) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "MAE", limits = c(0, 6)) +
  guides(fill = guide_colorbar(barheight = 15, barwidth = 0.5)) +
  scale_x_continuous(name = NULL, expand = c(0, 0)) +
  scale_y_continuous(name = NULL, expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))


# focal MAE summary
fivenum(getValues(r_ckm1$value$mw_mae))
fivenum(getValues(r_ckm2$value$mw_mae))


## local Moran's I
# NOTE: This part can take some time to run!

r_dat$value$locI  <- r_dat$value$count
r_ckm1$value$locI <- r_ckm1$value$count
r_ckm2$value$locI <- r_ckm2$value$count
r_dat$value$locI  <- setValues(r_dat$value$locI,  
                               get_localmoranI(r_dat$value$count, type = "queen",
                                               zero.policy = TRUE)[, 1])
r_ckm1$value$locI <- setValues(r_ckm1$value$locI, 
                               get_localmoranI(r_ckm1$value$count, type = "queen",
                                               zero.policy = TRUE)[, 1])
r_ckm2$value$locI <- setValues(r_ckm2$value$locI, 
                               get_localmoranI(r_ckm2$value$count, type = "queen",
                                               zero.policy = TRUE)[, 1])

# map view
r_locIfa <- rbind(as.data.frame(crop(r_dat$value$locI, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "original", fa = "focus area 1"),
                  as.data.frame(crop(r_dat$value$locI, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "original", fa = "focus area 2"),
                  as.data.frame(crop(r_dat$value$locI, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "original", fa = "focus area 3"),
                  as.data.frame(crop(r_ckm1$value$locI, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 1"),
                  as.data.frame(crop(r_ckm1$value$locI, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 2"),
                  as.data.frame(crop(r_ckm1$value$locI, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant I", fa = "focus area 3"),
                  as.data.frame(crop(r_ckm2$value$locI, fa1), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 1"),
                  as.data.frame(crop(r_ckm2$value$locI, fa2), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 2"),
                  as.data.frame(crop(r_ckm2$value$locI, fa3), xy = TRUE) %>%
                    mutate(x = x - min(x), y = y - min(y), 
                           variant = "variant II", fa = "focus area 3")) %>%
  na.omit()

ggplot(r_locIfa, aes(x, y, fill = abs(locI))) +
  geom_tile() +
  facet_grid(variant ~ fa) +
  scale_fill_viridis_c(direction = 1, option = "mako", name = "loc. I (abs.)", trans = "log10") +
  guides(fill = guide_colorbar(barheight = 15, barwidth = 0.5)) +
  scale_x_continuous(name = NULL, expand = c(0, 0)) +
  scale_y_continuous(name = NULL, expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        strip.text = element_text(size = 13),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))


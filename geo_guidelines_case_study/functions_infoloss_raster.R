##### R-Code for Case Study: Population grids with the Cell Key Method ####### #
#
# Functions to measure information loss in protected grid data
#
# contact: Martin Möhler, martin.moehler@destatis.de
# 2024-04-29
############################################################################## #


##### Get measures of distributional distance for perturbed raster data ###### #
# x:             a raster object (protected data)
# y:             a raster object (unprotected data)
# type:          one of "mse", "mae", "hd"
# rescale:       whether to transform cell values to [0,1]
# include_empty: whether to include cell values of 0 in the measure
############################################################################## #

distr_distance <- function(x, y,
                           type = "mse", 
                           rescale = TRUE, 
                           include_empty = FALSE) {
  
  # extract cell values
  v_x <- raster::getValues(x)
  v_o <- raster::getValues(y)
  
  if(!include_empty) {
    
    # check which cells are non-empty in either raster
    nempt <- !(is.na(v_x) & is.na(v_o))
    
    v_x <- v_x[nempt]
    v_o <- v_o[nempt]
  }
  
  # set (remaining) NAs to 0
  v_x[is.na(v_x)] <- 0
  v_o[is.na(v_o)] <- 0
  
  # rescale values to [0,1] if required
  if(rescale) {
    v_x <- v_x / sum(v_x)
    v_o <- v_o / sum(v_o)
  }
  
  if(type == "mse") {
    
    dd <- mean((v_x - v_o)^2)
      
  } else if(type == "mae") {
    
    dd <- mean(abs(v_x - v_o))
    
  } else if(type == "hd") {
    
    dd <- sqrt(sum((sqrt(v_x) - sqrt(v_o))^2)) * (1/sqrt(2))
  }
  
  dd
}


##### Get Kantorovich-Wasserstein distance for perturbed raster data ######### #
# x:    a raster object (protected data)
# y:    a raster object (unprotected data)
# type: how to handle a mass gap - one of "rescale", "floor", or NULL
# ...   other parameters to be passed to SpatialKWD::CompareOneToOne()   
############################################################################## #

get_KWD <- function(x, y, type, ...) {
  
  # get grid & weights
  xy <- raster::xyFromCell(x, 1:ncell(x))
  v_x <- getValues(x)
  v_y <- getValues(y)
  
  # set NAs to 0
  v_x[is.na(v_x)] <- 0
  v_y[is.na(v_y)] <- 0
  
  # assess mass gap
  s_x <- sum(v_x, na.rm = TRUE)
  s_y <- sum(v_y, na.rm = TRUE)
  
  # Options of dealing with Type-I mass gap
  if(type == "rescale") {
    
    # proportional rescaling for multiplicative error
    v_x <- v_x * (s_y / s_x)
    
  } else if(type == "floor") {
    
    # disproportional rescaling for additive error
    v_x[v_x > 0] <- v_x[v_x > 0] + (s_y - s_x) / length(v_x[v_x > 0])
  }
  
  kwd <- SpatialKWD::compareOneToOne(Coordinates = xy, Weights = cbind(v_y, v_x), ...)
  kwd
}


##### Get local (grid cell-level) deviation measures ######################### #
# x:    a raster (protected data)
# y:    a raster (unprotected data)
# type: one of "ad", "sd", "sdsr"
############################################################################## #

il_local <- function(x, y, type = "ad") {
  
  x[is.na(x)] <- 0
  y[is.na(y)] <- 0
  
  if(type == "ad") {
    
    z <- abs(x - y)
    
  } else if(type == "sd") {
    
    z <- (x - y)^2
    
  } else if(type == "sdsr") {
    
    z <- (sqrt(x) - sqrt(y))^2
  }
  
  z
}


##### Determine Variance-Mean-Ratio and variants ############################# #
# x:             a raster object 
# type:          one of "vmr", "mmr", "mmr1", "mmmr2
# include_empty: whether to include cell values of 0 in the measure
############################################################################## #

vmr <- function(x, type = "vmr", include_empty = TRUE) {
  
  v_x <- getValues(x)
  
  v_x[is.na(v_x)] <- 0 # NAs to 0
  if(!include_empty) {
    v_x <- v_x[v_x != 0]
  }
  
  if(type == "vmr") {
    
    vmr <- var(v_x) / mean(v_x)
    
  } else if(type == "mmr") {
    
    vmr <- median(abs(v_x - median(v_x))) / median(v_x)
  
  } else if(type == "mmr1") {
    
    vmr <- median(abs(v_x - mean(v_x))) / median(v_x)
    
  } else if(type == "mmr2") {
    
    vmr <- mean(abs(v_x - median(v_x))) / median(v_x)
    
  }
  
  vmr
}


##### Determine Moran's I #################################################### #
# x:    a raster object 
# type: spatial weights definition - one of "queen" (default) or "rook"
# ...   other parameters to be passed to spdep::moran.test()
############################################################################## #

get_moranI <- function(x, type = "queen", ...) {
  
  nb <- spdep::cell2nb(nrow(x), ncol(x), type = type)
  y <- getValues(x)
  y[is.na(y)] <- 0
  spdep::moran.test(y, nb2listw(nb, style = "W"), ...)
}

get_localmoranI <- function(x, type = "queen", ...) {
  
  nb <- spdep::cell2nb(nrow(x), ncol(x), type = type)
  y <- getValues(x)
  y[is.na(y)] <- 0
  spdep::localmoran(y, nb2listw(nb, style = "W"), ...)
}


##### Moving window / focal error ############################################ #
# x:     a raster object (of cell-level distance measures)
# wsize: window size (odd number, default: 5)
# fun:   focal function (typically mean)
# ...    other parameters to be passed to raster::focal()
############################################################################## #

mw_loss <- function(x, wsize = 5, fun = mean, ...) {
  
  W <- matrix(1, wsize, wsize)
  x[is.na(x)] <- 0
  raster::focal(x, w = W, fun = fun, pad = TRUE, padValue = 0, ...)
}


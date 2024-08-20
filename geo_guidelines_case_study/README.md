Code for Caste Study: Population grids with the Cell Key Method

These R scripts are related material for the 'Guidelines for Statistical Disclosure Control Methods Applied on Geo-Referenced Data' (2024), specifically Chapter 6.

Run in order:

01_prepare_data.R
  - download source data files
  - crop source data 
  - create synthetic data for analysis

02_protect_and_analyse.R
  - apply Cell Key Method to data
  - assess protection efficacy
  - measure information loss w.r.t. distributional distance
  - measure information loss w.r.t. spatial association
  - measure information loss w.r.t. local indicators

functions_ckm_raster.R is not run directly, but sourced in 02_protect_and_analyse.R; it contains helper functions for using the Cell Key Method on population grids.

functions_infoloss_raster.R is also not run directly, but sourced in 02_protect_and_analyse.R; it contains helper functions for measuring information loss on grid aggregates.

Contact: Martin Möhler (martin.moehler@destatis.de)

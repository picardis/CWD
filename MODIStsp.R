# Download NDVI using MODIStsp

# Load packages ####
library(MODIStsp)
library(terra)
library(sf)
library(tidyverse)

# Find product ####

MODIStsp_get_prodnames()
# Jerod's code uses "MOD09Q1"
# That is not listed here. Could try that first, otherwise
# "Vegetation Indexes_16Days_250m (M*D13Q1)"

MODIStsp_get_prodlayers("M*D09Q1")$bandnames
MODIStsp_get_prodlayers("M*D13Q1")$bandnames

MODIStsp_get_prodlayers("M*D13Q1")$quality_bandnames

MODIStsp_get_prodlayers("M*D13Q1")$indexes_bandnames


# # CRS conversion ####
#
# modis_sinusoidal <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m"
#
# sa <- ext(c(343428.7, 831738.7, 4029947, 4704497))
# sa_rast <- rast(extent = sa, crs = "epsg:32612")
# sa_sinus <- project(sa_rast, modis_sinusoidal)
# ext(sa_sinus)

# Download data ####

## NDVI ####

# This did not work; extent does not come out right

# MODIStsp(gui             = FALSE,
#          out_folder      = 'output/processed_layers/MODIS',
#          out_folder_mod  = 'output/processed_layers/MODIS',
#          selprod         = 'Vegetation Indexes_16Days_250m (M*D13Q1)',
#          bandsel         = 'NDVI',
#          sensor          = 'Terra',
#          user            = 'simona.picardi' , # your username for NASA http server
#          password        = '1990AserejeNASA!',  # your password for NASA http server
#          start_date      = '2016.01.01',
#          end_date        = '2021.01.12',
#          spatmeth        = 'bbox',
#          bbox            =  c(343428.7, 831738.7, 4029947, 4704497),
#          out_format      = 'GTiff',
#          compress        = 'LZW',
#          out_projsel     = 'User Defined',
#          output_proj     = 'epsg:32612',
#          verbose         = TRUE,
#          delete_hdf      = TRUE,
#          parallel        = TRUE
# )

MODIStsp(gui             = FALSE,
         out_folder      = 'output/processed_layers/MODIS',
         out_folder_mod  = 'output/processed_layers/MODIS',
         selprod         = 'Vegetation Indexes_16Days_250m (M*D13Q1)',
         bandsel         = 'NDVI',
         sensor          = 'Terra',
         user            = 'simona.picardi' , # your username for NASA http server
         password        = '1990AserejeNASA!',  # your password for NASA http server
         start_date      = '2016.01.01',
         end_date        = '2021.01.12',
         spatmeth        = 'file',
         spafile         = 'output/processed_layers/elevation.tiff',
         out_format      = 'GTiff',
         compress        = 'LZW',
         out_projsel     = 'User Defined',
         output_proj     = 'epsg:32612',
         verbose         = TRUE,
         delete_hdf      = TRUE,
         parallel        = TRUE
)

## Snow ####

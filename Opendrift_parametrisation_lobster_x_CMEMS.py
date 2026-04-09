#!/usr/bin/env python
"""
CMEMS current components
========================

The global CMEMS/copernicus ocean model contains several different surface current fields.
OpenDrift will by default use the variables with "standard_name" attribute equal to what is requested from a module, typically "x_sea_water_velocity"/"y_sea_water_velocity" for currents.
Note that "east/north" counterparts will also be detected, and eventual rotation will be performed automatically.

This example illustrates how a "standard_name_mapping" can be added to the generic netCDF reader to chose alternative variables.
The example also illustrates the alternative (experimental) mechanism of summing two readers.
"""

import numpy as np
import sys,  os
from datetime import datetime, timedelta
import cf_xarray
import logging
import copernicusmarine
from opendrift.models.oceandrift import OceanDrift, Lagrangian3DArray
from numpy.random import RandomState
from opendrift.readers.reader_netCDF_CF_generic import Reader
import matplotlib.pyplot as plt
import csv
from opendrift.models.IBM_BITER_NEPHROPS_CMEMS import PelagicShrimpDrift



#%%
# Set Parameters to modify
startTime=datetime(2011,1,1,0,0)
outputFilename='IBM_Nep_cmems_01_01_2011_t.nc' # Modify to fit with particle release date.
time=["2011-01-01","2011-02-15"] # Modify to make December, January, February releases and get 45 days of hydrodynamics


# Get an Xarray dataset from copernicusmarine client
ds = copernicusmarine.open_dataset(dataset_id='med-cmcc-cur-rean-d',username = "IDENTIFIANT",password="PWD", minimum_longitude = -1, maximum_longitude = 7, minimum_latitude = 38, maximum_latitude = 44, start_datetime = time[0], end_datetime = time[1], minimum_depth = 1.0182366371154785, maximum_depth = 3000)
ds_temp = copernicusmarine.open_dataset(dataset_id='med-cmcc-tem-rean-d',username = "IDENTIFIANT",password="PWD", minimum_longitude = -1, maximum_longitude = 7, minimum_latitude = 38, maximum_latitude = 44, start_datetime = time[0], end_datetime = time[1], minimum_depth = 1.0182366371154785, maximum_depth = 3000)
ds_sal = copernicusmarine.open_dataset(dataset_id='med-cmcc-sal-rean-d',username = "IDENTIFIANT",password="PWD",  minimum_longitude = -1, maximum_longitude = 7, minimum_latitude = 38, maximum_latitude = 44, start_datetime = time[0], end_datetime = time[1], minimum_depth = 1.0182366371154785, maximum_depth = 3000)
ds_bathy = copernicusmarine.open_dataset(dataset_id='cmems_mod_med_phy_my_4.2km_static',username = "IDENTIFIANT",password="PWD",  minimum_longitude = -1, maximum_longitude = 7, minimum_latitude = 38, maximum_latitude = 44, minimum_depth = 1.0182366371154785, maximum_depth = 3000)
ds_mld = copernicusmarine.open_dataset(dataset_id='med-cmcc-mld-rean-d',username = "IDENTIFIANT",password="PWD", minimum_longitude = -1, maximum_longitude = 7, minimum_latitude = 38, maximum_latitude = 44, start_datetime = time[0], end_datetime = time[1], minimum_depth = 1.0182366371154785, maximum_depth = 3000)


print(ds)     # Default Xarray output
print(ds.cf)  # Output from cf-xarray

reader_curr = Reader(ds, name='CMEMS current')
reader_temp = Reader(ds_temp, name='CMEMS temp')
reader_sal = Reader(ds_sal, name='CMEMS sal')
reader_bathy = Reader(ds_bathy, name='CMEMS bathy')
#reader_bathy.activate_environment_mapping('land_binary_mask_from_ocean_depth')
reader_mld = Reader(ds_mld, name='CMEMS mld')

o = PelagicShrimpDrift(loglevel=0)

# Configuration
o.set_config('general:seafloor_action', 'lift_to_seafloor')
o.set_config('general:coastline_action','stranding')
o.set_config('drift:deactivate_south_of',38.0)
o.set_config('drift:advection_scheme','runge-kutta4')
o.set_config('drift:vertical_advection',True)
o.set_config('drift:vertical_mixing', True)
o.set_config('vertical_mixing:TSprofiles', False)
o.set_config('vertical_mixing:timestep', 60) # seconds
o.set_config('vertical_mixing:diffusivitymodel','constant')
#o.set_config('drift:horizontal_diffusivity', 0.0001)
o.set_config('drift:vertical_swimming',True)

# Readers
o.add_reader(reader_curr  ,variables=['x_sea_water_velocity', 'y_sea_water_velocity'])
o.add_reader(reader_temp ,variables=['sea_water_temperature'])
o.add_reader(reader_sal ,variables=['sea_water_salinity'])
o.add_reader(reader_bathy)
o.add_reader(reader_mld ,variables=['ocean_mixed_layer_thickness'])


# Seeding
N=566115

# Below is for seeding particles in the MPAs from a text file
lon=[]
lat=[]

with open('~/path_to_folder/IBM_NEP_Med.txt') as f: ### Modify Path to Particle initialisation
      reader=csv.reader(f, delimiter='\t')
      for row in reader:
          lon=np.append(lon,float(row[0]))
          lat=np.append(lat,float(row[1]))
          if len(lon)== N:
              break
o.seed_elements(lon, lat, z='seafloor + 2.0',time=startTime)


# Run the model
#########################
o.run(time_step=3*3600,time_step_output=3600*24,duration=timedelta(days=45), outfile=outputFilename)


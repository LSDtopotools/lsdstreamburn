"""Main module."""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#import lsdtopytools as lsd
import lsdviztools.lsdbasemaptools as bmt
from lsdviztools.lsdplottingtools import lsdmap_gdalio as gio
import lsdviztools.lsdmapwrappers as lsdmw

def get_dem(OT_api_key_fname="my_OT_api_key.txt", source = "COP30",
            lower_left = [-16.88,-65.76],
            upper_right = [-16.50,-65.37],
            prefix = "Bolivia"):
    """Downloads the DEM

    Args:
        OT_api_key_fname (str): the filename where your opentopography api key is stored
        source (str): the type of DEM data (see opentopography.org for options)
        lower_left (list): a list of two floats that have the latitude and longitude of the lower left corner (you can copy this from google maps)
        upper_right (list): a list of two floats that have the latitude and longitude of the upper right corner (you can copy this from google maps)
        prefix (str): the prefix if the DEM name
        
    Return:
        the full filename of the DEM

    Author: SMM
    """  
    
    test_DEM = bmt.ot_scraper(source = source,
                        lower_left_coordinates = lower_left,
                        upper_right_coordinates = upper_right,
                        prefix = prefix,
                        api_key_file = OT_api_key_fname)

    fullfname,DataDirectory,BurnedDEM = test_DEM.download_pythonic() 
    return fullfname

def stream_burn(DataDirectory="./", DEM_prefix ="Bolivia", DEM_source = "COP30", location_year = "Bolivia2020",  
                burn_sediment=True, depth1=40, depth2=5):
    """Burns streams
 
    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.
        
    Return: 
        dem_utm (xarray raster): the conditioned DEM in UTM coordinates

    Author: QC and SMM
    """        


    # First, load the class map using xarray 
    print("Loading the class map")
    classmap_path = DataDirectory+location_year+"_example_ClassMap.tif"
    channel_mask = xr.open_dataarray(classmap_path)

       
    # Now load the DEM using xarray
    print("Loading the DEM")
    dem_path = DEM_prefix+"_"+DEM_source+'.tif'
    dem = xr.open_dataarray(dem_path)    
   
    # Reproject channel mask to match the CRS system of channel mask
    channel_mask_reprojected = dem.rio.reproject_match(dem)     
    
    # Set the output directory
    out_path = './burned_dem'
    isExist = os.path.exists(out_path)
    if not isExist:
        os.makedirs(out_path)
        print(f"I made a new directory called {out_path}")
        
    out_full_fname = out_path
 
    # Create a new DEM using the water and sediment mask
    dem1 = xr.where(channel_mask_reprojected == 1, dem - depth1, dem)
    if burn_sediment:
        dem1 = xr.where(channel_mask_reprojected == 2, dem1 - depth2, dem1)
        dem2 = dem1.rio.set_crs(dem.rio.crs)
        dem_utm = dem2.rio.reproject(dem2.rio.estimate_utm_crs()) # nodata type is default: positive inf

        dem_utm.rio.to_raster(os.path.join(out_path, f'{location_year}_{DEM_source}_Burned{depth1}m{depth2}m_UTM.tif'))
        out_full_fname = os.path.join(out_path, f'{location_year}_{DEM_source}_Burned{depth1}m{depth2}m_UTM.tif')
    else:
        dem2 = dem1.rio.set_crs(dem.rio.crs)
        dem_utm = dem2.rio.reproject(dem2.rio.estimate_utm_crs()) # nodata type is default: positive inf
        
        dem_utm.rio.to_raster(os.path.join(out_path, f'{location_year}_{DEM_source}_Burned{depth1}m_UTM.tif'))
        out_full_fname = os.path.join(out_path, f'{location_year}_{DEM_source}_Burned{depth1}m_UTM.tif')

    return [out_full_fname,dem_utm]


def extract_network_lsdtopytools(DataDirectory, BurnedDEM, area_thresh=15000):
    """
    mydem = lsd.LSDDEM(path = DataDirectory, 
                       file_name = BurnedDEM, 
                       already_preprocessed = False) # remove_seas=True, sea_level = 0.01, can also set it to the minimum elevation of DEM
    mydem.PreProcessing()

    # If getting HDF5EXTError, add a line in ~/.bashrc file: export HDF5_USE_FILE_LOCKING='FALSE'

    mydem.CommonFlowRoutines()
    mydem.ExtractRiverNetwork(method = "area_threshold", area_threshold_min = area_thresh)


    # Export river network
    network = mydem.df_base_river.set_index(['y','x']).to_xarray().set_coords(['x','y'])  #set_index() here will assign the dimensions for xarray dataset
    
    network_outpath = './network'
    isExist = os.path.exists(network_outpath)
    if not isExist:
        os.makedirs(network_outpath)
        print(f"A new directory {network_outpath} is created.")

    # Assign the crs of burned DEM to the river network
    dem_burned = xr.open_dataarray(os.path.join(DataDirectory, BurnedDEM))
    network.rio.write_crs(dem_burned.rio.crs, inplace=True)

    name, ext = os.path.splitext(BurnedDEM)
    new_name = "{name}_{area_thresh}_network{ext}".format(name=name, area_thresh=area_thresh, ext=ext)

    network['elevation'].rio.to_raster(os.path.join(network_outpath, new_name))
    print(f"The extracted river network Tiff file is exported to {network_outpath} folder.")

    drainage_area = network['drainage_area']
    print("Max_drainage_area:", drainage_area.max())

    return network    
    """    
    
    

def extract_network(DataDirectory, BurnedDEM, area_thresh=15000):
    """Extracts network

    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.

    Author: SMM
    """     
 

    # first figure out if lsdtopotools is available
    from shutil import which

    if (which("lsdtt-basic-metrics") is None):
        print("I did not find lsdtt-basic-metrcs. You need to install that.")
        exit(0)
    else: 
        print("Good news, I found lsdtt-basic-metrics. Lets go!")
        

    gio.convert4lsdtt(DataDirectory,BurnedDEM)

    area_thresh_string = str(area_thresh)   


    ## Get the basins
    lsdtt_parameters = {"remove_seas" : "true",
                    "threshold_contributing_pixels" : area_thresh_string,
                    "print_channels_to_csv" : "true"}
    r_prefix = DataDirectory+"_"+BurnedDEM +"_UTM"
    w_prefix = DataDirectory+"_"+BurnedDEM +"_UTM"
    
    print("The read prefix is: "+r_prefix)
    lsdtt_drive = lsdmw.lsdtt_driver(command_line_tool = "lsdtt-basic-metrics",
                                 read_prefix = r_prefix,
                                 write_prefix= w_prefix,
                                 read_path = "./",
                                 write_path = "./",
                                 parameter_dictionary=lsdtt_parameters)
    lsdtt_drive.print_parameters()
    lsdtt_drive.run_lsdtt_command_line_tool()    
    

def plot_network():
    print("Not finished yet")


def burning_driver(DataDirectory = "./", 
                   location_year = 'Bolivia2020', 
                   area_thresh=15000):
    """Extracts network

    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.

    Author: QC
    """  
    

    
    # Reproject DEM to match the CRS system of channel mask using rasterio
    dem_reprojected = dem.rio.reproject_match(channel_mask)

    # Burn water 20 meters and fluvial sediment 15 meters on the reprojected DEM
    dem_burned = stream_burn(location_year, channel_mask, 'COP30', dem_reprojected, True, 20, 15)
    print("The minimum value of the burned DEM is: ", np.nanmin(dem_burned.data))

    # Define a drainage area threshold for the river network
    threshold_area = 15000

    # Located burned DEM
    burned_dem_path = './burned_dem/'

    dem_list = [f for f in os.listdir(burned_dem_path) if f.endswith('.tif')]
    print("Available burned DEM: ", dem_list)

    # Extract the network
    Extract_Network(burned_dem_path, dem_list[0], threshold_area)

    
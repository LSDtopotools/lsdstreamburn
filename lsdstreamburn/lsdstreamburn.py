"""Main module."""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import lsdtopytools as lsd
import lsdviztools.lsdbasemaptools as bmt

def get_dem(OT_api_key_fname="my_OT_api_key.txt", source = "COP30",
            lower_left = [-16.88,-65.76],
            upper_right = [-16.50,-65.37],
            prefix = "Bolivia"):
    """Downloads the DEM

    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.

    Author: SMM
    """  
    
    test_DEM = bmt.ot_scraper(source = source,
                        lower_left_coordinates = lower_left,
                        upper_right_coordinates = upper_right,
                        prefix = prefix,
                        api_key_file = OT_api_key_fname)

    fullfname,DataDirectory,BurnedDEM = test_DEM.download_pythonic() 
    return fullfname

def stream_burn(location_year, water, 
                dem_type, dem, 
                burn_sediment=True, depth1=40, depth2=5):
    """Burns streams
api_key 
    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.
        
    Return: 
        dem_utm (xarray raster): the conditioned DEM in UTM coordinates

    Author: QC and SMM
    """        
    
    # Set the output directory
    out_path = './burned_dem'
    isExist = os.path.exists(out_path)
    if not isExist:
        os.makedirs(out_path)
        print(f"A new directory {out_path} is created.")
 
    # Create a new DEM using the water and sediment mask
    dem1 = xr.where(water == 1, dem - depth1, dem)
    if burn_sediment:
        dem1 = xr.where(water == 2, dem1 - depth2, dem1)
        dem2 = dem1.rio.set_crs(dem.rio.crs)
        dem_utm = dem2.rio.reproject(dem2.rio.estimate_utm_crs()) # nodata type is default: positive inf

        dem_utm.rio.to_raster(os.path.join(out_path, f'{location_year}_{dem_type}_Burned{depth1}m{depth2}m_UTM.tif'))
    else:
        dem2 = dem1.rio.set_crs(dem.rio.crs)
        dem_utm = dem2.rio.reproject(dem2.rio.estimate_utm_crs()) # nodata type is default: positive inf
        
        dem_utm.rio.to_raster(os.path.join(out_path, f'{location_year}_{dem_type}_Burned{depth1}m_UTM.tif'))

    return dem_utm


def extract_network(DataDirectory, BurnedDEM, area_thresh=15000):
    """Extracts network

    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.

    Author: QC
    """     
    
    mydem = lsd.LSDDEM(path = DataDirectory, 
                       file_name = BurnedDEM, 
                       already_preprocessed = False) # remove_seas=True, sea_level = 0.01, can also set it to the minimum elevation of DEM
    mydem.PreProcessing()

    # If getting HDF5EXTError, add a line in ~/.bashrc file: export HDF5_USE_FILE_LOCKING='FALSE'

    mydem.CommonFlowRoutines()
    mydem.ExtractRiverNetwork(method = "area_threshold", area_threshold_min = area_thresh)

    catch_coord = mydem.DefineCatchment(method = "main_basin")
    mydem.GenerateChi(theta = 0.35, A_0 = 1)

    # Plot rivernetwork on DEM
    fig, ax = lsd.quickplot.get_basemap(mydem , figsize = (15,15), cmap = "gist_earth", hillshade = True, 
	alpha_hillshade = 1, cmin = None, cmax = None,
	hillshade_cmin = 0, hillshade_cmax = 1, colorbar = False, 
	fig = None, ax = None, colorbar_label = None, colorbar_ax = None, fontsize_ticks = 16, normalise_HS = True)
    size_array = lsd.size_my_points(np.log10(mydem.df_base_river.drainage_area), 1, 5)

    ax.scatter(mydem.df_base_river.x, mydem.df_base_river.y, lw=0, c= "b",  zorder = 2, s=size_array, alpha=0.8)
    lsd.quickplot_utilities.add_basin_outlines(mydem, fig, ax, size_outline = 20, zorder = 5, color = "k")

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





def burning_driver(DataDirectory = "./", 
                   location_year = 'Bolivia2020', 
                   area_thresh=15000):
    """Extracts network

    Args:
        ax_list (axes objects): the list of axis objects
        axis_style (string): The syle of the axis. See options below.

    Author: QC
    """  
    
    # First, load the class map using xarray 
    classmap_path = DataDirectory+location_year+"_example_ClassMap.tif"
    channel_mask = xr.open_dataarray(classmap_path)
    
    # Now load the DEM using xarray
    dem_path = 'Bolivia_DEM30_COP30.tif'
    dem = xr.open_dataarray(dem_path)
    
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

    
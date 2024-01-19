============
lsdstreamburn
============

Welcome to the "lsdstreamburn" repository, focused on efficient river network extraction using advanced stream burning techniques. Ready for a quick tour?

Purpose
-------

This space is dedicated to sharing insights into river network extraction. Consult our stream burning paper and refer to an example shown in "stream_burn_example.ipynb" for detailed instructions and code. Explore, experiment, and adapt the code to suit your needs.

Dive In!
--------

The "lsdstreamburn.py" file under "lsdstreamburn" directory contains code for seamless river network extraction. The "stream_burn_example.ipynb" demonstrates the usage of this package. The "gee_codes" directory contains Google Earth Engine scripts to obtain a supervised land cover classification map.

Prerequisites
-------------

1. **Register on Earth Engine Code Editor:** [Create an account](https://code.earthengine.google.com/) on Earth Engine Code Editor to get started.

2. **Obtain OpenTopography API Key:** Follow [this guide](https://opentopography.org/blog/introducing-api-keys-access-opentopography-global-datasets) to register an OpenTopography account and obtain an API key. This key is essential for customizing your DEM.

3. **Python Environment Setup:** In a Python3 environment, make sure to install the following packages:
    - `numpy`
    - `xarray`
    - `rioxarray`
    - `matplotlib`

   You can install them using the following command:

   pip install numpy xarray rioxarray matplotlib

4. **LSDTopoTools & LSDVizTools Installation:**  In the same Python environment, install 'LSDTopoTools2' and 'lsdviztools' from their [GitHub repository](https://github.com/LSDtopotools) using the following commands:

    conda install -c conda-forge lsdtopotools
    conda install -y lsdviztools

    (This step can be a bit complicated. Please follow their installing instruction very closely!)

Let's Get Rolling!
------------------

1. **Clone & Download:** Grab a copy of this repository.

2. **Map Land Cover:** Before you start running stream burning, make sure you have a land cover classification map ready. It's a map with class values: 1 for water pixels and 2 for fluvial sediment. Need help creating one? Try out the classification script "GEE_supervised_classification_Sentinel_1_2.js". Open and Copy the codes, and paste them to a new and empty file in the Earth Engine Code Editor: https://code.earthengine.google.com/. The script cannot be run successfully without your input of training data and validating data. At line 93-100, you are asked to create classes of vegetation, water surface, fluvial sediment and a joint class of urban and agricultural lands (named as "vege", "water", "sedi", "urban_agri" in the script). Then you need to draw polygons for each class on the map following section 1.2 from this tutorial: https://kraaijenbrink.github.io/earthengine-workshop/image-classification.html. At line 155-163, you need to repeat the data collection process for validating data, with each class named "Vvege", "Vwater", "Vsedi", "Vurban_agri". After you finish collecting data for both training and validating, this script should be able to run successfully.

3. **Denoise Land Cover Map:** This step will help you to denoise the land cover classification map from last step. Open the script "GEE_post_classification.js" and copy everything. Open another new file in the Earth Engine code editor, paste all the codes in. Run it and download the exported clustered map.
Make sure your land cover classification map of a river basin is ready. If not, we prepare you an example here called "Bolivia2020_example_ClassMap.tif" that can be downloaded from [a Dropbox storage](https://www.dropbox.com/scl/fi/cqshbenb7xjp62c1h4u6o/Bolivia2020_example_ClassMap.tif?rlkey=m9jejwt4qaqdqpjgv3wtd906d&dl=0).

4. **Stream burn in a Jupyter Notebook:** Open "stream_burn_example.ipynb". Follow the notebook's guidance to run cells.

5.. **Check the Result** Download a Sentinel-2 Red-Green-Blue image from the same area of interest, and overlay the extracted river network on the Sentinel-2 image to check the alignment.

Downloading DEM
---------------

Want to tweak the DEM (Digital Elevation Model)? The notebook has a cool cell to download a DEM from OpenTopography API. Customize it by playing around with longitude, latitude, and "source" parameters! It can download all the DEM data that is ported by OpenTopography. Be aware that the FABDEM is not included, which you have to download manually or through Google Earth Engine.

License
-------

This place follows the [MIT License](LICENSE). Grab the code, play with it, and make it yours.

Happy stream burning! 🌊✨
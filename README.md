# Relationship between Arctic Sea-ice Morphology and Greenhouse Gases with SENTINEL-2 Images
This program plots flux footprints (using the simple parameterization flux footprint program of Kljun et al, 2015).

This package consist of 
1. Footprint_All_Image_overlay.m 
  This is the main MATLAB program to compute footprints over Satellite Image (SENTINEL-2). This program exports the flux footprint model of Kljun et al (2015) to shapefiles. It also plots SENTINEL-2 images under the footprints along with a rose diagram for wind direction, and plots for wind speed (u), Friction velocity (ustar), Latent Heat (H), sensible heat (E), and CO2 flux (fco2).  
2. donut_makge_new_single.py
  This is a PYTHON program that uses the ESRI proprietary Arcpy module to convert shapefiles to calculate the contribution of land, sea-ice, and open water to the footprint models. 

##################################################
##     ##     ##    ## ### ##     ##    ###     ##
#### #### ### ## ## ##  ## ## ### ## ### ## ### ##
#### #### ### ##   ### # # ## ### ## ### ## ### ##
#### #### ### ## ## ## ##  ##     ## ### ## ### ##
#### ####     ## ## ## ### ## ### ##    ###     ##
##################################################

## Acclimation to astropy and fits (flexible image transport system)
## Start Date: 2023, January, 30

## Links to the Astropy official website tutorials
# https://learn.astropy.org/tutorials/FITS-tables.html # FITS table handling
# https://learn.astropy.org/tutorials/celestial_coords1.html # Celestial coordinates

## TABLE OF CONTENT

## 1-) Opening FITS file
## 2-) Retrieving headers and managing it
## 3-) Plotting image data from FITS file
## 4-) Working with WCS, creating it, retrieving it, and plotting according to it

##### Part 0: Library Loading

import numpy as np # Numerical python

from astropy.io import fits # fits management
from astropy.table import Table # header/Table management from fits
from astropy.wcs import WCS # Working with celestial coordinates
from astropy.utils.data import download_file # ?? Astropy utility to download fits

from matplotlib.colors import LogNorm # ???
import matplotlib.pyplot as plt # Importing pyplot

##### Part 0.5: Tutorial specific file

# This file is from Chandra Observation, X ray with coordinate information

# cache parameter in the following line allows one to use the file in cache when present
# preventing unnecessary downloads
event_filename = download_file('http://data.astropy.org/tutorials/FITS-tables/chandra_events.fits', cache=True)

##### Part 1: Opening the FITS file

hdu_list = fits.open(event_filename, memmap=True) # ?memmap is memory related tweak, hdu: header unit
hdu_list.info() # This command summarized the information in fits header in a table

print(hdu_list[1].columns) # This prints the columns of the 1st header of the hdu_list

##### Part 2: Conversion to AstropyTable, convenience

evt_data = Table(hdu_list[1].data) # function of astropy.table pack
print(evt_data)

##### Part 3: Histogram creation, 1D and 2D

#1D
plt.figure() # Prepares an instance to show figure after plt operation
energy_hist = plt.hist(evt_data['energy'], bins='auto') # hist function of PyPlot, energy column of event data
plt.show() # actually plots the figure to the screen

#2D

# np.in1d function generates a boolean array with the dimensions of the inpt array
# which return True when the second array also contain the element, and 0 otherwise.
ii = np.in1d(evt_data['ccd_id'], [0,1,2,3]) # 0,1,2,3 ids correspond to the main chips, ACIS-I

NBINS = (100,100) # number of bins for two dimensions

# np.histogram2d is 2D histogram, will separately histogram the 'x' and 'y' parts of the ii, with NBINS (nx,ny) consideration
img_zero, yedges, xedges = np.histogram2d(evt_data['x'][ii], evt_data['y'][ii], NBINS)

extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]] # The image extent
plt.figure()
plt.imshow(img_zero, extent=extent, 
                     interpolation='nearest', 
                     cmap='gist_yarg', 
                     origin='lower') # Way to plot a 2D "raster"
                     
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.figure()
NBINS = (100,100)
img_zero_mpl = plt.hist2d(evt_data['x'][ii], evt_data['y'][ii], NBINS, 
                          cmap='viridis', norm=LogNorm())

cbar = plt.colorbar(ticks=[1.0,3.0,6.0])
cbar.ax.set_yticklabels(['1','3','6'])

plt.xlabel('x')
plt.ylabel('y')
plt.show()

##### Part 4: Closing the fits file

hdu_list.close()

##### Part 5: Working with Celestial Coordinates - Intro

### Important general keywords: 

# pixel: part of a detector corresponding to a real 2D area
# voxel: A 3D counterpart of pixel, indicating "volume"
# spaxel: spatial or spectral pixel
# gnomonic: tangent-plane

### WCS keywords:

# CRVALn: Actual coordinate of the information present in the pixel/voxel
# CRPIXn: Location in terms of pixel order, CRPIX1=1 and CRPIX2=1 represents a corner in a 2D image
# CDELTn: Delta (coordinate difference) of distance from one pixel to its neighbor
# CTYPEn: Eight-char string for axis type, e.g., 'RA---TAN' 'DEC--TAN' 
# CUNITn: String for unit, default is degree
# NAXISn: Number of pixels in each axises

##### Part 6: Working with Celestial Coordinates - Creating a WCS object with a dict

wcs_input_dict = {
    'CTYPE1': 'RA---TAN', # The first axis is RA, right ascension
    'CUNIT1': 'deg', 
    'CDELT1': -0.0002777777778, # This corresponds to 1 x 1 arcsec
    'CRPIX1': 1, 
    'CRVAL1': 337.5202808, 
    'NAXIS1': 1024,
    'CTYPE2': 'DEC--TAN', # The second axis is DEC, declination
    'CUNIT2': 'deg', 
    'CDELT2': 0.0002777777778, 
    'CRPIX2': 1, 
    'CRVAL2': -20.833333059999998, 
    'NAXIS2': 1024
}
wcs_helix_dict = WCS(wcs_input_dict)
print(wcs_helix_dict)

# We may also create an empty WCS object and fill it later

wcs_helix_list = WCS(naxis=2)
wcs_helix_list.wcs.crpix = [1, 1]
wcs_helix_list.wcs.crval = [337.5202808, -20.833333059999998]
wcs_helix_list.wcs.cunit = ["deg", "deg"]
wcs_helix_list.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs_helix_list.wcs.cdelt = [-0.0002777777778, 0.0002777777778]
wcs_helix_list.array_shape = [1024, 1024]
print(wcs_helix_list)

##### Part 7: Working with Celestial Coordinates - Working with a real FITS file

# Retrieving practice data, originally accessed from Digitized Sky Survey, ESO

hdu_list_pr = fits.open('https://github.com/astropy/astropy-data/raw/6d92878d18e970ce6497b70a9253f65c925978bf/tutorials/celestial-coords1/tailored_dss.22.29.38.50-20.50.13_60arcmin.fits')

print(hdu_list_pr.info()) # header info
image = hdu_list_pr[0].data
header = hdu_list_pr[0].header

print(header)

# Reading WCS from the header
wcs_helix = WCS(header)
print(wcs_helix)

# Plotting the image with pixel numbers in axes

fig = plt.figure(figsize=(10, 10))
plt.imshow(image, origin='lower', cmap='cividis')
plt.show() # Do not forget to put this if you use command panel like me, unlike Jupyter notebook
# Now plotting it with RA - DEC in axes

fig = plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=wcs_helix)
plt.imshow(image, origin='lower', cmap='cividis', aspect='equal')
plt.xlabel(r'RA')
plt.ylabel(r'Dec')

## Overlay in International Celestial Reference System
# overlay = ax.get_coords_overlay('icrs') 
# overlay.grid(color='white', ls='dotted')

## Overlay in Galactic Coordinates
# overlay = ax.get_coords_overlay('galactic')
# overlay.grid(color='white', ls='dotted')

plt.show()

##### Part 7: Working with Celestial Coordinates - Placing a Scale Marker

fig = plt.figure(figsize=(10, 10), frameon=False) # removing frame around fig
ax = plt.subplot(projection=wcs_helix)
# arrow method in Axes.arrow
ax.arrow(337, -21.2, 0, 0.1, 
         head_width=0, head_length=0, 
         fc='red', ec='red', width=0.003, 
         transform=ax.get_transform('icrs'))
plt.text(337.05, -21.18, '0.1 deg', 
         color='red', rotation=90, 
         transform=ax.get_transform('icrs'))
plt.imshow(image, origin='lower', cmap='cividis', aspect='equal')
plt.xlabel(r'RA')
plt.ylabel(r'Dec')

plt.show()

# Now a vertical arrow in white

fig = plt.figure(figsize=(10, 10), frameon=False) # removing frame around fig
ax = plt.subplot(projection=wcs_helix)
# arrow method in Axes.arrow
ax.arrow(337, -21.2, 0.1, 0, 
         head_width=0, head_length=0, 
         fc='white', ec='white', width=0.003, 
         transform=ax.get_transform('icrs'))
plt.text(337.05, -21.18, '0.1 deg', 
         color='white', rotation=0, 
         transform=ax.get_transform('icrs'))
plt.imshow(image, origin='lower', cmap='cividis', aspect='equal')
plt.xlabel(r'RA')
plt.ylabel(r'Dec')

plt.show()
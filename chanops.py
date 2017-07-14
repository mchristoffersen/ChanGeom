# Michael Christoffersen
# python 2.7

# To Do
# Categorize non-data pixels as NaN, 9999, etc. to reduce file size

# File Summary
# This file performs two main operations on a geotiff of a river channel, it measures and
# records the width at each point along the channel and measures the total length of the
# channel

# Input Arguments
#
# These are all the same as the original MATLAB script by Fisher et al., and as such are copied
# directly from it
#
# inputtif = this is the binary rasterized channel polygon geotiff prepared elsewhere. Must have
# background values of 0 and channel area equal to 1 for the algorithm to work. All projection
# information will be maintained.
#
# cellsize = raster cell dimensions (usually in meters)
#
# exporttif = name of the exported channel width data geotiff. This file will have width values
# associated with each centerline point and the same projection information as the inputtif.
#
# NoData = the value of no data (NaN) values that you want when exporting the tif files
# NOTE: I usually use -9999 because I have found ArcGIS to have problems
# with NaN values and -9999 is the equivalent of NoData in ArcGIS. NOTE any
# 0 or greater values will cause actual values to be lost (ie if you use 0 the
# first point will be lost when you SetNull in ArcGIS).
#
# start_pt = stream segment beginning where you want the algorithm to begin
# calculating distance from (i.e. 0 meters) for each pixel along the computed centerline
# It is either the southern tail (i.e. the start point is more southern than the end point)
# (start_pt = 0) or the northern tail (i.e. the start point is more northern than the end point)
# (start_pt = 1) (This was previously a separate function called
# centerlinedist.m)

# Usage

from __future__ import print_function
from osgeo import gdal
from gdalconst import *
import numpy as np
#import skimage.morphology
import time, sys
from PIL import Image
from scipy.ndimage.morphology import distance_transform_edt
from chanops_funcs import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.misc

# Debugging flag
debug = False

#Start timer
starttime = time.time()

#Command line arguments
inputtif = sys.argv[1]
cellsize = float(sys.argv[2])
exporttif = sys.argv[3]
NoData = sys.argv[4]
start_pt = sys.argv[5]

# Open the source raster
src_raster = gdal.Open(inputtif, GA_ReadOnly)

# Print some relevant information
print('Driver: ',src_raster.GetDriver().LongName)
print('Raster Sixe: ',src_raster.RasterXSize,'x',src_raster.RasterYSize)
print('Origin: ',src_raster.GetGeoTransform()[0],',',src_raster.GetGeoTransform()[3])
print('Pixel Size: ',src_raster.GetGeoTransform()[1],',',src_raster.GetGeoTransform()[5])

# Copy the source raster to a numpy array and save its size
binar = np.array(src_raster.GetRasterBand(1).ReadAsArray())
(m,n) = binar.shape

# Index the channel (1) the not channel (0) and save the number of pixels in the channel
idx1 = np.where(binar==1)
idx0 = np.where(binar==0)
len1 = len(idx1[0])

# distance transform midline testing
#distxfrm = distance_transform_edt(binar)
#distxfrm = distxfrm/distxfrm.max()
#grad = np.gradient(distxfrm)
#grad = grad[0]/grad[0].max()/2 + grad[1]/grad[1].max()/2
#scipy.misc.imsave("outfile.jpg",distxfrm*255)
#scipy.misc.imsave("outfile_gradient.jpg",grad*255)

#sys.exit()
# Create channel skeleton and save a copy to prune spurs on
#binar_skel = medial_axis(binar)
binar_skel = thin(binar,0)

#scipy.misc.imsave(exporttif + 'outfile.png', binar_skel.astype("int8")*255)
dataset = gdal.GetDriverByName('GTiff').Create(exporttif + '_skel.tif',n,m,1,GDT_Byte)
dataset.SetGeoTransform(src_raster.GetGeoTransform())
dataset.SetProjection(src_raster.GetProjection())
dataset.GetRasterBand(1).WriteArray(binar_skel*255)

# Get rid of a problematic pattern
probpat2 = np.array([[-1,1,-1],[1,0,1],[-1,1,-1]])
idx_check = np.where(binar_skel==1)
patar = binhitmiss_shift(binar_skel,probpat2,idx_check,1)
patloc = np.where(patar==1)
for idx in range(len(patloc[0])):
	i = patloc[0][idx]
	j = patloc[1][idx]
	binar_skel[i][j+1] = 1

# Re-skeletonize the channel to eliminate the pattern
struct = np.array([[0,1,0],[1,1,1],[0,1,0]])
#binar_skel = medial_axis(binar_skel)

#binar_skel = skmp.skeletonize(binar)
binar_noends = np.zeros([m,n]) + binar_skel

if(debug):
	#output_chborder = gdal.GetDriverByName('GTiff').Create(exporttif + '_skeleton_orig.tif',n,m,1,gdal.GDT_Float32)
	#output_chborder.GetRasterBand(1).WriteArray(binar_skel)
	plt.imsave(exporttif + '_skeleton_orig.tif', binar_skel, cmap=cm.gray)

# Initially find all of the ends of the channel
idx_check = np.where(binar_noends==1) #check all points on the channel skeleton
ends = find_endpoints(binar_noends,idx_check)
numends = len(ends[0])
print('There are ' + str(numends) + ' ends left')

# Prune the ends of the channel until only two remain (removes spurs)
pruned_idxs = [],[]
while numends > 2:
	idx_check = neighborones(binar_noends, ends)
	for k in range(numends):
		i = ends[0][k]
		j = ends[1][k]
		binar_noends[i][j] = 0
		pruned_idxs[0].append(i)
		pruned_idxs[1].append(j)

	ends = find_endpoints(binar_noends,idx_check)
	numends = len(ends[0])
	print('There are ' + str(numends) + ' ends left')

print('')

plt.imsave(exporttif + '_skeleton_pruned.tif', binar_noends, cmap=cm.gray)

# Set all of the points not in the channel to zero?

# Iteratively reconstruct channel ends
nabe_info = findneighbors(ends,pruned_idxs)
nabe_idxs = nabe_info[0]
nabe_stat = nabe_info[1]
newchk = nabe_info[2]
binar_recon = np.array(binar_noends)
ctr = 0

while nabe_stat==True:
	print('Reconstructing Channel Ends' + '.'*ctr, end='\r')
	for idx in range(len(nabe_idxs[0])):
		i = nabe_idxs[0][idx]
		j = nabe_idxs[1][idx]
		binar_recon[i][j] = 1
	newends = nabe_idxs
	nabe_info = findneighbors(newends,newchk)
	nabe_idxs = nabe_info[0]
	nabe_stat = nabe_info[1]
	newchk = nabe_info[2]
	ctr = ctr+1

print('')

# Check to make sure a fork was not created during channel end reconstruction
idx_check = np.where(binar_recon==1)
ends = find_endpoints(binar_recon,idx_check)
if(len(ends[0]) > 2):
	print('End reconstruction created a fork - reverting to pruned channel')
	binar_recon = np.array(binar_noends)
	idx_check = np.where(binar_recon==1)
	ends = find_endpoints(binar_recon,idx_check)


#Removing a problematic pattern
probpat1 = np.array([[-1,1,-1],[-1,1,1],[-1,-1,-1]])
for k in range(4):
	idx_check = np.where(binar_recon==1)
	patar = binhitmiss(binar_recon,np.rot90(probpat1,k),idx_check)
	patloc = np.where(patar==1)
	for idx in range(len(patloc[0])):
		i = patloc[0][idx]
		j = patloc[1][idx]
		binar_recon[i][j] = 0


#Define southern termination and northern termination
if ends[0][0] < ends[0][1]:
	northend = (ends[0][0],ends[1][0])
	southend = (ends[0][1],ends[1][1])
else:
	southend = (ends[0][0],ends[1][0])
	northend = (ends[0][1],ends[1][1])

print('Northern end at ' + str(northend))
print('Southern end at ' + str(southend))

#Euclidian distance transform
skelones = np.where(binar_recon==1)
distxfrm = distance_transform_edt(binar)

# Create new array to assign channel width values from distance transform to
binar_width = np.array(binar_recon)

# Multiply the channel skeleton by the distance transform to extract
# only the width values along the centerline

binar_width = np.multiply(distxfrm,binar_recon)

# Trace length of channel
#Assign starting point
if start_pt == '1':
	cbeg = northend
elif start_pt == '0':
	cbeg = southend
else:
	print('Invalid start_pt')

distot = 0 # assign initial value to total distance variable

binar_dist = np.array(binar_recon) # array to hold the cumulative distance

# Initialize tracing variables
laspt = [cbeg[0]],[cbeg[1]]
curpt = [cbeg[0]],[cbeg[1]]
nexpt = [cbeg[0]],[cbeg[1]]

binar_dist[curpt[0][0]][curpt[1][0]] = 0 # assign value of zero to initial pixel
binar_dist = binar_dist.astype(float)

if(debug):
	output_chborder = gdal.GetDriverByName('GTiff').Create(exporttif + '_reconstructed.tif',n,m,1,gdal.GDT_Float32)
	output_chborder.GetRasterBand(1).WriteArray(binar_recon)

#Iterate up the channel recording the total distance
for idx in range(len(skelones[0])-1):
	nexpt_info = rivcount(curpt,laspt,binar_recon,distot)
	nexpt = nexpt_info[0]
	distot = nexpt_info[1]
	laspt = curpt
	curpt = nexpt
	if(debug):
		print('curpt= ' + str(curpt) + ' nexpt= ' + str(nexpt))

	binar_dist[curpt[0][0]][curpt[1][0]] = distot
	print('Total distance is ' + str(distot) + ' units', end='\r')

print('')

#Save the channel width tif, and the cumulative distance tif
dataset = gdal.GetDriverByName('GTiff').Create(exporttif + '_centerline.tif',n,m,1,gdal.GDT_Float32)
dataset.SetGeoTransform(src_raster.GetGeoTransform())
dataset.SetProjection(src_raster.GetProjection())
dataset.GetRasterBand(1).WriteArray((cellsize*binar_width).astype("float32"))

dataset = gdal.GetDriverByName('GTiff').Create(exporttif + '_centerline_cumdist.tif',n,m,1,GDT_Float32)
dataset.SetGeoTransform(src_raster.GetGeoTransform())
dataset.SetProjection(src_raster.GetProjection())
dataset.GetRasterBand(1).WriteArray(np.multiply(cellsize,binar_dist).astype("float32"))

#Close the source raster
src_raster = None

#Processing time
print(str(int(time.time() - starttime)) + ' seconds to process')

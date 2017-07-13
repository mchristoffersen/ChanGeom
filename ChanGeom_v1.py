# -*- coding: utf-8 -*-
#Prepare data for ChanGeom

# ChanGeom1
# Script prepares KML files and performs the following steps:
# 1. Take KML and re-project (e.g., into UTM or equal area projection)
# 2. store as polygon shapefile
# 3. Rasterize polygon shapefile
# 4. Run Matlab-based RiverWidth algorithm
# 5. combine outputs into a point Shapefile with two Fields: RWidth_m and CumDist_m
#
# Edit the variables below to match your system and set pixel_size
#
# The script will determine if the output file(s) exists and will not delete these files.
# That is, if a file has been previously processed, you want to delete it if
# you want to re-run the script.
#
# Created May 2014 by Bodo Bookhagen
# Edited September 2016 by Michael Christoffersen


# Import gdal and ogr modules
try:
    from osgeo import gdal, ogr, osr
except ImportError:
    import gdal, ogr, osr
import os, time, subprocess, sys, csv

try:
    import numpy as Numeric
except ImportError:
    import Numeric

###
# PATH variables
IN_KML_PATH="./SampleChannelMasks/kml/"
OUT_SHAPEFILE_PATH="./SampleChannelMasks/out_shp/"
OUT_TIF_PATH="./SampleChannelMasks/tif/"
OUT_CENTERLINE_TIF_PATH="./SampleChannelMasks/center_tif/"
OUT_CENTERLINE_SHAPE_PATH="./SampleChannelMasks/center_shp/"

###
# gdal/ogr command path (in case you have muldiple versions installed)
ogr2ogr_command="/usr/bin/ogr2ogr"
gdal_rasterize_command="/usr/bin/gdal_rasterize"
gdal_translate_command="/usr/bin/gdal_translate"
gdal_polygonize_command='/usr/bin/gdal_polygonize.py'

###
# Matlab
matlab_command='/Applications/MATLAB_R2012a.app/bin/matlab'

###
# Processing parameters
pixel_size = 1 #set pixel size in projection units (m)
LETIBET_PROJ4="+proj=laea +lat_0=35 +lon_0=85 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

NoData = -9999
start_pt = 1

###
# Code starts here - no editing below this line

class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self

def raster2point(srcfile1, srcfile2, dstfile): #Write tif to CSV, ignoring 0s
    # based on gdal2xyz.py
    # reads two tif files and collects data from non-0 values
    srcwin = None
    skip = 1
    delim = ','

    # Open source file.
    srcds1 = gdal.Open( srcfile1 )
    if srcds1 is None:
        print('Could not open %s.' % srcfile1)
        sys.exit( 1 )

    band = srcds1.GetRasterBand(1)
    if band is None:
        print('Could not get band %d' % band_num)
        sys.exit( 1 )
    bands = []
    bands.append(band)
    gt = srcds1.GetGeoTransform()

    srcds2 = gdal.Open( srcfile2 )
    if srcds2 is None:
        print('Could not open %s.' % srcfile2)
        sys.exit( 1 )

    band2 = srcds2.GetRasterBand(1)
    if band2 is None:
        print('Could not get band %d' % band_num)
        sys.exit( 1 )
    bands2 = []
    bands2.append(band)
    gt = srcds2.GetGeoTransform()

    # Collect information on all the source files.
    srcwin1 = (0,0,srcds1.RasterXSize,srcds1.RasterYSize)

    # Open the output file.
    if dstfile is not None:
        dst_fh = open(dstfile,'wt')
    else:
        print "Can not open destination file: " + dstfile
        sys.exit( 1 )

    band_format = (("%g" + delim + "%g" + delim) * len(bands)).rstrip(delim) + '\n'

    # Setup an appropriate print format.
    if abs(gt[0]) < 180 and abs(gt[3]) < 180 \
       and abs(srcds1.RasterXSize * gt[1]) < 180 \
       and abs(srcds1.RasterYSize * gt[5]) < 180:
        format = '%.10g' + delim + '%.10g' + delim + '%s' + delim + '%s'
    else:
        format = '%.3f' + delim + '%.3f' + delim + '%s' + delim + '%s'
    dst_fh.write( 'X,Y,RWidth_m,CumDist_m\n' )

    # Loop emitting data.
    for y in range(srcwin1[1],srcwin1[1]+srcwin1[3],skip):
        data1 = []
        data2 = []
        for band in bands:
            band_data1 = band.ReadAsArray( srcwin1[0], y, srcwin1[2], 1 )
            band_data1 = Numeric.reshape( band_data1, (srcwin1[2],) )
            band_data2 = band2.ReadAsArray( srcwin1[0], y, srcwin1[2], 1 )
            band_data2 = Numeric.reshape( band_data2, (srcwin1[2],) )

            data1.append(band_data1)
            data2.append(band_data2)

        for x_i in range(0,srcwin1[2],skip):
            x = x_i + srcwin1[0]
            if data1[0][x_i] == 0:
                continue
            geo_x = gt[0] + (x+0.5) * gt[1] + (y+0.5) * gt[2]
            geo_y = gt[3] + (x+0.5) * gt[4] + (y+0.5) * gt[5]

            x_i_data1 = []
            x_i_data2 = []
            for i in range(len(bands)):
                x_i_data1.append(data1[i][x_i])
                x_i_data2.append(data2[i][x_i])
            band_str1 = "%g" % tuple(x_i_data1)
            band_str2 = "%g\n" % tuple(x_i_data2)
            line = format % (float(geo_x),float(geo_y), band_str1, band_str2)
            dst_fh.write( line )

def csv2shapefile( csvfile, shapefile, spatialreference ):
    # use a dictionary reader so we can access by field name
    reader = csv.DictReader(open(csvfile,"rb"),
        delimiter=',',
        quoting=csv.QUOTE_NONE)

    # set up the shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(shapefile)
    # create the layer
    layer = data_source.CreateLayer("RWidth", spatialreference, ogr.wkbPoint)

    # Add the fields we're interested in
    layer.CreateField(ogr.FieldDefn("Y", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("X", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("RWidth_m", ogr.OFTReal))
    layer.CreateField(ogr.FieldDefn("CumDist_m", ogr.OFTReal))

    # Process the text file and add the attributes and features to the shapefile
    for row in reader:
      feature = ogr.Feature(layer.GetLayerDefn())
      # Set the attributes using the values from the delimited text file
      feature.SetField("X", row['X'])
      feature.SetField("Y", row['Y'])
      feature.SetField("RWidth_m", row['RWidth_m'])
      feature.SetField("CumDist_m", row['CumDist_m'])

      # create the WKT for the feature using Python string formatting
      wkt = "POINT(%f %f)" %  (float(row['X']) , float(row['Y']))

      # Create the point from the Well Known Txt
      point = ogr.CreateGeometryFromWkt(wkt)

      # Set the feature geometry using the point
      feature.SetGeometry(point)
      # Create the feature in the layer (shapefile)
      layer.CreateFeature(feature)
      # Destroy the feature to free resources
      feature.Destroy()

    # Destroy the data source to free resources
    data_source.Destroy()

t_start = time.time()
for file in os.listdir(IN_KML_PATH):
     if file.endswith('.kml'):
          vector_file = os.path.join(IN_KML_PATH, file)
          print "Converting and Reprojecting " + vector_file

          #project KML into equal-area projection (UTM, LETIBET, LESAMRR)
          #define output projection
          outSpatialRef = osr.SpatialReference()
          #if there is a *.prj file in the PATH directory, take the projection information from that file
          #Can read ESRI and PROJ4 projection information
          for projfile in os.listdir(IN_KML_PATH):
              if projfile.endswith('.prj'):
                  prj_file = os.path.join(IN_KML_PATH, projfile)
                  #import projection file. FIRST try if this is an ESRI prj file
                  prj_text = open(prj_file, 'r').read()
                  try:
                      err = outSpatialRef.ImportFromESRI([prj_text])
                  except err != 0:
                      err = outSpatialRef.ImportFromProj4(prj_text)
                      if err != 0:
                          raise ValueError("Error importing ESRI or PROJ4 projection information from: %s" % prj_file)
                          break
                  if err == 0:
                      #print outSpatialRef.ExportToProj4()
                      print
          # Set Projection name to name of prj file
          outSpatialRef.SetProjCS(prj_file.split('/')[-1][0:-4]);

          #verify if OUT_SHAPEFILE_PATH exists
          if not os.path.exists(OUT_SHAPEFILE_PATH):
              os.makedirs(OUT_SHAPEFILE_PATH)
          if os.path.exists(OUT_SHAPEFILE_PATH) == False:
              print "Can not create directory: " + OUT_SHAPEFILE_PATH

          out_shapefile=OUT_SHAPEFILE_PATH + file.split('.')[0] + "_projected.shp"
          if os.path.exists(out_shapefile):
              print out_shapefile + " exists, skipping to next file..."
          else:
              #call ogr2ogr to convert KML
              out_shapefile_txt = out_shapefile + '.ogr2ogr.out'
              ogr2ogr_subprocess_command = ogr2ogr_command + ' -a_srs EPSG:4326 -t_srs ' + '"' + outSpatialRef.ExportToProj4() + '"' + ' -f "ESRI Shapefile" ' + out_shapefile + ' ' + vector_file + '>' + out_shapefile_txt
              os.system(ogr2ogr_subprocess_command)

# Process: Polygon to Raster
# Load projected shapefile and convert to grid
for file in os.listdir(OUT_SHAPEFILE_PATH):
     if file.endswith('.shp'):
          vector_file = os.path.join(OUT_SHAPEFILE_PATH, file)
          print "Rasterizing " + vector_file

          #verify if OUT_TIF_PATH exists
          if not os.path.exists(OUT_TIF_PATH):
              os.makedirs(OUT_TIF_PATH)
          if os.path.exists(OUT_TIF_PATH) == False:
              print "Can not create directory: " + OUT_TIF_PATH

          out_tiffile=OUT_TIF_PATH + file.split('.')[0] + "_r" + str(pixel_size) + "m.tif"

          if os.path.exists(out_tiffile):
              print out_tiffile + " exists, skipping to next file..."
          else:
              #call gdal_rasterize to rasterize shapefile
              out_tiffile_txt = out_tiffile.split('.')[0] + '.gdal_rasterize.out'
              gdal_rasterize_subprocess_command = gdal_rasterize_command + ' -co "NBITS=1" -burn 1 -tap -a_nodata 0 -ot Byte -tr ' + str(pixel_size) + ' ' + str(pixel_size) + ' ' + vector_file + ' ' + out_tiffile + '>'  + out_tiffile_txt
              os.system(gdal_rasterize_subprocess_command)

#start Matlab script for each TIF file in OUT_TIF_PATH
for file in os.listdir(OUT_TIF_PATH):
     if file.endswith('.tif'):
          in_tiffile = os.path.join(OUT_TIF_PATH, file)
          print "Finding Centerline " + in_tiffile

          #verify if OUT_CENTERLINE_TIF_PATH exists
          if not os.path.exists(OUT_CENTERLINE_TIF_PATH):
              os.makedirs(OUT_CENTERLINE_TIF_PATH)
          if os.path.exists(OUT_CENTERLINE_TIF_PATH) == False:
              print "Can not create directory: " + OUT_CENTERLINE_TIF_PATH

          out_tiffile1=OUT_CENTERLINE_TIF_PATH + file.split('.')[0] + "_centerline.tif"
          out_tiffile2=OUT_CENTERLINE_TIF_PATH + file.split('.')[0] + "_centerline_cumdist.tif"

          if os.path.exists(out_tiffile1) and os.path.exists(out_tiffile2):
              print out_tiffile1.split('/')[-1] + " and " + out_tiffile2.split('/')[-1] + " exist, skipping to next file..."
          else:



              #call Matlab to create centerline and cumulative distance TIFs
              #subprocess.call([matlab_command+" -nosplash -nodisplay -r \"chanextract(\'%s\',\'%s\',0)\"" % (in_tiffile, out_tiffile1)],shell=True);
              subprocess.call(['python','chanops.py',in_tiffile,str(pixel_size),OUT_CENTERLINE_TIF_PATH + file.split('.')[0],str(NoData),str(start_pt)]) #Michael Christoffersen

#post-process channel width from Matlab outputs
for file in os.listdir(OUT_CENTERLINE_TIF_PATH):
     if file.endswith('centerline.tif'):
          in_centerline_tiffile = os.path.join(OUT_CENTERLINE_TIF_PATH, file)
          print "Vectorizing results: " + in_centerline_tiffile
          cumdist_filename = file.split('.')[0] + '_cumdist.tif'
          in_cumdist_tiffile = os.path.join(OUT_CENTERLINE_TIF_PATH, cumdist_filename)

          #verify if OUT_CENTERLINE_SHAPE_PATH exists
          if not os.path.exists(OUT_CENTERLINE_SHAPE_PATH):
              os.makedirs(OUT_CENTERLINE_SHAPE_PATH)
          if os.path.exists(OUT_CENTERLINE_SHAPE_PATH) == False:
              print "Can not create directory: " + OUT_CENTERLINE_SHAPE_PATH

          out_centerline_tiffile=OUT_CENTERLINE_TIF_PATH + file.split('.')[0] + "2.tif"
          out_cumdist_tiffile=OUT_CENTERLINE_TIF_PATH + cumdist_filename.split('.')[0] + "2.tif"

          #csv file contains RWidth_m and CumDist_m
          out_csvfile=OUT_CENTERLINE_TIF_PATH + file.split('.')[0] + ".csv"
          out_shapefile=OUT_CENTERLINE_SHAPE_PATH + file.split('.')[0] + ".shp"

          if os.path.exists(out_shapefile):
              print out_shapefile + " exists, skipping to next file..."
          else:
              #call gdal to set background value (0 from Matlab-tif files)
              out_tiffile_txt = out_centerline_tiffile.split('.')[0] + '.gdal_translate.out'
              gdal_translate_subprocess_command = gdal_translate_command + ' -a_nodata 0 ' + in_centerline_tiffile + ' ' + out_centerline_tiffile + '>'  + out_tiffile_txt
              os.system(gdal_translate_subprocess_command)
              gdal_translate_subprocess_command = gdal_translate_command + ' -a_nodata 0 ' + in_cumdist_tiffile + ' ' + out_cumdist_tiffile + '>>'  + out_tiffile_txt
              os.system(gdal_translate_subprocess_command)
              #delete centerline and cumdist files with NoDATA and rename *2.tif to previous filename
              os.remove(in_centerline_tiffile)
              os.rename(out_centerline_tiffile, in_centerline_tiffile)
              os.remove(in_cumdist_tiffile)
              os.rename(out_cumdist_tiffile, in_cumdist_tiffile)

              #convert tif file to CSV file (use only values > 0)
              raster2point(in_centerline_tiffile, in_cumdist_tiffile, out_csvfile)

              #convert to point shapefile
              csv2shapefile( out_csvfile, out_shapefile, outSpatialRef )

print "Processing time " + str(prettyfloat(time.time() - t_start)) + "s"

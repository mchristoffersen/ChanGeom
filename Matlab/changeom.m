%ChanGeom v0.3
%Updated 2/8/14
    %Combined centerline and channel width algorithms that were separate in v.1
    %Also put an extra filter in to prevent errors when getting the
        %centerline_cumdist that causes it to not process all the way
        %through (see manual for more info)
    
%This is the general workflow to extract channel centerlines and widths
    %from binary tif files that are channel (1) and no channel (0). You can
    %execute the progression of function calls on multiple files in this M-file.
%-------------------------------------------------------------------------------------------------------------
%HELPFUL HINTS  
    %Make sure you have Matlab 2012b or later because the geotiffwrite
        %function on earlier versions has a limit of exporting x or y
        %dimensions of <20,000. This is not an issue starting with Matlab 2012B
    %If you are working with big files make sure you go to
        %Preferences->General->MAT-Files and switch to -v7.3 so you can
        %save big files. Otherwise it won't write the rwidth variable to
        %the MAT file created
%-------------------------------------------------------------------------------------------------------------   
%PART 1 - (OPTIONAL) Convert tiff file to a logical if it is a big file. Just saves 
        %space on your hard drive by shrinking the inputtif
        
%Just put the directory where all the tif files are located and then run the following

myFolder = '/Users/burch/Desktop/WidthData/Processing';
filePattern = fullfile(myFolder, '*.tif');
tifFiles = dir(filePattern);
for k = 1:length(tifFiles)
    baseFileName = tifFiles(k).name;
    fullFileName = fullfile(myFolder, baseFileName);
    [data,header] = geotiffread(fullFileName);
    info = geotiffinfo(fullFileName);
    datal = logical(data);
    %EXPORT Geotiff of Rwidth Result
    key=info.GeoTIFFTags.GeoKeyDirectoryTag;
    geotiffwrite(baseFileName,datal, header, 'GeoKeyDirectoryTag', key);
end
clear all
%------------------------------------------------------------------------------------------------------------- 
%PART 2 - CENTERLINE AND WIDTH EXTRACTION

chanextract('kg_trib1_let2m.tif',2,'kg_trib1_let2m_rwidth1.tif',-9999,0);clear all
chanextract('kg_trib2_let2m.tif',2,'kg_trib2_let2m_rwidth1.tif',-9999,0);clear all
chanextract('mars_trib1_let2m.tif',2,'mars_trib1_let2m_rwidth1.tif',-9999,0);clear all

%FUNCTION--> chanextract(inputtif,cellsize,exporttif,NoData,start_pt) 

%VARIABLES
%inputtif = this is the binary rasterized channel polygon geotiff prepared elsewhere. Must have 
    %background values of 0 and channel area equal to 1 for the algorithm to work. All projection 
    %information will be maintained.
%cellsize = raster cell dimensions (usually in meters)
%exporttif = name of the exported channel width data geotiff. This file will have width values 
    %associated with each centerline point and the same projection information as the inputtif. 
%NoData = the value of no data (NaN) values that you want when exporting the tif files
        %NOTE: I usually use -9999 because I have found ArcGIS to have problems
        %with NaN values and -9999 is the equivalent of NoData in ArcGIS. NOTE any 
        %0 or greater values will cause actual values to be lost (ie if you use 0 the 
        %first point will be lost when you SetNull in ArcGIS).
%start_pt = stream segment beginning where you want the algorithm to begin
        %calculating distance from (i.e. 0 meters) for each pixel along the computed centerline
        %It is either the southern tail (i.e. the start point is more southern than the end point)
        %(start_pt = 0) or the northern tail (i.e. the start point is more northern than the end point)
        %(start_pt = 1) (This was previously a separate function called
        %centerlinedist.m)

%IMPORTANT!!!!
    %Must have the following M-files in the same directory.
        %endpoints.m
        %endpoints_fcn.m
        %riv_msk_fil.m
        %riv_msk.m

%The end result should be a MAT file and an image plot to verify that it
%worked. Also geotiff files with channel width and cumulative upstream distance along
%the centerline will be exported

%ISSUES - the only issue I have noticed with this algorithm is if there are
    %places where three pixels are adjacent instead of two (which can occur based on the
    %thinning algorithm in chanextract.m) the algorithm will enter an endless loop and stop at that spot.
    %You will need to go in and find the extra pixel and delete its value 
    %and then rerun the algorithm. This is detailed further in the pdf manual.  
    
%VARIABLES
    %rwidth1 contains the channel widths
    %centerline_numbers contains upward counting numbers for each pixel along the centerline
    %centerline_distance contains distances between each pixel
    %centerline_cumdist contains distances from first pixel along the centerlne

%------------------------------------------------------------------------------------------------------------- 
%FINISHING UP
%Just check and make sure stuff looks ok
load('kg_trib2_let2m_rwidth1.mat')    

%These are a check to make sure it processed all the way through - if you 
    %see a -9999 for any of the cumdist points there was a mistake.
rwidth1(south_point(1,1), south_point(1,2))
rwidth1(north_point(1,1), north_point(1,2))
centerline_cumdist(south_point(1,1), south_point(1,2))
centerline_cumdist(north_point(1,1), north_point(1,2)) 

figure; imagesc(rwidth1)
figure; imagesc(centerline_cumdist); colorbar

%-------------------------------------------------------------------------------------------------------------       
%TO FIX A PROBLEM IN PART 2 - FIND THE PROBLEM POINT AND THEN RUN THE FOLLOWING
%see manual for example
%This has largely been fixed in v0.3

load('seti_trib1_let2m_rwidth1.mat')
y=1715; x=4400;
imagesc(centerline((y-1):(y+1),(x-1):(x+1)));
rwidth1(y,x)=-9999;
centerline(y,x)=0;
clear x y 
save(MATname)
clear all
%now rerun the second part of chanextract.m manually and 
    %also reexport the rwidth1 variable as a .tif using Part 1 (above)
%EXPORT Geotiff of Rwidth Result
key=info.GeoTIFFTags.GeoKeyDirectoryTag;
geotiffwrite('spit4_let3m_rwidth1.tif',rwidth1, header, 'GeoKeyDirectoryTag', key);
%-------------------------------------------------------------------------------------------------------------
%OPTIONAL - This script combines the riv_msk and distance_qe files and
    %gives you the real halfwidth distance for each pixel in the channel 
    %for use in figures. Null value is -9999 upon export.

%msk = single(riv_msk);
%dist = single(riv_msk_distance_qe);
%idf = find(dist == 0);
%dist(idf)=-1;
%dist=dist+1;
%rwidth_algorithm_pt_dist=msk+dist-1;
%idf = find(rwidth_algorithm_pt_dist == -1);
%rwidth_algorithm_pt_dist(idf) = -9999; 
%rwidth_algorithm_m_dist = rwidth_algorithm_pt_dist.*gridcell_size;
%id0 = find(rwidth_algorithm_m_dist==0);
%rwidth_algorithm_m_dist(id0)=0.1;
%rwidth_algorithm_m_dist(idf) = -9999; 
%clear msk dist idf
    
  




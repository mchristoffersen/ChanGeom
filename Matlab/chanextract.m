%ChanExtract
    %The program reads TIF files and calculates a centerline and the
    %across-channel river width as given by the polygon boundary.

function chanextract(inputtif,cellsize,exporttif,NoData,start_pt) 
%Example of function call --> %chanextract('coloradoriver.tif',3,'coloradoriver_rwidth1.tif',NoData)

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

%----------------------------------------------------------------------------------------
%IF MANUALLY RUNNING ENTER DATA--->
%Bring geotiff file with polygon fill equal to 1 and background equal to 0
[riv_msk_fil,header] = geotiffread(inputtif);
info = geotiffinfo(inputtif);
riv_msk_cellsize = cellsize
gridcell_size= cellsize
%--------->HERE
%---------------------------------------------------------------------------------------- 

%Need to make outline polygon from tiff file
riv_msk = bwmorph(riv_msk_fil, 'remove');
[idx1, idx0, nr1, riv_msk_fil] = riv_mask_fil(riv_msk);

%Create centerline
fprintf(1,'filtering binary grid...\n');
centerline = bwmorph(riv_msk_fil, 'thin', 500);
centerline = bwmorph(centerline, 'clean', 10);

%Prune with endpoints until only two endpoints exist (start & beginning)
no_end_points = centerline;
nr_of_endpoints = 3; 
nr_of_iterations = 0;

while(nr_of_endpoints > 2)
    nr_of_iterations = nr_of_iterations + 1; 
    end_points = endpoints(no_end_points); 
    no_end_points = no_end_points & ~end_points; 
    nr_of_endpoints = length(find(end_points > 0));
    if nr_of_iterations > 500, fprintf(1,'more than 500 iterations in finding endpoints - returning\n'); return; end
end

centerline = no_end_points; 
clear no_end_points clear nr_of_endpoints nr_of_iterations

%Alternative approach to prune spurs, but slow:
%centerline = bwmorph(centerline, 'spur', 100); %this takes a long time and may have to be altered - it removes start and endpoints, too

centerline = bwmorph(centerline, 'thin');

%Do This
idx_centerline = find(centerline == 1); idx_no_centerline = find(centerline == 0);

%find southernmost point
centerline_xy = size(centerline); 
end_points = endpoints(centerline); [endpoint_idxx, endpoint_idxy] = find(end_points > 0);
%first coorindate in centerline_xy is Y extend - find coordinate in
%endpoint_idxx that is closest to centerline_xy(1)
idx_south = find(min(centerline_xy(1) - endpoint_idxx) == (centerline_xy(1) - endpoint_idxx));
idx_north = find(max(centerline_xy(1) - endpoint_idxx) == (centerline_xy(1) - endpoint_idxx));
south_point = [endpoint_idxx(idx_south) endpoint_idxy(idx_south)];
north_point = [endpoint_idxx(idx_north) endpoint_idxy(idx_north)];
centerline_idx = bwtraceboundary(centerline, south_point, 'N', 8);
if isempty(centerline_idx), fprintf(1,'Can not find bwtraceboundary of centerline'); return; end

%centerline_idx no contains the pixel coordinates of the filtered centerline
clear idx_south idx_north end_points endpoint_idxx endpoint_idxy

%Quasi-Euclidean Distance -   
    %calculate distance from edge to centerline and multiplies x 2    
    riv_msk_distance_qe = single(bwdist(riv_msk, 'quasi-euclidean')); 
    riv_msk_distance_qe(idx0) = 0; 
    riv_msk_distance_qe_idx = find(riv_msk_distance_qe == 1);
    rwidth1 = riv_msk_distance_qe; 
    rwidth1(idx_no_centerline) = 0; 
    rwidth1 = single(rwidth1 .* (2 * gridcell_size)+gridcell_size); 
    rwidth1(idx_no_centerline) = 0;
    %clear riv_msk_distance
    %This adds one pixel to the width and then resets the background values to 0
    %save riverwidth values
    rwidth1_col = single(NaN(size(centerline_idx,1),1));
    for i = 1:size(centerline_idx,1), rwidth1_col(i) = rwidth1(centerline_idx(i,1), centerline_idx(i,2)); end
    %rwdith1 contains riverwidth in meters along the centerline from south to north

%------------------------------------------------------------------------------------------------------------------
index0 = find(rwidth1==0);
rwidth1(index0)= NoData;

clear centerline_idx centerline_xy i idx0 idx1 idx_centerline idx_no_centerline index0 nr1...
    riv_msk_distance_qe_idx rwidth1_col gridcell_size riv_msk_cellsize

%EXPORT a .MAT file with all the variables from the function (name will be
    %the same as the exporttif without the tif)
MATname = strtok(exporttif,'.')
save(MATname) 

%View the result
imagesc(rwidth1); colorbar

%EXPORT Geotiff of Rwidth Result
key=info.GeoTIFFTags.GeoKeyDirectoryTag;
geotiffwrite(exporttif,rwidth1, header, 'GeoKeyDirectoryTag', key);
%------------------------------------------------------------------------------------------------------------------  
%------------------------------------------------------------------------------------------------------------------  
%------------------------------------------------------------------------------------------------------------------  
%RUN FROM HERE MANUALLY IF THERE IS AN ISSUE WITH CENTERLINE_CUMDIST VARIABLE

%NOW Find the centerline distances
if start_pt == 0; start = south_point;
else start = north_point;
end    

%-------------------------------------------------------------------------
%INPUT VALUES HERE if doing manually
%To adjust the distance for the appropriate pixel size - currently only
    %operates as 1 and sqrt of 2 so this is a multiplier

%Define first_pixelx and first_pixely, the algorithm goes in any direction
%but you must have the coordinates of the last or first point along the
%centerline
first_pixely= start(1,2); %actually x value
first_pixelx = start(1,1); %actually y value 
%-------------------------------------------------------------------------

centerline_numbers = single(NaN(size(centerline)));
centerline_distance = single(NaN(size(centerline)));
centerline_cumdist = single(NaN(size(centerline)));
centerline_processed = zeros(size(centerline)); centerline_processed = logical(centerline_processed);

runningdistance = 0;
counter = 1;
whilecounter = 1;
stop = 0;

while (stop == 0)
    if (mod(whilecounter, 100) == 0), fprintf('%d, ', whilecounter); end
    if (mod(whilecounter, 2000) == 0), fprintf('\n'); end
    %start at first_pixelx, first_pixely
    if whilecounter == 1
        centerline_numbers(first_pixelx, first_pixely) = counter; centerline_distance(first_pixelx, first_pixely) = 0;
        centerline_cumdist(first_pixelx, first_pixely) = 0; centerline_processed(first_pixelx, first_pixely) = 1;
        if ((centerline(first_pixelx,first_pixely-1) == 1) & (centerline_processed(first_pixelx,first_pixely-1) == 0)), next_x = first_pixelx; next_y = first_pixely-1; step_y = 1; step_x = 0;
        elseif ((centerline(first_pixelx,first_pixely+1) == 1) & (centerline_processed(first_pixelx,first_pixely+1) == 0)), next_x = first_pixelx; next_y = first_pixely+1; step_y = 1; step_x = 0;
        elseif ((centerline(first_pixelx-1,first_pixely-1) == 1) & (centerline_processed(first_pixelx-1,first_pixely-1) == 0)), next_x = first_pixelx-1; next_y = first_pixely-1; step_y = 1; step_x = 1;
        elseif ((centerline(first_pixelx-1,first_pixely) == 1) & (centerline_processed(first_pixelx-1,first_pixely) == 0)), next_x = first_pixelx-1; next_y = first_pixely; step_y = 0; step_x = 1;
        elseif ((centerline(first_pixelx-1,first_pixely+1) == 1) & (centerline_processed(first_pixelx-1,first_pixely+1) == 0)), next_x = first_pixelx-1; next_y = first_pixely+1; step_y = 1; step_x = 1;
        elseif ((centerline(first_pixelx+1,first_pixely-1) == 1) & (centerline_processed(first_pixelx+1,first_pixely-1) == 0)), next_x = first_pixelx+1; next_y = first_pixely-1; step_y = 1; step_x = 1;
        elseif ((centerline(first_pixelx+1,first_pixely) == 1) & (centerline_processed(first_pixelx+1,first_pixely) == 0)), next_x = first_pixelx+1; next_y = first_pixely; step_y = 0; step_x = 1;
        elseif ((centerline(first_pixelx+1,first_pixely+1) == 1) & (centerline_processed(first_pixelx+1,first_pixely+1) == 0)), next_x = first_pixelx+1; next_y = first_pixely+1; step_y = 1; step_x = 1;
        else stop = 1; %no adjacent pixel found
        end
        distance = sqrt(step_x^2 + step_y^2);
        runningdistance = runningdistance + distance;
        centerline_distance(next_x,next_y) = distance;
        centerline_cumdist(next_x,next_y) = runningdistance;
        counter = counter + 1;
        centerline_processed(next_x, next_y) = 1;
    else
        %now find adjacent pixel
        if ((centerline(next_x,next_y-1) == 1) & (centerline_processed(next_x,next_y-1) == 0)), next_x = next_x; next_y = next_y-1; step_y = 1; step_x = 0;
        elseif ((centerline(next_x,next_y+1) == 1) & (centerline_processed(next_x,next_y+1) == 0)), next_x = next_x; next_y = next_y+1; step_y = 1; step_x = 0;
        elseif ((centerline(next_x-1,next_y-1) == 1) & (centerline_processed(next_x-1,next_y-1) == 0)), next_x = next_x-1; next_y = next_y-1; step_y = 1; step_x = 1;
        elseif ((centerline(next_x-1,next_y) == 1) & (centerline_processed(next_x-1,next_y) == 0)), next_x = next_x-1; next_y = next_y; step_y = 0; step_x = 1;
        elseif ((centerline(next_x-1,next_y+1) == 1) & (centerline_processed(next_x-1,next_y+1) == 0)), next_x = next_x-1; next_y = next_y+1; step_y = 1; step_x = 1;
        elseif ((centerline(next_x+1,next_y-1) == 1) & (centerline_processed(next_x+1,next_y-1) == 0)), next_x = next_x+1; next_y = next_y-1; step_y = 1; step_x = 1;
        elseif ((centerline(next_x+1,next_y) == 1) & (centerline_processed(next_x+1,next_y) == 0)), next_x = next_x+1; next_y = next_y; step_y = 0; step_x = 1;
        elseif ((centerline(next_x+1,next_y+1) == 1) & (centerline_processed(next_x+1,next_y+1) == 0)), next_x = next_x+1; next_y = next_y+1; step_y = 1; step_x = 1;
        else stop = 1; %no adjacent pixel found
        end
        
        %centerline_numbers(next_x,next_y) =
        %calculate distance
        distance = sqrt(step_x^2 + step_y^2);
        runningdistance = runningdistance + distance;
        centerline_distance(next_x,next_y) = distance;
        centerline_cumdist(next_x,next_y) = runningdistance;
        centerline_processed(next_x, next_y) = 1;
        centerline_numbers(next_x, next_y) = counter;
        counter = counter + 1;
    end
    whilecounter = whilecounter + 1;
end %while
fprintf('\n')

%Multiplies results times cellsize
centerline_distance = centerline_distance.*cellsize;
centerline_cumdist = centerline_cumdist.*cellsize;

%Converts NaN values in background to whatever the NoData value is
    %specified as in the function
poop = isnan(centerline_cumdist);
centerline_cumdist(poop)= NoData;
centerline_numbers(poop)= NoData;

%Images the centerline_cumdist to make sure it didn't get truncated
imagesc(centerline_cumdist); colorbar

%Uses the name of the .mat file and concatenates with the numbers and
    %cumdist files for export
name = strtok(MATname,'.')
cent_numb = strcat(name,'_cent_numbers.tif')
cent_cumdist = strcat(name,'_cent_cumdist.tif')


%EXPORT Geotiffs of Results
key=info.GeoTIFFTags.GeoKeyDirectoryTag;
%geotiffwrite(cent_numb,centerline_numbers, header, 'GeoKeyDirectoryTag', key);
geotiffwrite(cent_cumdist,centerline_cumdist, header, 'GeoKeyDirectoryTag', key);

clear counter distance next_x next_y poop runningdistance step_x step_y stop...
    whilecounter x y centerline_distance first_pixelx first_pixely centerline_processed...
    chanextractMATfile exporttif name cent_numb cent_cumdist start
save(MATname)

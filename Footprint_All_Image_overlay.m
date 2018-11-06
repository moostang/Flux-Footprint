% Plot Flux footprint
% -------------------------
% Author: Gyanendra Gurung
% Date: 2018/09/11
% 
% DIRECTIONS FOR USE
% ------------------
% 1. Input the directory location where the satellite images are stored,
%    i.e. the locatin of    TCI_CROPPED   directory
% 2. Select directory of your choice for OUTPUT
% 3. Click on satellite image for location of Tower
% 4. Select the "filename part" that represents the DATE for image
%    acquisition.
% 5. DONE
% https://gis.stackexchange.com/questions/241183/sentinel-2-reflectance-calculation
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% All vectors need to be of equal length (one value for each time step)
% 
%    zm       = Measurement height above displacement height (i.e. z-d) [m]
%               usually a scalar, but can also be a vector 
%    z0       = Roughness length [m] - enter [NaN] if not known 
%               usually a scalar, but can also be a vector 
%    umean    = Vector of mean wind speed at zm [ms-1] - enter [NaN] if not known 
%               Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint
%    h        = Vector of boundary layer height [m]
%    ol       = Vector of Obukhov length [m]
%    sigmav   = Vector of standard deviation of lateral velocity fluctuations [ms-1]
%    ustar    = Vector of friction velocity [ms-1]
%    wind_dir = Vector of wind direction in degrees (of 360) for rotation of the footprint    
%    BoundingBox = Minimum and maximumm X and Y
%    X           = Coordinates of vector along x-axis
%    Y           = Coordinates of vector along y-axis
%    ID          = ID of the vector
%    CONTOUR     = Value of the vector
% 
% Station Coordinate Information
% ------------------------------
% Use UTM values as station's x and y coordinates when exporting shapefiles
% for use use in GIS software. Since, there have been problems reading
% geo-reference GeoTiff files in Matlab, use number of pixel counts to the
% station's location as x and y coordinate. 
% 
%    station_x = x-coordinate value station location in UTM or pixel
%    station_y = y-coordinate value station location in UTM or pixel
% 
% For exporting to shapefiles, please use 
% Coordinate System: WGS 1984 UTM ZONE 13N
%    station_x = 466519; station_y = 7652699;
% 
% This is for Coordinate System: WGS 1984 UTM ZONE 13N
% 
% Otherwise, for the current image, the station is located at 
%    station_x = 4096; station_y = 4318;
% 
% or as selected by user. 
% 
% 
% NOTE: PIXEL COORDINATES WILL VARY AMONG DIFFERENTLY-CROPPED IMAGES
% 
% -------------------------------------------------------------------------
clc; clear; close all
% 
DirTCI = 'C:/GURUNG/MGIS-FINAL-PROJECT-WORK/DATA/IMAGERY/SENTINEL/TCI_CROPPED';
DirSAV = 'output';
station_x = 4107; station_y = 3903;
tHH = 19; 
tMM = 20;
tAcqui = tHH + (tMM/60)
iselect = 3;

% ----------------------------------------------------
% Get Image Location and Station's Coordinates (pixel)
% ----------------------------------------------------
% Select Directory with Satellite Image
if isempty(DirTCI) == 1
    disp('Please select directory with Image files')
    CurDir = pwd;
    DirTCI = uigetdir(CurDir);
end
disp(['Image Directory: ' DirTCI])
    
if isempty(DirSAV) == 1
    disp('Please select directory to save OUTPUT')
    CurDir = pwd;
    DirSAV = uigetdir(CurDir);
end
disp(['OUTPUT Directory: ' DirSAV])

% Input Image Type, e.g. TIFF, jp2, jog, png, etc. 
jp2files = fullfile(DirTCI, '*.jp2');
files = dir(jp2files);

% Show Image to pick station location
if ( isempty(station_x) == 1 ) || ( isempty(station_y) )
    tmpImage = imread([files(1).folder '\' files(1).name]);
    imshow(flip(tmpImage))
    set(gca,'Ydir','normal')
    disp('Please CLICK on the location of the TOWER')
    disp('Does not have to be precise; only for viewing purposes')
    text(0,-120,'Please CLICK on the location of the TOWER', ...
        'Color','red','FontSize',14)
    text(0,60,'(Does not have to be precise; only for viewing purposes', ...
        'Color','red','FontSize',10)
    
% Get Station Coordinates (pixel) from Image
    [station_x, station_y] = ginput(1);
end
disp(['Station coordinates are x = ' num2str(station_x) ' y = ' ...
    num2str(station_y)]);

% Get Dates from Image files
nSplit = split(files(1).name, '_');
disp('Please select the filename [part] that represents the')
disp('acquisition DATE of the Satellite Image')
disp(' ')
for i = 1:length(nSplit)
    tmp = strcat('[',num2str(i),'] ', nSplit(i));
    disp(tmp)
end
disp('     [0] EXIT PROGRAM')

% Get user input
indDate = 7;
while indDate > length(nSplit) || indDate < 0
    prompt = 'Please choose one [Default is 3] = ';
    indDate = input(prompt);
    if(indDate == 0)
        disp('Program Exit')
        break
    end
end
datePart = nSplit{indDate};
disp(strcat('DATE parameter is stored in :', datePart))
% 
% ---------------------------------
% Input plotting range (in meters)
d = [-2000 -2000 2000 2000];
% Bounding box of Map
distFromStation = 2500;
x_min = station_x - distFromStation; x_max = station_x + distFromStation;
y_min = station_y - distFromStation; y_max = station_y + distFromStation;
% 
% Input Boundary Layer Height (for 24 hours and in meters)
blHeight = 300;
% True Wind Direction of Tower
true_wind = 330;
% ---------------------------------
% 
% 
% ----------------------------------------------- 
% Read input MAT file and prepare input variables
% ----------------------------------------------- 
load('cambay.mat')
z0All     = z0';
variance  = wind_stats(5,:)';
sigmavAll = sqrt(variance); clear variance;
wDirAll   = wind_direction'; 
ustarAll  = ustar';
olAll     = L';
hAll      = blHeight*ones(1,length(olAll));
hAll      = hAll';
zm        = sonic_height;
umean     = NaN;
fco2All   = fco2';
hlatAll   = hlatent';
hsensiAll = hsensible';
speedAll  = wind_speed';
% 
% -------------------------------
% Reorientation of Wind Direction
% -------------------------------
% The input MATfile has wind directions for which the true North direction
% is along 330 degrees. This needs to be corrected. This is done by (a)
% first % changing the wind direction to 360 format and reiorienting the 
% (a) clockwise positive values from 330 to 150 degrees, and (c) 
% anticlockwise -ve values from 330 to 150 degrees.
% 
% Reorient clockwise positive values from 330 to 150 degrees
wDirSave   = wDirAll;
wd         = wDirAll;
% Add 330 degrees to ONLY positive values
m          = wDirAll    >   0;
wDirAll(m) = wDirAll(m) + 330; 
k          = wDirAll    > 360;
wDirAll(k) = wDirAll(k) - 360;
p          = wDirAll ==   360;
wDirAll(p) =                0;
% Reorient anticlockwise -ve values from 330 to 150 degrees
n          = wDirAll    <   0;
wDirAll(n) = wDirAll(n) + 330;
clear m, clear n, clear k, clear p 
%
% 
% Batch-process footprint plotting process
tic
for i = 1:length(files)
    disp(files(i).name)
end
indDate = 3;
% Collect windspeed for the day
% Loop over images
for ifiles = 1:length(files)
    Isat = imread([files(ifiles).folder '\' files(ifiles).name]);
    nSplit = split(files(ifiles).name, '_');
    datePart = nSplit{indDate};
    disp(datePart)
% 
%   Select start and end times for footprint plots
    sYY = datePart(1:4); sMM = datePart(5:6); sDD = datePart(7:8);
    dnum      = datenum(str2num(sYY),str2num(sMM),str2num(sDD),0,0,0); % Start Date-time
    dend      = datenum(str2num(sYY),str2num(sMM),str2num(sDD),23,40,0); % End   Date-time (12 hours)    
%
    tStart    = find(tmaster==dnum);
    tEnd      = find(tmaster==dend);
%
    disp([num2str(tStart ) ' ' num2str(tEnd)])
    disp(['for date: ' datestr(datenum([sYY sMM sDD  ], 'yyyymmdd'))])
    j = tStart:tEnd;
%     
    z0        = z0All(    j);
    h         = hAll(     j);
    ol        = olAll(    j);
    sigmav    = sigmavAll(j);
    ustar     = ustarAll( j);
    wind_dir  = wDirAll(  j);
    wDir_rad  = wind_dir.*(pi/180); % In Radians for polar histogram only
    fco2      = fco2All(  j);
    hlatent   = hlatAll(  j);
    hsensible = hsensiAll(j);
    wSpeed    = speedAll( j);
    t         = tmaster(  j);
    tj        = uint16(tAcqui*(60/20));
%     
    figure(ifiles)
% 
    subplot(6,3,[1 2 4 5 7 8 10 11 13 14 16 17]); 
    imshow(flip(Isat))
    hold all;
    set(gca,'Ydir','normal')
    
    subplot(6,3,[3 6]);polarhistogram(wDir_rad,18); 
    title('{Wind Direction (24 hrs)}'); set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    
    subplot(6,3,9);plot(t,wSpeed); hold on
    plot([t(tj) t(tj)], [min(wSpeed) max(wSpeed)], 'Color', [1 0 0]);
    set(gca, 'xticklabel', []); ylabel('{\itu} (m/s)'); hold off
    axis tight
    
    subplot(6,3,12);plot(t,ustar); hold on
    plot([t(tj) t(tj)], [min(ustar) max(ustar)], 'Color', [1 0 0]);
    set(gca, 'xticklabel', []); ylabel('{\itu}^*'); hold off
    axis tight       
    
    subplot(6,3,15);    
    plotyy(t,hsensible,t,hlatent); hold on
    a = get(gca, 'YLim');
    plot([t(tj) t(tj)], a , 'Color', [1 0 0]);
    ylabel('{\itH_{sensible}}, {\itH_{latent}} (Red) ');    
    set(gca, 'xticklabel', []); hold off
    axis tight
    
    subplot(6,3,18);plot(t,fco2); hold on      
    plot([t(tj) t(tj)], [min(fco2) max(fco2)], 'Color', [1 0 0]);
    ylabel('{\itF_c}'); 
    xlabel('Time (HH:MM)'); datetick('x','HH:MM'); 
    axis tight        
    
    % ---------------------------------------------- %     
    % Calculate Footprint climatology for entire day
    % ---------------------------------------------- %   
    [FFP,~] = calc_footprint_FFP_climatology(zm,z0,umean,h,ol,...
        sigmav,ustar,wind_dir,'domain',[-5000 5000 -5000 5000],'nx',...
        3000,'r',10:10:90,'smooth_data',1);
    % 
    % Adjust Station Coordinates
    FFP(1).x_2d = FFP(1).x_2d + station_x;
    FFP(1).y_2d = FFP(1).y_2d + station_y;
    
    % surf(FFP(1).x_2d, FFP(1).y_2d, FFP(1).fclim_2d); shading flat; view(2);

    % Assign Subplot for Footprint 
    subplot(6,3,[1 2 4 5 7 8 10 11 13 14 16 17]); 
    disp([num2str(j), ' Plotting Footprint for ', datestr(dnum) ' UTC'])
    hold all
                    
    for i=1:9
        FFP(i).xr = station_x + FFP(i).xr;
        FFP(i).yr = station_y + FFP(i).yr;
        z = FFP(i).fr.*100000.*ones(size(FFP(i).yr));
        plot3(FFP(i).xr,FFP(i).yr,z,'green'); hold on  

        %
        S(i).Geometry = 'Line';
%         minX = min(FFP(i).xr); minY = min(FFP(i).yr);
%         maxX = max(FFP(i).xr); maxY = max(FFP(i).yr);
%         S(i).BoundingBox = [minX, minY; maxX, maxY];
        S(i).BoundingBox = [x_min, y_min; x_max, y_max];            
        if(isnan(FFP(i).xr(1)))
            S(i).X = fliplr(FFP(i).xr);
            S(i).Y = fliplr(FFP(i).yr);
        end
        S(i).ID = i;
        S(i).CONTOUR = i*10;    
    end
%
    xlabel('Easting (m)')
    ylabel('Northing (m)')

    % --------------------------------------------------------------- %     
    % Footprint for Obs. time corresponding to Image Acquisition time
    % --------------------------------------------------------------- % 
    tmp = datenum(str2num(sYY),str2num(sMM),str2num(sDD),19,00,0); 
    j   = find(tmaster == tmp);
    searchFlg = 1;
    
    disp(["Calculating Footprint for " datestr(tmaster(j))])
    z0        = z0All(    j);
    h         = hAll(     j);
    ol        = olAll(    j);
    sigmav    = sigmavAll(j);
    ustar     = ustarAll( j);
    wind_dir  = wDirAll(  j);
    
    [FFP,~] = calc_footprint_FFP_climatology(zm,z0,...
        umean,h,ol,sigmav,ustar,wind_dir,'domain',...
        [-5000 5000 -5000 5000],'nx',3000,'r',10:10:90,...
        'smooth_data',1);
    
    while isempty(FFP(1).xr) == 1
        if j > tEnd
            searchFlg = 0;
        end
        
        if searchFlg == 1 && j <= tEnd
            j = j + 1;
        elseif searchFlg == 0 && j > tStart
            j = j - 1;
        end
        
        if j < tStart
            disp(['There are not calculable Footprint for date ' datestr(tmaster(j))])
            disp('Please choose another day to calculate footprint')
        end
            
        disp(['Footprint not found for ' datestr(tmaster(j-1))])
        disp(['Calculating next Footprint for ' datestr(tmaster(j))])
        s1 = datestr(dnum, 30);    
        % 
        z0        = z0All(    j);
        h         = hAll(     j);
        ol        = olAll(    j);
        sigmav    = sigmavAll(j);
        ustar     = ustarAll( j);
        wind_dir  = wDirAll(  j);
        
        [FFP,~] = calc_footprint_FFP_climatology(zm,z0,...
                umean,h,ol,sigmav,ustar,wind_dir,'domain',...
                [-5000 5000 -5000 5000],'nx',3000,'r',10:10:90,...
                'smooth_data',1);
    end
    disp(['Footprint found for ' datestr(tmaster(j))])
    subplot(6,3,[1 2 4 5 7 8 10 11 13 14 16 17]); 
    hold all
    for i=1:9
        FFP(i).xr = station_x + FFP(i).xr;
        FFP(i).yr = station_y + FFP(i).yr;
        z = FFP(i).fr.*100000.*ones(size(FFP(i).yr));
        plot3(FFP(i).xr,FFP(i).yr,z,'red'); hold on  

        %
    %                       minX = min(FFP(i).xr); minY = min(FFP(i).yr);
    %                       maxX = max(FFP(i).xr); maxY = max(FFP(i).yr);
        S(i).Geometry = 'Line';
    %                       S(i).BoundingBox = [minX, minY; maxX, maxY];
        S(i).BoundingBox = [x_min, y_min; x_max, y_max];            
        if(isnan(FFP(i).xr(1)))
            S(i).X = fliplr(FFP(i).xr);
            S(i).Y = fliplr(FFP(i).yr);
        end
        S(i).ID = i;
        S(i).CONTOUR = i*10;    
    end
    
    title(['Footprint (UTC): ' datestr(tmaster(j)) ' (Red)'],'FontSize', 8)
    %    Insert raster outputs from ArcGIS. 

%    Unsupervised classification in matlab using k-means (isodata)
%    clustering technique.    
        
    set(gcf, 'units', 'pixels', 'position', [10, 10, 1920, 1080])
    subplot(6,3,18);datetick('x','HH:MM'); axis tight 
    s1 = datestr(dnum, 30);
    saveas(  gcf,[char(DirSAV) '/footprint_overlay' s1 '.png'])
    close all
    clear FFP, clear S, clear z0, clear h, clear ol, clear sigmav
    clear ustar, clear wind_dir
end     

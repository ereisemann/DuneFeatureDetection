% ---- From Lauren Dunkin 02/2017
% ---- edited by Eve Eisemann 03/2017
% ---- edited by Michael Hartman 04/2017
% ---- edited by Eve Eisemann 04/2017, 05/2017
clear,clc

addpath(genpath('E:\Eisemann\Projects\SWG_Feature_Extraction\data\'))

%%%%%%%%%%%%%%%%%  USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%
% --- Adding folder and subfolders to search path --- %

% ---- Save a file? y/n ----%
savefile = 'y';

% --- Change this to skip transects - meters
pSkip = 10;

% ---- Input Smoothing Parameter ---- %
S = 10;
% ---- S = 5-10 helps to minimize jaggedness 

% ---- Input Minimum Crest Prominence ----%
P = 0.2;
% ---- in vertical units, minimum peak prominance. 
% ---- P = 0.2 is good for detecting small foredune but avoiding berm &
% ---- shoreline. Increasing S beyond 5 or 10 will not improve prominance of peaks 

% --- location of scripts etc. 
% --- indir: Source of Grids --- %
indir = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\';
SLbase = '2009_ncmp_tx_';
outdir = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\ExtractedFeatures\';

TESTING = 'E:\Eisemann\Projects\SWG_Feature_Extraction\TESTING\OUTPUT\';
tSLbase = 'TEST_2009_ncmp_tx_';
Segment1 = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\ExtractedFeatures\ShorelineSections\Seg1\';
Segment2 = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\ExtractedFeatures\ShorelineSections\Seg2\';
Segment3 = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\ExtractedFeatures\ShorelineSections\Seg3\';
Segment4 = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\ExtractedFeatures\ShorelineSections\Seg4\';

TRANSECTS = 'E:\Eisemann\Projects\SWG_Feature_Extraction\data\ExtractedFeatures\TransectLines\';
BASE = Segment4;

addpath(genpath(indir))

%% CREATING LOOP FOR FILES

GRIDSstruct = dir([indir 'NCMP_TX_BareEarth_1mGrid_2009\*.tif']);
GRIDS = {GRIDSstruct(:).name}; % making cell array of file names
nn = length(GRIDS);
close all

% for jj = 1:nn % All
% for jj = 1:13 % SEG 1 (South)
% for jj = 14:40 % SEG 2
% for jj = 41:74 % SEG 3
% for jj = 75:118 % SEG 4 (North)

% figure
Transects = [];
for jj = 1:13
   
   tiff = char(GRIDS(jj));
   blocknum = tiff(14:16);
   blocknum
   tiffpath = strcat(indir, 'NCMP_TX_BareEarth_1mGrid_2009\' , tiff);
   
   % import tiff
   SISMfindfile = strcat(indir, 'SISM\', SLbase, blocknum, '*.shp');
   SISMfilestruct = dir(SISMfindfile);
   SISMfile = char(SISMfilestruct(1).name);
   SISMfilepath = strcat(indir, 'SISM\', SISMfile);
   
   % import shoreline
   SMfindfile = strcat(indir, 'SMOOTH\', SLbase, blocknum, '*.shp' );
   SMfilestruct = dir(SMfindfile);
   SMfile = char(SMfilestruct(1).name);
   SMfilepath = strcat(indir, 'SMOOTH\', SMfile);
   
   % Loading TIFF
   [z, r, bbox] = geotiffread(tiffpath);
   info = geotiffinfo(tiffpath);
   
   % Loading simplified-smoothed shoreline (SISM)for shoreline slope
   readSISM = shaperead(SISMfilepath);
   
   % Loading smoothed shoreline (SM) for beachwidth calc 
   readSM = shaperead(SMfilepath);
   
   % concatenating the segments into single vertical X,Y vectors per box
   SS = size(readSISM);
   NumSegs = SS(1);
   slx = [];
   sly = [];
   for ss = 1:NumSegs
       tempX = readSISM(ss).X;
       tempY = readSISM(ss).Y;      
       slx = vertcat(tempX', slx);
       sly = vertcat(tempY', sly);

   end
   % creating linear fit for shoreline
   % feed slx and yfit into ShoreNormalTransects.m
   if abs(max(slx)-min(slx))<abs(max(sly)-min(sly)) % greater change in latitude(Y)
     nans = isnan(slx)|isnan(sly);
     PF = polyfit(sly(~nans),slx(~nans),1);
     xfit = PF(1).*sly+PF(2);
     slx = xfit;
     yfit = sly;
   elseif abs(max(slx)-min(slx))>abs(max(sly)-min(sly)) % greater change in longitude(X)
     nans = isnan(slx)|isnan(sly);
     PF = polyfit(slx(~nans),sly(~nans),1);
     yfit = PF(1).*slx+PF(2);
   end
   
%     nanind = find(isnan(slx));
%     slx(nanind) = [];
%     yfit(nanind) = [];
    [yfit, ind]=sort(yfit,'ascend');
    slx = slx(ind);
%% Clipping Grids 

[nR, nC] = size(z);
noDataValue = min( min( z ) );

% ---- Make location vector
x1 = ( r(3,1) : r(2,1) : r(3,1) + (nC-1)*r(2,1) );
y1 = ( r(3,2) : r(1,2) : r(3,2) + (nR-1)*r(1,2) )';

% ---- Make location matrices
x = repmat(x1,nR,1);
y = repmat(y1,1,nC);

% ---- Find empty space in top, bottom and sides of grids
rInd = [];
dirs = {'fwd','bwd'};
for iDir = 1:2
    
    switch dirs{iDir}
        case 'fwd'
            rind = 1:nR;
            cind = 1:nC;
        case 'bwd'
            rind = nR:-1:1;
            cind = nC:-1:1;
    end
    
    rInd = [];
    for R = rind
        dataInd = find( z(R,:) ~= noDataValue, 1, 'first' );
        if isempty( dataInd )
            rInd = [ rInd, R];
        end
    end
    
    cInd = [];
    for c = cind
        dataInd = find( z(:,c) ~= noDataValue, 1, 'first' );
        if isempty( dataInd )
            cInd = [ cInd, c];
        end
    end
end

% % ---- Clip the grids
% z( rInd, : ) = [];
% x( rInd, : ) = [];
% y( rInd, : ) = [];
% z( :, cInd ) = [];
% x( :, cInd ) = [];
% y( :, cInd ) = [];

% ---- Replace no data values with nan
z( z == noDataValue ) = NaN; 
z = double(z);   

% ---- Feed to ShoreNormalTransects.m
[Xsn,Ysn,Zsn] = ShoreNormalTransects(slx, yfit, pSkip, z, r);


% --- PLOT
% figure
% hold on
% daspect([1 1 1])
% pcolor(x,y,z)
% shading flat
% title(['Block ' blocknum])
% plot(slx, yfit, 'r-')
% --- PLOT

% --- Pre-allocating
ntrans = length(Zsn(:,1));

toe_x = zeros(1,length(ntrans)); 
toe_y = zeros(1,length(ntrans));
toe_z = zeros(1,length(ntrans));

dune_c = zeros(1,length(ntrans));
dune_x = zeros(1,length(ntrans));
dune_y = zeros(1,length(ntrans));

duneHigh_c = zeros(1,length(ntrans));
duneHigh_x = zeros(1,length(ntrans));
duneHigh_y = zeros(1,length(ntrans));
                       
%% Extracting Features

pInd = 1;
for iProfile = 1 : ntrans % --- ntrans = number of rows from function output
%  for iProfile = iProfile % Testing
    iProfile;
    pInd = pInd +1;

X = Xsn(iProfile, :);
Y = Ysn(iProfile, :);
Z = Zsn(iProfile, :);

% ---- Plots transects on map
% hold on
% plot(X,Y,'k-')
% ---- Save transects to structure array
Trans = struct('ID',iProfile,'Geometry','Line','X',X,'Y',Y);
Transects = [Transects,Trans];

% % ---- replacing Nan with empty set
nanind = find(isnan(Z));
Z(nanind) = [];
X(nanind) = [];
Y(nanind) = [];
if length(Z)<20 % --- Remove if little data
    Z = [];
    X = [];
    Y = [];
continue % --- exits loops if empty
end

% ---- Creating meters cross-shore vector
format long
lat1 = Y(1); % -- Most shoreward point
lon1 = X(1);
lat2 = Y(end); % -- Most inland point
lon2 = X(end);
[arclen, ~ ] = distance (lat1, lon1, lat2, lon2); % ---- Find arc distance, degrees
pdist =  distdim (arclen, 'degrees', 'meters'); % ---- Convert to meters
CSM = (0 : pdist/(length(X)-1) : pdist); % --- Cross Shore Meters
% ---- zero is seaward edge of data

% ---- Smooth each profile z values 
S = 10;
z2 = smoothn(Z, S);

    % ---- DUNE PEAKS ---- %
    [pks, locs] = findpeaks(z2, 'MinPeakProminence', P);
    dune_loc_i = find(pks > 1.2);      % ---- find location where peaks > 1
    dune_loc = locs(dune_loc_i);
    test = isempty(dune_loc);
    if test == 1
        dune_c(pInd) = nan;
        dune_x(pInd) = nan;
        dune_y(pInd) = nan;
    else
        if max(pks(dune_loc_i)) < 1.2
            dune_c(pInd) = nan;
            dune_x(pInd) = nan;
            dune_y(pInd) = nan;
        else
            % ---- find dune crest ---- %
            dune_c(pInd) = pks(dune_loc_i(1));
            dune_x(pInd) = X(dune_loc(1));
            dune_y(pInd) = Y(dune_loc(1));

            % --- find highest crest ---%
            duneHigh_c(pInd) = max(pks);
            maxloc = find(pks == max(pks));
            duneHigh_x(pInd) = X(locs(maxloc));
            duneHigh_y(pInd) = Y(locs(maxloc));

            % ---- DUNE TOE  ---- %
            crestDist = CSM(dune_loc(1));% --- distance of crest from shore in meters
            ToeWindow = 40; 
            if crestDist > ToeWindow % -- if the distance from the crest to shore less than window
                ToeWindow = 40; 
            else
                ToeWindow = crestDist; 
            end 
            toeSearchInd = find(abs(CSM - (crestDist - ToeWindow))<1, 1, 'first'); % --- find distance from shore to stop looking for toe
            shore2crestX  = X(toeSearchInd:dune_loc(1)); % --- location of dune to end
            shore2crestY  = Y(toeSearchInd:dune_loc(1));  % --- location of dune to end
            shore2crestZ  = z2(toeSearchInd:dune_loc(1)); % --- uses smoothn z values

            % ---- Second Derivative of crest2shoreZ 
            d2CSZ = del2(shore2crestZ);

            [d2pks, d2pklocs] = findpeaks(d2CSZ); % --- finding peaks in d2
            [d2trs, d2trlocs] = findpeaks(-d2CSZ); % --- finding troughs

            for ii = (0.9 : - 0.0001 : 0.001)
                toe_loc_i = find(d2pks > ii, 1, 'last'); % --- direction change
                % --- next iteration if toe_loc_i empty --- %
                test = isempty (toe_loc_i);
                if test == 0
                    break % --- to keep current toe_loc_i
                else
                    continue
                end
            end

            % --- Finding nearest seaward trough --- %
            test = isempty (toe_loc_i);
            if test == 0
                % --- Incase there are some that have no large enough peaks in D2
                SWTRind = find(d2trlocs < d2pklocs(toe_loc_i), 1, 'last');
                SWTR = d2trlocs(SWTRind); %seaward trough loc
                % --- next iteration if trough is certain size --- %
                % --- relative to D2 peak found for toe --- %
                if toe_loc_i > 1 % --- if not already most seaward
                    if (abs(d2CSZ(SWTR))/d2CSZ(d2pklocs(toe_loc_i))>0.5)
                        toe_loc_i = (toe_loc_i - 1); 
                    end
                end
                
                toe_loc = d2pklocs(toe_loc_i);
                test = isempty(toe_loc);
                
                if test == 1 % --- if toe_loc empty
                    toe_x(pInd) = NaN;
                    toe_y(pInd) = NaN;
                    toe_z(pInd) = NaN;
                else
                    toe_x(pInd) = shore2crestX(toe_loc);
                    toe_y(pInd) = shore2crestY(toe_loc);
                    toe_z(pInd) = shore2crestZ(toe_loc); 

                    % --- Beach Width ---%
                    [arclen, ~] = distance(shore2crestY(toe_loc), shore2crestX(toe_loc), Y(1), X(1));
                    BeachWid = distdim(arclen, 'deg', 'm'); % --- Length of ave shoreline (meters)
                    BeachWidth(pInd) = BeachWid;
                end
                    else % --- if toe_loc_i empty
                        toe_x(pInd) = NaN;
                        toe_y(pInd) = NaN;
                        toe_z(pInd) = NaN;
            end
        end 
    end
end
%             end % --- ending if z2 not empty keep going

dune_x(dune_x == 0) = NaN;
dune_y(dune_y == 0) = NaN; % For some reason there are zeros
toe_x(toe_x == 0) =  NaN;
toe_y(toe_y == 0) = NaN;
duneHigh_x(duneHigh_x == 0) = NaN;
duneHigh_y(duneHigh_y == 0) = NaN;

% plotting on map
plot(dune_x,dune_y,'bo')

%% Writing Files

if savefile == 'y'
%     
%     mpCrest = mappoint(dune_x',dune_y','Dune_Crest', dune_c');
%     shapewrite(mpCrest, [outdir 'DuneCrest\' SLbase blocknum '_DuneCrest'])
% 
%     mpToe = mappoint(toe_x', toe_y', 'Dune_Toe', toe_z');
%     shapewrite(mpToe, [outdir 'DuneToe\' SLbase blocknum '_DuneToe'])
% 
%     mpHigh = mappoint( duneHigh_x', duneHigh_y', 'Dune_High', duneHigh_c');
%     shapewrite(mpHigh, [outdir 'DuneHigh\' SLbase blocknum '_DuneHigh'5])
% 
%     % --- Making .prj file --- %
%     GEOGCS = 'GEOGCS["GCS_NAD_1983_2011",DATUM["D_NAD_1983_2011",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433],AUTHORITY["EPSG",6318]]';
%     save([ outdir 'DuneCrest\' SLbase blocknum '_DuneCrest.prj'],'GEOGCS')
%     save([ outdir 'DuneToe\' SLbase blocknum '_DuneToe.prj'],'GEOGCS')
%     save([ outdir 'DuneHigh\' SLbase blocknum '_DuneHigh.prj'],'GEOGCS')
    
% TESTING or SEGMENTS
    mpCrest = mappoint(dune_x',dune_y','Dune_Crest', dune_c');
    shapewrite(mpCrest, [BASE  SLbase blocknum '_DuneCrest'])

    mpToe = mappoint(toe_x', toe_y', 'Dune_Toe', toe_z');
    shapewrite(mpToe, [BASE  SLbase blocknum '_DuneToe'])

    mpHigh = mappoint( duneHigh_x', duneHigh_y', 'Dune_High', duneHigh_c');
    shapewrite(mpHigh, [BASE  SLbase blocknum '_DuneHigh'])
   
    shapewrite(Transects,[ TRANSECTS 'Seg4Transects.shp'])
    
%     C = [dune_x;dune_y;dune_c];
%     CrestID = fopen([BASE SLbase blocknum 'Dune_Crest.txt'],'w');
%     fprintf(CrestID, '%12.8f %12.8f %12.8f\n',C);
%     
%     T = [toe_x;toe_y;toe_z];
%     ToeID = fopen([BASE SLbase blocknum 'Dune_Toe.txt'],'w');
%     fprintf(ToeID, '%12.8f %12.8f %12.8f\n', T);
%     
%     H = [duneHigh_x;duneHigh_y;duneHigh_c];
%     HighID = fopen([BASE SLbase blocknum 'Dune_High.txt'],'w');
%     fprintf(HighID, '%12.8f %12.8f %12.8f\n', H);

    % --- Making .prj file --- %
%     GEOGCS = 'GEOGCS["GCS_NAD_1983_2011",DATUM["D_NAD_1983_2011",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433],AUTHORITY["EPSG",6318]]';
%     save([ BASE  SLbase blocknum '_DuneCrest.prj'],'GEOGCS')
%     save([ BASE  SLbase blocknum '_DuneToe.prj'],'GEOGCS')
%     save([ BASE  SLbase blocknum '_DuneHigh.prj'],'GEOGCS')

end % --- end if yes to save file
end


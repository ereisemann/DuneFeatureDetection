% ---- From Lauren Dunkin 02/2017
% ---- edited by Eve Eisemann 03/2017
% ---- edited by Michael Hartman 04/2017
% ---- edited by Eve Eisemann 04/2017, 05/2017
clear,clc

% ---- Save a file? y/n ----%
savefile = 'y';

% --- Change this to skip transects - meters
% pSkip = 10;
% --- Minimum crest height ---%
minelev = 2.0;

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
indir = 'F:\DUNEX\AirborneLidar\';
outdir = 'F:\DUNEX\Features_Output\';
SLbase = 'NC_postSandy_USGS_2012_2_';
addpath(genpath(indir))
addpath(genpath(outdir))

%--TRANSECTS--%
% % TRANSpath =([indir 'TransectSegs\']);
% TRANSpath = 'F:\DUNEX\VolChangeAssessments';
% TRANSstruct = dir([TRANSpath '*.shp']);
% TRANS = {TRANSstruct(:).name};

transpath = ('F:\DUNEX\VolChangeAssessments\NOBX_transects_geographic.shp');
trans = shaperead(transpath,'UseGeoCoords',true);
Ntrans = length(trans);
addpath(genpath('F:\DUNEX\VolChangeAssessments\'))

%--GRIDS--%
GRIDSpath = ([indir '2012_PostSandy_USGS_NC_BareEarth_1mGrid\postMatthew_2012_1kmBlocks\']);
GRIDSstruct = dir([GRIDSpath '*.TIF']);
GRIDS = {GRIDSstruct(:).name}'; % making cell array of file names


%% CREATING LOOP FOR FILES

nn = length(GRIDS);
close all

% figure
Transects = [];


 for jj = 1:nn
% for jj = 2
   
   tiff = char(GRIDS(jj));
   blocknum = tiff(23:24); % CHANGE for different file names
%    blocknum = jj;
   blocknum
   
   tiffpath = strcat(GRIDSpath, tiff);
   
%    % import tiff
%    SISMfindfile = strcat(indir, 'SISM\', SLbase, blocknum, '*.shp');
%    SISMfilestruct = dir(SISMfindfile);
%    SISMfile = char(SISMfilestruct(1).name);
%    SISMfilepath = strcat(indir, 'SISM\', SISMfile);
   
%    % import shoreline
%    SMfindfile = strcat(indir, 'SMOOTH\', SLbase, blocknum, '*.shp' );
%    SMfilestruct = dir(SMfindfile);
%    SMfile = char(SMfilestruct(1).name);
%    SMfilepath = strcat(indir, 'SMOOTH\', SMfile);
   
   %--Loading TIFF--%
   %%%%%%%%%%%%%%%%%%
   [z, r, bbox] = geotiffread(tiffpath);
   info = geotiffinfo(tiffpath);

   %--Loading TRANSECTS--%
   %--MAKE SURE FILE NAME IS CORRECT!--%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    transpath = strcat(TRANSpath,char(TRANS(jj)));

%    transpath = strcat(TRANSpath,char(TRANS(blocknum)));
%    trans = shaperead(transpath,'UseGeoCoords',true);
   
   
%    % Loading simplified-smoothed shoreline (SISM)for shoreline slope
%    readSISM = shaperead(SISMfilepath);
%    
%    % Loading smoothed shoreline (SM) for beachwidth calc 
%    readSM = shaperead(SMfilepath);
   
   % concatenating the segments into single vertical X,Y vectors per box
%    SS = size(readSISM);
%    NumSegs = SS(1);
%    slx = [];
%    sly = [];
%    for ss = 1:NumSegs
%        tempX = readSISM(ss).X;
%        tempY = readSISM(ss).Y;      
%        slx = vertcat(tempX', slx);
%        sly = vertcat(tempY', sly);
% 
%    end
%    % creating linear fit for shoreline
%    % feed slx and yfit into ShoreNormalTransects.m
%    if abs(max(slx)-min(slx))<abs(max(sly)-min(sly)) % greater change in latitude(Y)
%      nans = isnan(slx)|isnan(sly);
%      PF = polyfit(sly(~nans),slx(~nans),1);
%      xfit = PF(1).*sly+PF(2);
%      slx = xfit;
%      yfit = sly;
%    elseif abs(max(slx)-min(slx))>abs(max(sly)-min(sly)) % greater change in longitude(X)
%      nans = isnan(slx)|isnan(sly);
%      PF = polyfit(slx(~nans),sly(~nans),1);
%      yfit = PF(1).*slx+PF(2);
%    end
%    
% %     nanind = find(isnan(slx));
% %     slx(nanind) = [];
% %     yfit(nanind) = [];
%     [yfit, ind]=sort(yfit,'ascend');
%     slx = slx(ind);
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

%    ptAX = [];
%    ptAY = [];
%    ptBX = [];
%    ptBY = [];
   
    Zsn = [];
    Xsn = [];
    Ysn = [];
   
% -- looping through transects in segment ---%    
pInd = 1;    
toex = []; 
toey = [];
toez = [];

dunec = [];
dunex = [];
duney = [];

duneHighc = [];
duneHighx = [];
duneHighy = [];

% --- Finding transects within segment
in_1 = [];
in = [];
for t  = 1:Ntrans
     transLon = trans(t).Lon(:);
     transLat = trans(t).Lat(:);
     transLon(isnan(transLon)) = [];
     transLat(isnan(transLat)) = [];
    in_1 = inpolygon(transLon(2), transLat(2), bbox(:,1),bbox(:,2));
    in = vertcat(in, in_1);
end 

trans_ind = find(in  == 1);

% --- Cycling through only those transects 

for tt = min(trans_ind):max(trans_ind)
% for tt = 1100
% for tt = 93
    
%        [maxLon, ind] = max(trans(tt).Lon(:));
%        ptAX = horzcat(ptAX, maxLon);
%        ptAY = horzcat(ptAY, trans(tt).Lat(ind));
%        [minLon, ind] = min(trans(tt).Lon(:));
%        ptBX = horzcat(ptBX, minLon);
%        ptBY = horzcat(ptBY, trans(tt).Lat(ind));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%below: removed and replaced 5/30/2019%%%%%%
%       [ptAX, ind] = max(trans(tt).Lon(:));
%        ptAY = trans(tt).Lat(ind);
%       [ptBX, ind] = min(trans(tt).Lon(:));
%        ptBY = trans(tt).Lat(ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     transLon = trans(tt).Lon(:);
     transLat = trans(tt).Lat(:);
     transLon(isnan(transLon)) = [];
     transLat(isnan(transLat)) = [];

     ptAX = transLon(1);
     ptAY = transLat(1);
     ptBX = transLon(end);
     ptBY = transLat(end);
        % ---- Looking at shore-perpendicular transects
        disp([num2str(tt),'/',num2str(max(trans_ind))])
        [ZZ, ~, YY, XX] = mapprofile(z, r, [ptAY, ptBY], [ptAX, ptBX]);
        test = isempty(ZZ);
        if test == 1
            continue
        else
       % --- convert column vector to row vector
            Yi = YY';
            Xi = XX';
            Zi = ZZ';
            
%        % --- append row vector to bottom of storage matrix
%             Ysn = vertcat(Ysn, Yi);
%             Xsn = vertcat(Xsn, Xi);
%             Zsn = vertcat(Zsn, Zi);
        end
%    end
   
    
%     % --- create vector of nan-separated segments 
%     Yseg_mat = [ptAY(1:Ntrans); ptBY(1:Ntrans); nan(1,Ntrans)]; % could loop, but using indexing trick 
%     Xseg_mat = [ptAX(1:Ntrans); ptBX(1:Ntrans); nan(1,Ntrans)]; % could loop, but using indexing trick
%     
%     Yseg_vec = Yseg_mat(:);
%     Xseg_vec = Xseg_mat(:);
%     
%     % --- make a single call to mapprofile (went from 575 seconds to 3.6 seconds)
%     [ZZ, ~, YY, XX] = mapprofile(z, r, Yseg_vec, Xseg_vec);
%     
%     % Note: the input to mapprofile is a vector of nan-separated segments,
%     % then XX, YY, and ZZ will have NaNs separating the output segments. 
%     
%     % If the interpolated transects are different lengths, then we might 
%     % exclude them, although this case hasn't been encountered yet.
%     nan_idx = find(isnan(YY));
%     idxDiff = diff(nan_idx);
%     uniqueDiff = unique(idxDiff);
%         
%     if length(uniqueDiff) == 1
%           
%         YY(nan_idx) = [];
%         XX(nan_idx) = [];
%         ZZ(nan_idx) = [];
%         
%         Ysn = reshape(YY,[],Ntrans);
%         Xsn = reshape(XX,[],Ntrans);
%         Zsn = reshape(ZZ,[],Ntrans);
%         
%         Ysn = Ysn';
%         Xsn = Xsn';
%         Zsn = Zsn';
%         
%    else
%        return
%    end
%     
%     if flag_testing == 1
%         % load method_1 results
%         method1 = load('Method1.mat');
%     end
% end % ---END transect extraction 

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
% ntrans = length(Zsn(:,1));
% ntrans = Ntrans;
% 
% toe_x = []; 
% toe_y = [];
% toe_z = [];
% 
% dune_c = [];
% dune_x = [];
% dune_y = [];
% 
% duneHigh_c = [];
% duneHigh_x = [];
% duneHigh_y = [];
                       
%% Extracting Features for each transect

% pInd = 1;
% for iProfile = 1 : ntrans % --- ntrans = number of rows from function output
%  for iProfile = iProfile % Testing
    iProfile = tt;
%     pInd = pInd +1;
    pInd = tt;

% X = Xsn(iProfile, :);
% Y = Ysn(iProfile, :);
% Z = Zsn(iProfile, :);

pos = find(Zi>0); % finding elevations below zero 

X = Xi(pos); % excluding bathy data from X Y Z
Y = Yi(pos);
Z = Zi(pos);
% ---- Plots transects on map
% hold on
% plot(X,Y,'k-')
% ---- Save transects to structure array
% Trans = struct('ID',iProfile,'Geometry','Line','X',X,'Y',Y);
% Transects = [Transects,Trans];

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
    dune_loc_i = find(pks > minelev);      % ---- find location where peaks > minelev
    dune_loc = locs(dune_loc_i);
    test = isempty(dune_loc);
    if test == 1
%         dune_c(pInd) = nan;
%         dune_x(pInd) = nan;
%         dune_y(pInd) = nan;
        dune_c = nan;
        dune_x = nan;
        dune_y = nan;
    else
        if max(pks(dune_loc_i)) < 1.2
            dune_c = nan;
            dune_x = nan;
            dune_y = nan;
        else
            % ---- find dune crest ---- %
            dune_c = pks(dune_loc_i(1));
            dune_x = X(dune_loc(1));
            dune_y = Y(dune_loc(1));
            
            dunec = horzcat(dunec,dune_c);
            dunex = horzcat(dunex,dune_x);
            duney = horzcat(duney, dune_y);

            % --- find highest crest ---%
            duneHigh_c = max(pks);
            maxloc = find(pks == max(pks));
            duneHigh_x = X(locs(maxloc));
            duneHigh_y = Y(locs(maxloc));
            
            duneHighc = horzcat(duneHighc, duneHigh_c);
            duneHighx = horzcat(duneHighx, duneHigh_x);
            duneHighy = horzcat(duneHighy, duneHigh_y);

            % ---- DUNE TOE  ---- %
            crestDist = CSM(dune_loc(1));% --- distance of crest from shore in meters
            ToeWindow = 50; 
            if crestDist > ToeWindow % -- if the distance from the crest to shore less than window
                ToeWindow = 50; 
            else
                ToeWindow = crestDist; 
            end 
            toeSearchInd = find(abs(CSM - (crestDist - ToeWindow))<10, 1, 'first'); % --- find distance from shore to stop looking for toe
            shore2crestX  = X(toeSearchInd:dune_loc(1)); % --- location of dune to end
            shore2crestY  = Y(toeSearchInd:dune_loc(1));  % --- location of dune to end
            shore2crestZ  = z2(toeSearchInd:dune_loc(1)); % --- uses smoothn z values
            
            testz = isempty(shore2crestZ);
            if testz == 1
               continue  
            end

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
            if test == 0 % I think this should actually be 1 - TEST ME
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
                    toe_x = NaN;
                    toe_y = NaN;
                    toe_z = NaN;
                else
                    toe_x = shore2crestX(toe_loc);
                    toe_y = shore2crestY(toe_loc);
                    toe_z = shore2crestZ(toe_loc); 
                    
                    toex = horzcat(toex,toe_x); 
                    toey = horzcat(toey,toe_y);
                    toez = horzcat(toez,toe_z);

%                     % --- Beach Width ---%
%                     [arclen, ~] = distance(shore2crestY(toe_loc), shore2crestX(toe_loc), Y(1), X(1));
%                     BeachWid = distdim(arclen, 'deg', 'm'); % --- Length of ave shoreline (meters)
%                     BeachWidth(pInd) = BeachWid;
                end
                    else % --- if toe_loc_i empty
                        toe_x = NaN;
                        toe_y = NaN;
                        toe_z = NaN;
            end
        end 
    end

end
%             end % --- ending if z2 not empty keep going

dunex(dunex == 0) = NaN;
duney(duney == 0) = NaN; % For some reason there are zeros
toex(toex == 0) =  NaN;
toey(toey == 0) = NaN;
duneHighx(duneHighx == 0) = NaN;
duneHighy(duneHighy == 0) = NaN;

% plotting on map
% plot(dunex,duney,'bo')

%% Writing Files
% if savefile == 'y'

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

    if isempty(dunec) % if no dune crests were found in this segment, move on!
        ('EMPTY')
        continue
    end

    mpCrest = mappoint(dunex',duney','Dune_Crest', dunec');
    shapewrite(mpCrest, [outdir 'Crest\' SLbase blocknum '_DuneCrest'])

    mpToe = mappoint(toex', toey', 'Dune_Toe', toez');
    shapewrite(mpToe, [outdir 'Toe\' SLbase blocknum '_DuneToe'])

    mpHigh = mappoint( duneHighx', duneHighy', 'Dune_High', duneHighc');
    shapewrite(mpHigh, [outdir 'High\' SLbase blocknum '_DuneHigh'])
   
    
    
    
%     shapewrite(Transects,[ TRANSECTS 'Seg4Transects.shp'])
    
%     C = [dune_x;dune_y;dune_c];
%     CrestID = fopen([BASE SLccbase blocknum 'Dune_Crest.txt'],'w');
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

% end % --- end if yes to save file
end


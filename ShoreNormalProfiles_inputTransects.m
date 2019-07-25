function [Xsn, Ysn, Zsn] = ShoreNormalProfiles_inputTransects(slX, slY, pSkip, z, r)

% --- INPUT: 'SLX' and 'SLY' are the X and Y coords for linear shoreline
% --- pSkip is the interval in meters at which transects should be made
% --- 'z' is the elevation data grid. 
% --- 'r' is the referencing vector from the original import of the z data
% --- OUTPUT: X, Y, and Z grids oriented such that transects can be
% --- produced perpendicular to the shoreline input.
% --- Each shore-perpendicular transect is represented as a row 
% --- The number of rows therefore indicates the number of transects.
% --- The number of points along the transect is determined by the
% --- referencing vector 'r' which can be changed in this code. 
% 
% --- Eve Eisemann 2017

% --- Removing NaNs
Xnanind = isnan(slX);
Ynanind = isnan(slY);
slX(Xnanind) = [];
slY(Ynanind) = [];

% --- Finding max and min values
maxX = max(slX);
minX = min(slX);
maxY = max(slY);
minY = min(slY);

% --- Calculating change in x and change in y
[arclenX, ~] = distance(maxY, maxX, maxY, minX); % -- Y/lat constant
[arclenY, ~] = distance(maxY, maxX, minY, maxX); % -- X/lon constant

delXkm = distdim(arclenX, 'degrees', 'kilometers'); 
delYkm = distdim(arclenY, 'degrees', 'kilometers');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- ORIENTING SHORELINE --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- ordered South to North or 
% --- West to East

%%%%%%%%%%%%%%%%%%%%%%%
% --- East - West --- %
%%%%%%%%%%%%%%%%%%%%%%%

if delXkm >= delYkm % --- diff in x or diff in y larger?
    minind = find(slX == min(slX));
    % --- If min value does not exist in the Western 1/4 of shoreline, flip
    if minind >= floor(length(slX)/4)
        SLX = fliplr(slX);
        SLY = fliplr(slY);
        
    else
        if minind < floor(length(slX)/4)
            SLX = slX;
            SLY = slY;
        end
    end
else
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % --- North-South --- %
    %%%%%%%%%%%%%%%%%%%%%%%
    if delYkm > delXkm % --- In case of ~N-S shoreline
        minind = find(slY == min(slY));
        % --- If min value does not exist in the Southern 1/4 of shoreline, flip
        if minind >= floor(length(slY)/4)
            SLX = fliplr(slX);
            SLY = fliplr(slY);
        else
            if minind < floor(length(slY)/4)
                SLX = slX;
                SLY = slY;
            end
        end
    end
end % --- End shoreline orientation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Calculating slopes --- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = (SLY(2:end) - SLY(1:end-1)) ./ (SLX(2:end) - SLX(1:end-1)); % vectorized
% could be changed now that feeding linear fit instead. 

Mmovmean = movmean(M, 20); % --- buffer extreme values
meanslope = nanmean(Mmovmean);

% --- mean slope and ~central point to create ave shoreline for segment
[arclen, ~] = distance(SLY(1), SLX(1), SLY(end), SLX(end));
lenSL = distdim(arclen, 'deg', 'm'); % --- Length of ave shoreline (meters)
Ntrans = floor(lenSL/pSkip); % --- number of transects to produce. 

% --- Slope for transects
Mtrans = (-1/meanslope);

% --- X Y coords for seaward transect points
ptX = (SLX(1) : (SLX(end)-SLX(1))/Ntrans : SLX(end));
ptY = (SLY(1) : (SLY(end)-SLY(1))/Ntrans : SLY(end));

% --- 1 kilometer inland 
% LenIn = distdim(0.5, 'km', 'deg');%maybe change to 0.5 km instead
LenIn = distdim(1, 'km', 'deg');

% --- X Y coords for landward transect points
% --- Needs adjustment to be universal, here distance is negative, 
% --- So only works when shoreline is to the right.  
ptBX = ptX - (LenIn/sqrt((LenIn+Mtrans.^2)));
ptBY = ptY - ((LenIn.*Mtrans)/sqrt((LenIn+Mtrans.^2)));

% test = isreal(ptBX);
% if test == 0
%     ptBX = ptX + (-LenIn/sqrt(-LenIn+Mtrans.^2));
%     ptBY = ptY + ((-LenIn.*Mtrans)/sqrt(-LenIn+Mtrans.^2));
% end

% --- 0.5 km seaward
LenOut = distdim(1, 'km', 'deg');

% --- X Y coords for seaward point of transect, seaward of ave shoreline
ptAX = ptX + (LenOut/sqrt(LenOut+Mtrans.^2));
ptAY = ptY + ((LenOut.*Mtrans)/sqrt(LenOut+Mtrans.^2));

%%% SPECIAL EXCEPTION FOR BLOCKS 116-118
% A to B goes inland to shore for these blocks, so switch
% ptAXr = ptAX;
% ptAYr = ptAY;
% ptBXr = ptBX;
% ptBYr = ptBY;
% 
% ptAX = ptBXr;
% ptAY = ptBYr;
% ptBX = ptAXr;
% ptBY = ptAYr;
%%% END SPECIAL EXCEPTION, TURN OFF FOR ALL OTHER BLOCKS

flag_testing = 0;
profile_method = 2; % 1 - Original; 2 - Testing 
if profile_method == 1
    
    Zsn = [];
    Xsn = [];
    Ysn = [];
    
    for ii = 1 : Ntrans
        % ---- Looking at shore-perpendicular transects
        disp([num2str(ii),'/',num2str(Ntrans)])
        [ZZ, ~, YY, XX] = mapprofile(z, r, [ptAY(ii), ptBY(ii)], [ptAX(ii), ptBX(ii)]);
        test = isempty(ZZ);
        if test == 1
            continue
        else
            % --- convert column vector to row vector
            Yi = YY';
            Xi = XX';
            Zi = ZZ';
            
            % --- append row vector to bottom of storage matrix
            Ysn = vertcat(Ysn, Yi);
            Xsn = vertcat(Xsn, Xi);
            Zsn = vertcat(Zsn, Zi);
        end
    end
    
elseif profile_method == 2
    
    % --- pre-allocate, but might be unnecessary
    Zsn = [];
    Xsn = [];
    Ysn = [];
    
    % --- create vector of nan-separated segments 
    Yseg_mat = [ptAY(1:Ntrans); ptBY(1:Ntrans); nan(1,Ntrans)]; % could loop, but using indexing trick 
    Xseg_mat = [ptAX(1:Ntrans); ptBX(1:Ntrans); nan(1,Ntrans)]; % could loop, but using indexing trick
    
    Yseg_vec = Yseg_mat(:);
    Xseg_vec = Xseg_mat(:);
    
    % --- make a single call to mapprofile (went from 575 seconds to 3.6 seconds)
    [ZZ, ~, YY, XX] = mapprofile(z, r, Yseg_vec, Xseg_vec);
    
    % Note: the input to mapprofile is a vector of nan-separated segments,
    % then XX, YY, and ZZ will have NaNs separating the output segments. 
    
    % If the interpolated transects are different lengths, then we might 
    % exclude them, although this case hasn't been encountered yet.
    nan_idx = find(isnan(YY));
    idxDiff = diff(nan_idx);
    uniqueDiff = unique(idxDiff);
        
    if length(uniqueDiff) == 1
          
        YY(nan_idx) = [];
        XX(nan_idx) = [];
        ZZ(nan_idx) = [];
        
        Ysn = reshape(YY,[],Ntrans);
        Xsn = reshape(XX,[],Ntrans);
        Zsn = reshape(ZZ,[],Ntrans);
        
        Ysn = Ysn';
        Xsn = Xsn';
        Zsn = Zsn';
        
    else
        return
    end
    
    if flag_testing == 1
        % load method_1 results
        method1 = load('Method1.mat');
    end
end 
end

clearvars; close all; clc;

% select project directory
% projDir = fullfile(filesep, 'Users', 'fabian', 'Work', 'Projects', 'PEANUTS_PRIZE');
projDir = fullfile(filesep, 'media', 'aidan', 'STORAGE', 'Work', 'microARC', 'from_Fabian', 'test_case_FramStrait');

% select forcing model
frcModel = 'sinmod';

% select regional setup name
regName = 'FramStrait_dep0-100';

% select experiment name and years
expName = 'AWI-Hausgarten-YEAR';
years = 2017:2018;

% select observation data to be used for calculation of target region
obsDir = fullfile(projDir, 'DATA', 'AWI_Hausgarten');
obsFile = 'AWI-Hausgarten.xlsx';

% distance (in km) by which region around stations should be expanded
radius = 25;

% select trajectory file (with place holders) to be used
traDir = fullfile(projDir, 'DATA', 'particle_trajectories', frcModel, regName);
traFile = 'particles_MODEL_REGION_EXPERIMENT_YEAR_t*_iSub03.mat';

% define topography file to be used
% NOTE: last letter before file ending in 'topoFile' determines resolution
%       - (f)ull, (h)igh, (i)ntermediate, (l)ow, (c)rude
%       - lower resolution significantly reduces run time
topoDir = fullfile(projDir, 'MATLAB', 'map_tools', 'topo_mat');
topoFile = 'AWI-Hausgarten_topo_i.mat';

% set map's center longitude, and lon-lat range
ctrLon = 0;
mapLon = [-50, 50];
mapLat = [60, 90];

% set output directory
outDir = obsDir;

% set figure visibility ('on' recommended for code development/debugging)
vis = 'on';

%% START PROGRAM

%%
% load observation file
try
    [NUM,TXT] = xlsread(fullfile(obsDir,obsFile));
    % There doesn't appear to be any info difference between columns of
    % sampling event number and sampling station... something to ask Marcus
    % ... NOT TRUE, there are differences in sampling times, not positions
    
    obs.time = NUM(:,1);
    % correct time information
    obs.time = obs.time - 1; % one day offset
    obs.time = datevec(obs.time);
    obs.time(:,1) = obs.time(:,1) + 1900; % years not correct    
    obs.yearday = yearday(datenum(obs.time(:,1:3))); % store yearday as well as date
    obs.lat = NUM(:,2);
    obs.lon = NUM(:,3);
    obs.station = TXT(2:end,strcmp(TXT(1,:), 'HG Station name'));
    obs.event = TXT(2:end,strcmp(TXT(1,:), 'Pangaea event label'));
    clear NUM TXT
catch
    error('Error loading file:\n ==> %s', fullfile(obsDir, obsFile));
end

% separate observation data by year
nYears = length(years);
for iy = 1:nYears
    it = obs.time(:,1)==years(iy);
    yStr = sprintf('y%04i', years(iy));
    obs.(yStr).time = obs.time(it,1:3);
    obs.(yStr).yearday = obs.yearday(it);
    obs.(yStr).lat = obs.lat(it);
    obs.(yStr).lon = obs.lon(it);    
    obs.(yStr).station = obs.station(it);
    obs.(yStr).event = obs.event(it);
    % within years, refer to each station and to each event by index numbers
    obs.(yStr).stationNum = zeros(size(obs.(yStr).station));
    obs.(yStr).eventNum = zeros(size(obs.(yStr).event));
    stations = unique(obs.(yStr).station, 'stable');
    events = unique(obs.(yStr).event, 'stable');
    for jj = 1:length(stations)
        obs.(yStr).stationNum(strcmp(obs.(yStr).station, stations{jj})) = jj;
    end
    for jj = 1:length(events)
        obs.(yStr).eventNum(strcmp(obs.(yStr).event, events{jj})) = jj;
    end
    p = polygon_from_stations(obs.(yStr), radius);
    obs.(yStr).reg = [p.lon; p.lat];
    clear it yStr stations events jj p
end


% update trajectory file name
trajFile = replace(fullfile(traDir, traFile), {'MODEL', 'REGION', 'EXPERIMENT'}, ...
                                              {frcModel, regName, expName});
                                          
                                          
% load and store trajectories and associated BIOMAS output metrics
for iy = 1:nYears
    yStr = sprintf('y%04i', years(iy));
    fList = dir(strrep(trajFile, 'YEAR', num2str(years(iy))));
    frcFile = fullfile(fList.folder, fList.name);
    try
        load(frcFile);
    catch
        error('Error loading file:\n ==> %s', frcFile);
    end
    traj.(yStr).nt = size(PDat.x,1);
    traj.(yStr).nx = size(PDat.x,2);
    traj.(yStr).t = PDat.t;
    traj.(yStr).yearday = yearday(traj.(yStr).t);    
    traj.(yStr).x = PDat.x;
    traj.(yStr).y = PDat.y;
    traj.(yStr).profiles = PDat.profiles;
    clear PDat
end    



% For each sampling station: find position of all trajectories at the time
% (yearday) of sampling, then sort trajectories by distance from sampling
% station - the closest are the most relevant. This could be repeated for
% times within some 'reasonable' interval of the sampling time.
% Store the info within the traj structure, creating separate substructures
% for each station - this storage could be improved later...

for iy = 1:nYears
    yStr = sprintf('y%04i', years(iy));    
    for ie = 1:length(unique(obs.(yStr).eventNum)) % loop through events
        eStr = sprintf('event%i', ie);
        event = obs.(yStr).eventNum == ie;
        % I don't understand why the events are not unique (obviously the
        % labelling system on this cruise differs from my own experience of
        % event labelling...). Perhaps it's simply explained by
        % representing multiple different sample types from the same gear
        % deployment... because fortunately the multiple table rows sharing
        % event numbers have the same space-time coords, which is all that
        % really matters.
        sTime = unique(obs.(yStr).yearday(event)); % sample time
        lat = unique(obs.(yStr).lat(event)); % sample position
        lon = unique(obs.(yStr).lon(event));
        % Convert ot radians to calculate distances
        samplePos_rads = deg2rad([lon lat]);

        yd = traj.(yStr).yearday == sTime;
        
        trajLon = traj.(yStr).x(yd,:); % particle positions along trajectories
        trajLat = traj.(yStr).y(yd,:); % at time of sample event
        
        traj.(yStr).(eStr).yearday = sTime;
        traj.(yStr).(eStr).lat = trajLat;
        traj.(yStr).(eStr).lon = trajLon;
        
        trajPos_rads = deg2rad([trajLon; trajLat]');
        
        % distance (km) between sample site and particle at time of event
        Dist = distance_on_earth(samplePos_rads, trajPos_rads)' / 1000;
        
        [~,I] = sort(Dist);
        traj.(yStr).(eStr).distance = Dist; % distance between particle and station at time of sampling event        
        % rank trajectories by closeness to sampling events
        [~,rankDist] = sort(I);
        traj.(yStr).(eStr).rank = rankDist;
        
    end
    clear yStr eStr event sTime lat lon samplePos_rads yd trajLon trajLat trajPos_rads Dist I rankDist
end


% Find ratios  of large to small plankton
for iy = 1:nYears
    yStr = sprintf('y%i', years(iy));
    traj.(yStr).profiles.P_ratio = traj.(yStr).profiles.PL ./ traj.(yStr).profiles.PS;
end

% plotVar = 'PS';
depthLayer = 1;
labelVar = {'temp', 'salt', 'NO3', 'P_ratio';
    'temperature', 'salinity', 'NO_3', 'plankton L/S';
    [char(176) 'C'], 'psu', 'mmol N / m^3', 'large to small plankton ratio'};

nVars = size(labelVar,2); % number of metrics to plot


% Choose a min & max number of particles to use per sampling event. A
% minimum to preserve potential variability in forcing data, and a maximum
% to limit computational effort when optimising parameters. (The minimum
% will also help prevent bad behaviour due to thus-far unexplained nans in
% the trajectories...)
np_max = 100;
np_min = 10;
% Particles too far from sampling sites are not useful, so choose some
% maximum distance, dmax (km), from site at time of sampling event.
dmax = 50;

for iy = 1:nYears
    yStr = sprintf('y%i', years(iy));
    ne = length(unique(obs.(yStr).eventNum)); % number of sampling events
    for ie = 1:ne
        eStr = sprintf('event%i', ie);
        yd = traj.(yStr).(eStr).yearday; % event time
        % for each sampling event extract the best ranked (closest)
        % particles situated within the maximum distance from the sampling
        % site at time yd.
%         p_best = ismember(traj.(yStr).(eStr).rank, 1:np); % index np particles closest to event
%        p_near = traj.(yStr).(eStr).distance <= dmax;
%        p_extract = p_best & p_near;
        p_best = ismember(traj.(yStr).(eStr).rank, 1:np_max); % for now, just extract the 100 closest particles for each event, then consider distances later...
        p_extract = p_best;
        for iv = 1:nVars            
            vStr = labelVar{1,iv};
            
            % ONLY DIFFERENCE IN THIS CHANGING COLOURS SCRIPT IS NOT
            % FILTERING BY yd IN THE NEXT LINE
            
            traj.(yStr).(eStr).(vStr) = traj.(yStr).profiles.(vStr)(:,p_extract,depthLayer);
        end        
    end
end



% Use the tailedColourScale function, which uses the viridis colour
% palettes

colourMap = 'viridis'; % choose either 'viridis', 'plasma', 'magma' or 'inferno'
nCols = 100; % number of colours. 100 is the default of tailedColourScale
qr = [0.025 0.975]; % focus of colour scale lies within metric quantiles qr

% lo = min(traj.y2017.event1.temp(:));
% hi = max(traj.y2017.event1.temp(:));
% qq = quantile(traj.y2017.event1.temp(:), [0.025, 0.975]);
% hilo = qq(1); lohi = qq(2);
% cm = tailedColourScale(lo, hilo, lohi, hi);



% Set the colour scales before running the loop that generates plots.

% collect plotting variables and assign colour scale
% nCols = 100; % number of colours in scale
% cm = parula(nCols); % colour-map
% cm = viridis(nCols); % colour-map
% cm = inferno(nCols); % colour-map
% cm = magma(nCols); % colour-map
% cm = plasma(nCols); % colour-map

% There a couple of ways to set the colour scales -
% (1) Each plot could use the full range of colours, increasing
% contrast/interprebility of each plot, but changing scales between plots
% reducing their comparability.
% (2) Each year use the same colour scale across all samping events,
% potentially reducing contrast within separate plots, but allowing for
% easy comparison between sampling events.
% (3) Use the same colour scale across all sampling events and years,
% further reducing contrast within each individual plot, but allowing for
% easy comparison between all years and sampling events.
% - REMEMBER, in each case trajectories will be filtered by rank...


% Try creating all 3 colour schemes then using a switch in the plotting loop
% that controls which scheme is used...

rangesTot = nan(nVars,2);
for iy = 1:nYears
    yStr = sprintf('y%i', years(iy));
    ne = length(unique(obs.(yStr).eventNum)); % number of sampling events
    metricRanges.(yStr) = nan(2,ne,nVars);
    % colour scheme 1
    for ie = 1:ne
        eStr = sprintf('event%i', ie);
        for iv = 1:nVars
            vStr = labelVar{1,iv};
            % colour trajectories according to metric values
            
            % THE ONLY DIFFERENCE IN SETTING COLOURS FOR THIS CHANGING
            % COLOURS SCRIPT IS THE (:) IN THE NEXT 2 LINES
            rv = [min(traj.(yStr).(eStr).(vStr)(:)) ... % metric range
                max(traj.(yStr).(eStr).(vStr)(:))];
            metricRanges.(yStr)(:,ie,iv) = rv;
            traj.(yStr).(eStr).([vStr '_col1']) = 1 + round((nCols-1) * ...
                ((traj.(yStr).(eStr).(vStr) - rv(1)) ./ (rv(2) - rv(1))));
        end
    end
    % colour scheme 2
    % metric ranges over all events
    rv = [squeeze(min(metricRanges.(yStr)(1,:,:))) squeeze(max(metricRanges.(yStr)(2,:,:)))];
    for ie = 1:ne
        eStr = sprintf('event%i', ie);
        for iv = 1:nVars
            vStr = labelVar{1,iv};
            traj.(yStr).(eStr).([vStr '_col2']) = 1 + round((nCols-1) * ...
                ((traj.(yStr).(eStr).(vStr) - rv(iv,1)) ./ (rv(iv,2) - rv(iv,1))));
        end
    end
    % colour scheme 3
    % metric ranges over all events and years
    for iv = 1:nVars
        if ~isnan(rangesTot(iv,1)), rangesTot(iv,1) = min(rv(iv,1),rangesTot(iv,1));
        else, rangesTot(iv,1) = rv(iv,1); end
        if ~isnan(rangesTot(iv,2)), rangesTot(iv,2) = max(rv(iv,2),rangesTot(iv,2));
        else, rangesTot(iv,2) = rv(iv,2); end
    end    
    if iy == nYears
        for iyy = 1:nYears
            yyStr = sprintf('y%i', years(iyy));
            ne = length(unique(obs.(yyStr).eventNum)); % number of sampling events
            for ie = 1:ne
                eStr = sprintf('event%i', ie);
                for iv = 1:nVars
                    vStr = labelVar{1,iv};
                    traj.(yyStr).(eStr).([vStr 'col_3']) = 1 + round((nCols-1) * ...
                        ((traj.(yyStr).(eStr).(vStr) - rangesTot(iv,1)) ./ (rangesTot(iv,2) - rangesTot(iv,1))));
                end
            end
        end
        clear rangesTot
    end
end





% For each separate event make a multipanel plot showing the trajectories
% most relevant to the event, the separate panels should show temperature,
% salinity, nitrate and the ratio of large to small phytoplankton.

plotStartingLocations = 1;

transparency = 0.25;
% cm = cm(:,1:3); cm = [cm repmat(transparency, [size(cm,1),1])];

% Choose colour scheme 1, 2 or 3. Only 1 will work... for now
colourScheme = 1;

pointSize = 6; % Point size (marks station)
axisTextSize = 8;

for iy = 1:nYears
    % loop trough years    
    yStr = sprintf('y%i', years(iy));    
    ne = length(unique(obs.(yStr).eventNum)); % number of sampling events
    for ie = 1:ne
        % loop through events
        eStr = sprintf('event%i', ie);
        i_event = obs.(yStr).eventNum == ie;
        event_label = unique(obs.(yStr).event(i_event));
        event_label = event_label{1};
        event_date = datestr(unique(datenum(obs.(yStr).time(i_event,:))),'dd mmm yyyy');
        
        % Multipanel plots
        figure(); clf;
                
        nCol = ceil(sqrt(nVars)); nRow = ceil(nVars / nCol);
        rowH = 0.8 / nRow; colW = 0.85 / nCol;
        colX = 0.05 + linspace(0, 0.95, nCol+1); colX = colX(1:end-1);
        rowY = 0.1 + linspace(0.9, 0, nRow+1); rowY = rowY(2:end);

        for iv = 1:nVars
            % loop through plotting metrics
            plotVar = labelVar{1,iv};
            
            % index plotting panel
            rowID = ceil(iv / nCol);
            colID = iv - (rowID - 1) * nCol;
            
            axes('Position', [colX(colID), rowY(rowID), colW, rowH]);
            
            m_proj('lambert', 'clong', ctrLon, 'long', mapLon, 'lat', mapLat);
            
            if  ~exist(fullfile(topoDir,topoFile), 'file')
                res = topoFile(strfind(topoFile, '.mat') - 1);
                topo = m_gshhs(sprintf('%sc', res), 'patch', [.7, .7, .7], 'EdgeColor', 'none');
                if ~exist(topoDir, 'dir'), mkdir(topoDir); end
                m_gshhs(sprintf('%sc', res), 'save', fullfile(topoDir,topoFile));
            else
                m_usercoast(fullfile(topoDir,topoFile), 'patch', [.7, .7, .7], 'EdgeColor', 'none');
            end
%             m_grid('linewi', 2, 'FontSize', 24);
            m_grid('linewi', 2, 'FontSize', axisTextSize);
            
            hold on
            % plot trajectories
            best = ismember(traj.(yStr).(eStr).rank, 1:np_max); % most relevant trajectories for event ie
            lon = traj.(yStr).x(:,best);
            lat = traj.(yStr).y(:,best);
            % Colour scale
            
            col = traj.(yStr).(eStr).([plotVar '_col' num2str(colourScheme)]);
            
            % Use tailedColourScale - set tails using metric distribution
            % properties
            
            
            lo = min(traj.(yStr).(eStr).(plotVar)(:));
            hi = max(traj.(yStr).(eStr).(plotVar)(:));
            qq = quantile(traj.(yStr).(eStr).(plotVar)(:), [qr(1), qr(2)]);
            vmed = median(traj.(yStr).(eStr).(plotVar)(:));
            vstd = std(traj.(yStr).(eStr).(plotVar)(:));
            highTail = hi - qq(2) > vstd;
            lowTail = qq(1) - lo > vstd;
            
            if lowTail, hilo = qq(1); else, hilo = lo; end
            if highTail, lohi = qq(2); else, lohi = hi; end
            
            cm = tailedColourScale(lo, hilo, lohi, hi, colourMap);            
            cm = [cm repmat(transparency, [size(cm,1),1])];
            
            
            for ix = 1:np_max
                iMask = lon(:,ix) >= mapLon(1) & lon(:,ix) <= mapLon(2);
                missing = isnan(traj.(yStr).(eStr).(plotVar)(:,ix));
                if any(~iMask) || any(missing), continue, end
                % If there's no in-built to plot lines with changing
                % colours then we need to loop through time points...
                for it = 1:size(col,1)-1
                    col_it = round(sum(col(it:it+1,ix)) / 2);
                    m_plot(lon(it:it+1,ix), lat(it:it+1,ix), 'Linestyle', '-', ...
                        'Color', cm(col_it,:), 'LineWidth', 1);
                end
%                 m_plot(lon(:,ix), lat(:,ix), 'Linestyle', '-', 'Color', cm(col(ix),:), 'LineWidth', 1);
            end
            
            % plot trajectory starting locations
            if plotStartingLocations
                for ix = 1:np_max
                    iMask = lon(:,ix) >= mapLon(1) & lon(:,ix) <= mapLon(2);
                    missing = isnan(traj.(yStr).(eStr).(plotVar)(:,ix));
%                     missing = isnan(traj.(yStr).(eStr).(plotVar)(1,ix));
                    if any(~iMask) || any(missing), continue, end
%                     if any(~iMask) || missing, continue, end
                    m_plot(lon(1,ix), lat(1,ix), 'LineStyle', 'none', ...
                        'Marker', 'o', 'MarkerEdgeColor', 'none', ...
                        'MarkerFaceColor', brighten(cm(col(1,ix),1:3),-0.5), ...
                        'MarkerSize', ceil(pointSize/3));
                end
            end
            
            % Indicate location of sampling event            
            eventLat = unique(obs.(yStr).lat(i_event));
            eventLon = unique(obs.(yStr).lon(i_event));            
            m_plot(eventLon, eventLat, 'LineStyle', 'none', 'Marker', 'p', ...
                'MarkerSize', pointSize, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'r')
                        
            
            % Add colour bar - details depend on number of tails in colour map.
            cAxis = gca;
            % Colour bar position defined relative to axis
            cbP3 = 0.75 * cAxis.Position(3); % legend width
            cbP1 = cAxis.Position(1) - 0.5 * (cbP3 - cAxis.Position(3)); % legend x-position
            cbP4 = cAxis.Position(4) / 20; % legend height
            cbP2 = cAxis.Position(2); % legend y-position
            
            % Tick positions and labels
            nTicks = 5;
            if ~lowTail && ~highTail
                tickNames = linspace(lo,hi,nTicks);
                tickPositions = linspace(0,1,nTicks);
            elseif lowTail && ~highTail
                tickNames = [lo linspace(hilo,hi,nTicks-1)];
                tickPositions = [0 linspace(1/6,1,nTicks-1)];
            elseif ~lowTail && highTail
                tickNames = [linspace(lo,lohi,nTicks-1) hi];
                tickPositions = [linspace(0,1-1/6,nTicks-1) 1];
            elseif lowTail && highTail
                tickNames = [lo linspace(hilo,lohi,nTicks-2) hi];                
                tickPositions = [0 linspace(1/6,1-1/6,nTicks-2) 1];
            end
            
            if all(tickNames > 10)
                if min(diff(tickNames)) < 2
                    tickNames = round(tickNames,3,'significant');
                else
                    tickNames = round(tickNames,2,'significant');
                end
            else
                if min(diff(tickNames)) < 1.5
                    tickNames = round(tickNames,3,'significant');
                else
                    tickNames = round(tickNames,2,'significant');
                end
            end
            
            % Colour scale of colour bar
            
            if ~lowTail && ~highTail
                cAxis.Colormap = cm(:,1:3);
            elseif lowTail && ~highTail
                nc = size(cm,1);
                cmbar = nan(nc,3);
                ntail = floor(nc/6);
                nmain = nc - ntail;
                cmbar(ntail+1:end,:) = feval(colourMap,nmain);                
                lowCol = rgb2hsv(cmbar(ntail+1,:));
                hue = repmat(lowCol(1), [ntail 1]);
                sat = repmat(lowCol(2), [ntail 1]);
                val = linspace(4/5*lowCol(3), lowCol(3), ntail+1);
                val = val(1:end-1);
                cmbar(1:ntail,:) = hsv2rgb([hue(:) sat(:) val(:)]);
                cAxis.Colormap = cmbar(:,1:3);
            elseif ~lowTail && highTail
                nc = size(cm,1);
                cmbar = nan(nc,3);
                ntail = floor(nc/6);
                nmain = nc - ntail;
                cmbar(1:nmain,:) = feval(colourMap,nmain);                
                highCol = rgb2hsv(cmbar(nmain,:));
                hue = repmat(highCol(1), [ntail 1]);
                sat = repmat(highCol(2), [ntail 1]);
                val = linspace(highCol(3), 4/5*highCol(3), ntail+1);
                val = val(2:end);                
                cmbar(nmain+1:end,:) = hsv2rgb([hue(:) sat(:) val(:)]);
                cAxis.Colormap = cmbar(:,1:3);
            elseif lowTail && highTail
                nc = size(cm,1);
                cmbar = nan(nc,3);
                ntail = floor(nc/6);
                nmain = nc - 2*ntail;
                cmbar(ntail+1:end-ntail,:) = feval(colourMap,nmain);                                
                lowCol = rgb2hsv(cmbar(ntail+1,:));                
                hue = repmat(lowCol(1), [ntail 1]);
                sat = repmat(lowCol(2), [ntail 1]);
                val = linspace(4/5*lowCol(3), lowCol(3), ntail+1);
                val = val(1:end-1);
                cmbar(1:ntail,:) = hsv2rgb([hue(:) sat(:) val(:)]);                
                highCol = rgb2hsv(cmbar(end-ntail,:));
                hue = repmat(highCol(1), [ntail 1]);
                sat = repmat(highCol(2), [ntail 1]);
                val = linspace(highCol(3), 4/5*highCol(3), ntail+1);
                val = val(2:end);                
                cmbar(ntail+nmain+1:end,:) = hsv2rgb([hue(:) sat(:) val(:)]);
                cAxis.Colormap = cmbar(:,1:3);
            end
            
            cb = colorbar('southoutside','Ticks', tickPositions, 'TickLabels', tickNames, ...
                'Position', [cbP1 cbP2 cbP3 cbP4]);
            
            % label
            labVarFull = labelVar{2,strcmp(labelVar(1,:), plotVar)};
            labUnit = labelVar{3,strcmp(labelVar(1,:), plotVar)};            
            cb.Label.String = [labVarFull ' (' labUnit ')'];
            
%             hold off
        end
        
        % title axes and title.
        axes('Position', [0, 0.95, 1, 0.05] );
        set(gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None', 'TickLabelInterpreter', 'none');
        
        mTitle = ['AWI Hausgarten: ' event_date ' (sample ' event_label ')'];
        mTitle = strrep(mTitle, '_', '\_'); % change '_' to '\_' to prevent latex subscripts
        text(0.5, 0, mTitle, ...
            'FontSize', axisTextSize+4, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'Bottom')
        
        % save figure
        experiment = strrep(expName, 'YEAR', yStr);
        event = ['event' num2str(ie)];
        
        figName = ['trajectoriesMapPlot_' experiment '_' event '_colourChange.png'];        
%         figName = ['trajectoriesMapPlot_' experiment '_' event '.png'];
        
        figFile = fullfile(outDir, figName);
        
        fig = gcf;
        % Adjust figure window dimensions
        fig.Units = 'inches';
        fig.Position = [0 0 6 7];
        fig.PaperPositionMode = 'auto'; % save figure window as-is
        print(fig, figFile, '-r300', '-dpng');
        
        close(fig);
%         imcut(figFile);
    end
end




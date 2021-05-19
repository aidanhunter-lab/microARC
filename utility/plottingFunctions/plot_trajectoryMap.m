function plt = plot_trajectoryMap(Directories, Forc, varargin)

plt = figure;

extractStruct(Directories)
extractVarargin(varargin)

if ~exist('pointSize', 'var')
    pointSize = 10; % Point size of trajectory initial locations
end
if ~exist('axisTextSize', 'var')
    axisTextSize = 8;
end
if ~exist('alphaLine', 'var')
    alphaLine = 0.2; % transparancy of trajectory lines
end
if ~exist('projection', 'var')
    projection = 'lambert'; % map projection
end
if ~exist('EdgeColour', 'var')
    EdgeColour = 'none';
end
if ~exist('gridLineWidth', 'var')
    gridLineWidth = 2;
end
if ~exist('ytickSpace', 'var')
    ytickSpace = 5;
end
if ~exist('colAtlantic', 'var')
    colAtlantic = [1 0 0]; % colour of trajectory lines for Atlantic
end
if ~exist('colArctic', 'var')
    colArctic = [0 0 1]; % and Arctic
end
if ~exist('trajLineWidth', 'var')
    trajLineWidth = 1;
end
if ~exist('polygonLineWidth', 'var')
    polygonLineWidth = 2;
end
if ~exist('colPolygon', 'var')
    colPolygon = [0 0 0];
end
if ~exist('labelWaterOrigin', 'var')
    labelWaterOrigin = false; % include labels for trajectory groups
end
if ~exist('labelTextSize', 'var')
    labelTextSize = axisTextSize + 4;
end
if ~exist('labelFontWeight', 'var')
    labelFontWeight = 'bold';
end
if ~exist('labelSampleArea', 'var')
    labelSampleArea = false;    
end
if ~exist('plotTitle', 'var')
    plotTitle = [];
end
if ~exist('titleTextSize', 'var')
    titleTextSize = labelTextSize;
end
if ~exist('titleFontWeight', 'var')
    titleFontWeight = 'bold';
end


%% Coordinate ranges
ctrLon = 0; % Longitude of map centre
mapLon = [-max(abs(Forc.x(:))), max(abs(Forc.x(:)))];
mapLon = 0.02 .* diff(mapLon) .* [-1,1] + mapLon;
mapLat = [min(Forc.y(:)), max(abs(Forc.y(:)))];
mapLat = 0.06 .* diff(mapLat) .* [-1, 1] + mapLat;

%% Coastlines
m_proj(projection, 'clong', ctrLon, 'long', mapLon, 'lat', mapLat);
m_usercoast(fullfile(topoDir,topoName), 'patch', [.7, .7, .7], 'EdgeColor', EdgeColour);
m_grid('linewidth', gridLineWidth, 'fontsize', axisTextSize, 'ytick', ytickSpace);
hold on

%% Particle trajectories
for i = 1:Forc.nTraj
    lon = Forc.x(:,i);
    lat = Forc.y(:,i);
    type = Forc.waterMass(:,i);
    if strcmp(type, 'Atlantic'), Col = [colAtlantic alphaLine]; end
    if strcmp(type, 'Arctic'), Col = [colArctic alphaLine]; end
    m_plot(lon, lat, 'Linestyle', '-', ...
        'Color', Col, 'LineWidth', trajLineWidth);
end
% Highlight starting locations
for i = 1:Forc.nTraj
    lon = Forc.x(:,i);
    lat = Forc.y(:,i);
    type = Forc.waterMass(:,i);
    if strcmp(type, 'Atlantic'), Col = [colAtlantic alphaLine]; end
    if strcmp(type, 'Arctic'), Col = [colArctic alphaLine]; end
    m_plot(lon(1), lat(1), 'LineStyle', 'none', ...
        'Marker', 'o', 'MarkerEdgeColor', 'none', ...
        'MarkerFaceColor', brighten(Col(1:3),-0.5), ...
        'MarkerSize', ceil(pointSize/4));
end

%% Sample sites
if exist('Data', 'var') % Plotted only if 'Data' is passed as optional argument
    if ~exist('years', 'var')
        years = unique(Data.scalar.Year); % Years of data to include
    end
    ind = ismember(Data.scalar.Year, years);
    lon = Data.scalar.Longitude(ind);
    lat = Data.scalar.Latitude(ind);
    sampleCoords = table2struct(unique(table(lon, lat)), 'ToScalar', true);
    sampleArea = polygon_from_stations(sampleCoords, 25); % this 25 should be replaced by maxDist (the maximum distance of trajectories from the samples), which should be stored in Forc
    m_plot(sampleArea.lon, sampleArea.lat, 'Color', colPolygon, 'LineWidth', polygonLineWidth)
end
    
%% Labels
switch labelWaterOrigin, case true
    m_text(mapLon(1) + 0.1 * diff(mapLon), mapLat(2) - 0.4 * diff(mapLat), ...
        {'Arctic'; 'origin'}, 'Color', colArctic, ... 
        'FontSize', labelTextSize, 'FontWeight', labelFontWeight)
    m_text(mapLon(2) - 0.3 * diff(mapLon), mapLat(2) - 0.6 * diff(mapLat), ...
        {'Atlantic'; 'origin'}, 'Color', colAtlantic, ...
        'FontSize', labelTextSize, 'FontWeight', labelFontWeight)
end
switch labelSampleArea, case true
    m_text(mapLon(2) - 0.4 * diff(mapLon), mapLat(2) - 0.36 * diff(mapLat), ...
        {'sample'; 'area'}, 'Color', colPolygon, ... 
        'FontSize', labelTextSize, 'FontWeight', labelFontWeight)
end
if ~isempty(plotTitle) && ischar(plotTitle)
    m_text(mean(mapLon), mapLat(2) + 0.05 * diff(mapLat), plotTitle, ...
        'FontSize', titleTextSize, ...
        'FontWeight', titleFontWeight, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom');
end


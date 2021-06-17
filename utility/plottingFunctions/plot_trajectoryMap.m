function plt = plot_trajectoryMap(Directories, Forc, varargin)
% Plot map of Fram Strait overlayed with physical model trajectories

extractStruct(Directories)
extractVarargin(varargin)

if ~exist('newPlot', 'var')
    newPlot = true;
end

if newPlot
    plt = figure;
end

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

if ~exist('flowArrows', 'var')
    flowArrows = false;
end

if ~exist('includeLegend', 'var')
    includeLegend = true;
end
if ~exist('legendPosition', 'var')
    legendPosition = 'east';
end
if ~exist('legendTitle', 'var')
    legendTitle = 'Water origin';
end
if ~exist('legendTitleSize', 'var')
    legendTitleSize = 12;
end
if ~exist('legendTextSize', 'var')
    legendTextSize = 10;
end

if ~exist('stripedBorder', 'var')
    stripedBorder = false; % stipey polygon surrounds data area?
end
if ~exist('stripeColours', 'var')
    stripeColours = [0, 0, 0; 1, 1, 0]; % black and yellow
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
    if strcmp(type, 'Atlantic')
        Col = [colAtlantic alphaLine];
        pAtl =     m_plot(lon, lat, 'Linestyle', '-', ...
            'Color', Col, 'LineWidth', trajLineWidth);
    end
    if strcmp(type, 'Arctic')
        Col = [colArctic alphaLine];
        pArc = m_plot(lon, lat, 'Linestyle', '-', ...
            'Color', Col, 'LineWidth', trajLineWidth);
    end
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

%% Sample site area
if exist('Data', 'var') % Plotted only if 'Data' is passed as optional argument
    if ~exist('years', 'var')
        years = unique(Data.scalar.Year); % Years of data to include
    end
    ind = ismember(Data.scalar.Year, years);
    lon = Data.scalar.Longitude(ind);
    lat = Data.scalar.Latitude(ind);
    sampleCoords = table2struct(unique(table(lon, lat)), 'ToScalar', true);
    sampleArea = polygon_from_stations(sampleCoords, 25); % this 25 should be replaced by maxDist (the maximum distance of trajectories from the samples), which should be stored in Forc
    if ~stripedBorder
        m_plot(sampleArea.lon, sampleArea.lat, 'Color', colPolygon, 'LineWidth', polygonLineWidth)
    else
        ne = length(sampleArea.lon) - 1; % number of edges
        sampleArea_rad.lon = sampleArea.lon .* pi ./ 180; % convert to radians
        sampleArea_rad.lat = sampleArea.lat .* pi ./ 180;
        % Find the longest bounding edge
        edgeSize = nan(1, ne);
        for i = 1:ne
            lon = sampleArea_rad.lon(i:i+1);
            lat = sampleArea_rad.lat(i:i+1);
            edgeSize(i) = acos(sin(lat(1)) .* sin(lat(2)) + cos(lat(1)) .* cos(lat(2)) .* cos(abs(lon(2)-lon(1))));
        end
        longestEdge = edgeSize == max(edgeSize);
        % Coords of longest edge
        lon = sampleArea.lon(longestEdge | circshift(longestEdge, 1));
        lat = sampleArea.lat(longestEdge | circshift(longestEdge, 1));
        % Plot-scale coords
        lon_ = (lon - min(mapLon)) ./ diff(mapLon);
        lat_ = (lat - min(mapLat)) ./ diff(mapLat);
        el_ = (diff(lon_) .^ 2 + diff(lat_) .^ 2) .^ 0.5; % edge length
        stripeWidth_ = el_ / 10;
        spilloverfrac = 0;
        R = 6371;
        
        % Start with longest edge
        sampleArea.lon = unique(circshift(sampleArea.lon, -(find(longestEdge)-1)), 'stable');
        sampleArea.lat = unique(circshift(sampleArea.lat, -(find(longestEdge)-1)), 'stable');
        sampleArea.lon = [sampleArea.lon, sampleArea.lon(1)];
        sampleArea.lat = [sampleArea.lat, sampleArea.lat(1)];
        
        for i = 1:ne
            lon = sampleArea.lon(i:i+1);
            lat = sampleArea.lat(i:i+1);
            latIncrease = diff(lat) > 0;
            lon_ = (lon - min(mapLon)) ./ diff(mapLon);
            lat_ = (lat - min(mapLat)) ./ diff(mapLat);
            el_ = (diff(lon_) .^ 2 + diff(lat_) .^ 2) .^ 0.5; % edge length
            sn = el_ ./ stripeWidth_; % number of stripes
            % Need to rotate coords to get proper alignment
            p1 = [lon(1), lat(1)];
            p2 = [lon(2), lat(2)];
            p1_rad = pi ./ 180 .* p1;
            p2_rad = pi ./ 180 .* p2;
            p1_Cart = R .* [cos(p1_rad(2)) .* cos(p1_rad(1)), ...
                cos(p1_rad(2)) .* sin(p1_rad(1)), ...
                sin(p1_rad(2))];
            p2_Cart = R .* [cos(p2_rad(2)) .* cos(p2_rad(1)), ...
                cos(p2_rad(2)) .* sin(p2_rad(1)), ...
                sin(p2_rad(2))];
            % Rotate about z-axis
            ang = p1_rad(1); % longitude of 1st vertex
            K1 = [cos(ang), sin(ang), 0;
                -sin(ang), cos(ang), 0;
                0, 0, 1];
            p1_Cart = (K1 * p1_Cart(:))';
            p2_Cart = (K1 * p2_Cart(:))';
            % Rotate about y-axis
            ang = pi / 2 - p1_rad(2); % colatitude of 1st vertex
            if latIncrease
                ang = ang + pi; % place 1st vertex on south pole if it has lower latitude than 2nd vertex
            end
            K2 = [cos(ang), 0, sin(-ang);
                0, 1, 0;
                -sin(-ang), 0, cos(ang)];
            p1_Cart = (K2 * p1_Cart(:))';
            p2_Cart = (K2 * p2_Cart(:))';
            % Rotate about z-axis
            ang = atan(p2_Cart(2) / p2_Cart(1)); % longitude of transformed 2nd vertex
            K3 = [cos(ang), sin(ang), 0;
                -sin(ang), cos(ang), 0;
                0, 0, 1];
            p1_Cart = (K3 * p1_Cart(:))';
            p2_Cart = (K3 * p2_Cart(:))';
            % The transformed boundary edge is now aligned along
            % the prime meridian with 1st vertex at pole =>
            % curvature due to lat-longs should no longer exist...
            % Now create vectors of latitude values separated by
            % the stripe widths...
            % Convert our transformed Cartesian coords to
            % lat-longs, then split the latitudes, then transform
            % back to original coords.
            if latIncrease
                p1_rad = [0, -pi / 2]; % lat-longs (radians) of transformed vectors
            else
                p1_rad = [0, pi / 2];
            end
            p2_rad = [0, atan(p2_Cart(3) ./ p2_Cart(1))];
            latDiff = abs(p1_rad(2) - p2_rad(2)); % difference in transformed latitudes
            lat_inc = latDiff / (sn - spilloverfrac);
            if latIncrease
                lat_inc = -lat_inc;
            end
            lat_new = unique([p2_rad(2):lat_inc:p1_rad(2), p1_rad(2)]);
            if ~latIncrease
                lat_new = flip(lat_new);
            end
            lon_new = zeros(size(lat_new));
            % Inverse transform the new coords
            pnew_Cart = R .* [cos(lat_new(:)) .* cos(lon_new(:)), ...
                cos(lat_new(:)) .* sin(lon_new(:)), ...
                sin(lat_new(:))];
            pnew_Cart = (K1 \ (K2 \ (K3 \ pnew_Cart')))'; % original coordinate axes
            pnew_rad = [atan(pnew_Cart(:,2) ./ pnew_Cart(:,1)), ...
                atan(pnew_Cart(:,3) ./ (sum(pnew_Cart(:,1:2) .^ 2, 2) .^ 0.5))];
            pnew = pnew_rad .* 180 ./ pi;
            if i == 1
                clear sampleAreaNew
                sampleAreaNew.lon = pnew(:,1);
                sampleAreaNew.lat = pnew(:,2);
                stripeColours_ = repmat(stripeColours, [size(pnew, 1), 1]);
                sampleAreaNew.colour = stripeColours_(1:size(pnew, 1),:);
            else
                sampleAreaNew.lon = [sampleAreaNew.lon; pnew(2:end,1)];
                sampleAreaNew.lat = [sampleAreaNew.lat; pnew(2:end,2)];
                stripeColours_ = repmat(stripeColours, [size(pnew, 1), 1]);
                if spilloverfrac > 0
                    if ~all(stripeColours_(1,:) == sampleAreaNew.colour(end,:))
                        stripeColours_ = stripeColours_(2:end,:);
                    end
                else
                    if all(stripeColours_(1,:) == sampleAreaNew.colour(end,:))
                        stripeColours_ = stripeColours_(2:end,:);
                    end
                end
                sampleAreaNew.colour = [sampleAreaNew.colour; stripeColours_(1:size(pnew, 1)-1,:)];
            end
            spilloverfrac = 1 - (sn - floor(sn)); % proportion of last stripe in next line
        end
        np = length(sampleAreaNew.lon);
        for i = 1:np-1
            m_plot(sampleAreaNew.lon(i:i+1), sampleAreaNew.lat(i:i+1), 'Color', sampleAreaNew.colour(i,:), 'LineWidth', polygonLineWidth)
        end
    end
end


%% Arrows
switch flowArrows, case true
    % Build arrows manually -- automatic code may be tricky...
    n = 100;
    lon = linspace()
    
    lat = linspace(68, 80, n)
    
m_vec



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

switch includeLegend, case true
    pAtl.Color = pAtl.Color(1:3);
    pArc.Color = pArc.Color(1:3);
    leg = legend([pArc,pAtl], 'Arctic', 'Atlantic', 'Location', legendPosition);
    leg.FontSize = legendTextSize;
    leg.Title.String = legendTitle;
    leg.Title.FontSize = legendTitleSize;
end
    
    


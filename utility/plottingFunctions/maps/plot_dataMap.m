function plt = plot_dataMap(Directories, Data, varargin)
% Plot map of Fram Strait overlayed with in-situ sample locations

extractStruct(Directories)
extractVarargin(varargin)

if ~exist('newPlot', 'var')
    newPlot = true;
end

if newPlot
    plt = figure;
end

if ~exist('projection', 'var')
    projection = 'miller';
end
if ~exist('lonSpace', 'var')
    lonSpace = 0.3; % extra map space around data
end
if ~exist('latSpace', 'var')
    latSpace = 0.3;
end
if ~exist('EdgeColour', 'var')
    EdgeColour = 'none';
end
if ~exist('gridLineWidth', 'var')
    gridLineWidth = 2;
end
if ~exist('axisTextSize', 'var')
    axisTextSize = 8;
end
if ~exist('ytickSpace', 'var')
    ytickSpace = 5;
end

if ~exist('Year', 'var')
    Year = max(Data.scalar.Year);
end
if ~exist('pointSize', 'var')
    pointSize = 10;
end

if ~exist('Marker', 'var')
    Marker = 'o';
end

if ~exist('showWaterOrigin', 'var')
    showWaterOrigin = false;
end
if ~exist('polygonAlpha', 'var')
    polygonAlpha = 0.5;
end

if ~exist('polygonExpand', 'var')
    polygonExpand = 0; % extra days around sample points
end

if ~exist('polygonSmooth', 'var')
    polygonSmooth = false;
end

if ~exist('colourByDataType', 'var')
    colourByDataType = false;
end

if ~exist('pieSize', 'var')
    pieSize = 0.02; % radius of data pies relative to plot scale
end

if ~exist('includeLegend', 'var')
    includeLegend = true;
end

if ~exist('legendPosition', 'var')
    legendPosition = 'southeast';
end

if ~exist('omitMapGrid', 'var')
    omitMapGrid = false;
end

if ~exist('colSat', 'var')
    colSat = [];
end

if ~exist('fullAreaPolygon', 'var')
    fullAreaPolygon = false;
end
if ~exist('colPolygon', 'var')
    colPolygon = [0 0 0];
end
if ~exist('polygonLineWidth', 'var')
    polygonLineWidth = 2;
end

if ~exist('trimPolygons', 'var')
    trimPolygons = false;
end

if ~exist('legendTitle', 'var')
    legendTitle = 'Data';
end
if ~exist('legendTitleSize', 'var')
    legendTitleSize = 12;
end
if ~exist('legendTextSize', 'var')
    legendTextSize = 10;
end

if ~exist('stripedBorder', 'var')
    stripedBorder = false;
end
if ~exist('stripeColours', 'var')
    stripeColours = [0, 0, 0; 1, 1, 0]; % black and yellow
end


%% Coordinate ranges
ctrLon = 0; % Longitude of map centre
if ~exist('mapLon', 'var')
    mapLon = [min(Data.scalar.Longitude(:)), ...
        max(Data.scalar.Longitude(:))];
    mapLon = lonSpace .* diff(mapLon) .* [-1, 1] + mapLon;
%     mapLon = round(mapLon, 2, 'significant');
end
if ~exist('mapLat', 'var')
    mapLat = [min(Data.scalar.Latitude(:)), ...
        max(Data.scalar.Latitude(:))];
    mapLat = latSpace .* diff(mapLat) .* [-1, 1] + mapLat;
%     mapLat = round(mapLat, 2, 'significant');
end
if ~exist('ctrLon', 'var')
    ctrLon = mean(mapLon);
end
lonRange = range(mapLon);
latRange = range(mapLat);


% Colours for data sampled from different water masses
if ~exist('colArctic', 'var')
    colArctic = [0 0 1];
end
if ~exist('colAtlantic', 'var')
    colAtlantic = [1 0 0];
end
if ~exist('colMix', 'var')
    colMix = [1 0 1];
end
waterMass = unique(Data.scalar.waterMass);
nWaterMass = length(waterMass);
col = cell(nWaterMass, 1);
for i = 1:length(waterMass)
    switch waterMass{i}
        case 'Atlantic'
            col{i,1} = colAtlantic;
        case 'Arctic'
            col{i,1} = colArctic;
        case 'Arctic/Atlantic'
            col{i,1} = colMix;
    end
end
colTable = table(waterMass, col);
if ~isempty(colSat) && isnumeric(colSat)
    if 1 < colSat(1) || colSat(1) <= 0
        colSat = 1;
        warning('Colour saturation "colSat" must be a scalar in (0,1] -- value has been reset to 1.')
    end
    for i = 1:height(colTable)
        col_ = rgb2hsv(colTable.col{i});
        col_(2) = colSat;
        colTable.col{i} = hsv2rgb(col_);
    end
end


%% Coastlines
switch omitMapGrid
    case false
        m_proj(projection, 'clong', ctrLon, 'long', mapLon, 'lat', mapLat);
        m_usercoast(fullfile(topoDir,topoName), 'patch', [.7, .7, .7], 'EdgeColor', EdgeColour);
        m_grid('linewidth', gridLineWidth, 'fontsize', axisTextSize, 'ytick', ytickSpace);
        hold on
    case true
        m_proj(projection, 'clong', ctrLon, 'long', mapLon, 'lat', mapLat);
        m_grid('linestyle', 'none', 'xticklabels', [], 'yticklabels', [], ... 
            'ticklen', 0, 'box', 'none', 'backgroundcolor', 'none'); %, 'color', 'none');
        hold on
end

%% Water origin
% Use physical model trajectories to make polygons indicating water origin.

switch showWaterOrigin, case true
    if ~exist('Forc', 'var')
        warning('Water mass origins are not highlighted! If showWaterOrigin=true then "Forc" must be included within "varargin" so that plot_dataMap.m can access the physical model trajectories.')
    else
        waterMasses = unique(Forc.waterMass);
        
        switch trimPolygons
            
            case false
                
                for i = 1:length(waterMasses)
                    waterMass = waterMasses{i};
                    ind = strcmp(Forc.waterMass, waterMass); % index trajectories originating from Atlantic or Arctic
                    ind = repmat(ind, [size(Forc.t, 1), 1]);
                    % trim to plot boundaries
                    ind = ind & ...
                        mapLon(1) <= Forc.x & Forc.x <= mapLon(2) & ...
                        mapLat(1) <= Forc.y & Forc.y <= mapLat(2);
                    % trim to data sampling times
                    ind = ind & ...
                        min(Data.scalar.Yearday) - polygonExpand <= Forc.Yearday & ...
                        Forc.Yearday <= max(Data.scalar.Yearday) + polygonExpand;
                    % coordinates of all relevant trajectory points
                    coords.lon = Forc.x(ind);
                    coords.lat = Forc.y(ind);
                    col = colTable.col{strcmp(colTable.waterMass, waterMass)};
                    switch polygonSmooth
                        case true
                            edges = polygon_from_stations(coords, polygonExpand);
                            patch = m_patch(edges.lon, edges.lat, col);
                        case false
                            edges = boundary(coords.lon, coords.lat);
                            patch = m_patch(coords.lon(edges), coords.lat(edges), col);
                    end
                    patch.FaceAlpha = polygonAlpha;
                end
                
            case true
                
                ind = Data.scalar.Year == Year;
                dat_ = unique(table(Data.scalar.Longitude(ind), Data.scalar.Latitude(ind)));
                dat_.Properties.VariableNames = {'lon', 'lat'};
                sampleArea = polygon_from_stations(dat_, 0);
                for i = 1:length(waterMasses)
                    clear coords
                    waterMass = waterMasses{i};
                    ind = strcmp(Forc.waterMass, waterMass); % index trajectories originating from Atlantic or Arctic
                    ind = repmat(ind, [size(Forc.t, 1), 1]);
                    % trim to plot boundaries
                    ind = ind & ...
                        mapLon(1) <= Forc.x & Forc.x <= mapLon(2) & ...
                        mapLat(1) <= Forc.y & Forc.y <= mapLat(2);
                    % trim to data sampling times
                    ind = ind & ...
                        min(Data.scalar.Yearday) - polygonExpand <= Forc.Yearday & ...
                        Forc.Yearday <= max(Data.scalar.Yearday) + polygonExpand;
                    % coordinates of all relevant trajectory points
                    coords.lon = Forc.x(ind);
                    coords.lat = Forc.y(ind);
                    % omit points outside exterior boundary
                    ind = inpolygon(coords.lon, coords.lat, sampleArea.lon, sampleArea.lat);
                    coords.lon = coords.lon(ind);
                    coords.lat = coords.lat(ind);
                    
                    edges = boundary(coords.lon, coords.lat);
                    edgeCoords.(waterMass).lon = coords.lon(edges);
                    edgeCoords.(waterMass).lat = coords.lat(edges);
                end
                
                % edgeCoords contains lat-longs for all vertices defining
                % polygons representing water from Arctic and Atlantic.
                % We want to isolate the vertices separating the
                % water masses...
                coordsArc = [edgeCoords.Arctic.lon, edgeCoords.Arctic.lat];
                coordsAtl = [edgeCoords.Atlantic.lon, edgeCoords.Atlantic.lat];

                % As edgeCoords specify polygons, the coordinates are
                % sequential and form circuits around the shape. The main
                % split between water masses within the sample area is
                % longitudinal (east-west). If choose one of the polygons
                % and select its west-most vertex then there's two routes
                % to the east-most vertex. Extract vertices from the route
                % passing the water mass boundary. This will work better
                % with the Arctic polygon because it spans more of the
                % total width...
                coords = coordsArc;
                minLon = min(coords(:,1));
                maxLon = max(coords(:,1));
                iminLon = find(coords(:,1) == minLon); % index west-most and
                imaxLon = find(coords(:,1) == maxLon); % east-most vertices
                % Get indices of vertices from (1) forwards and (2)
                % backwards routes from west to east (end points included)
                if iminLon < imaxLon
                    i1 = iminLon:imaxLon;
                    i2 = [iminLon:-1:1, size(coords, 1):-1:imaxLon];
                else
                    i1 = [iminLon:size(coords, 1), 1:imaxLon];
                    i2 = iminLon:-1:imaxLon;
                end
                % The Arctic polygon was chosen => vertices on the water
                % mass boundary will occur on the west-east route with
                % lowest latitude.
                if max(coords(i1,2)) < max(coords(i2,2))
                    route = i1;
                else
                    route = i2;
                end
                % Extract all vertices along lower side of Arctic polygon
                vertices = coords(route,:);
                % Now use these points to construct 2 polygons that fit
                % neatly inside the full bounding polygon. If the west- and
                % east-most vertices are connected to the nearest point on
                % the surrounding polygon then we should have all we
                % need...
                
                coordsBound = [sampleArea.lon(:), sampleArea.lat(:)]; % vertex lat-longs of surrounding polygon
                westPoint = vertices(vertices(:,1) == min(vertices(:,1)),:);
                eastPoint = vertices(vertices(:,1) == max(vertices(:,1)),:);

                % Find which bounding line is nearest to westPoint and
                % eastPoint. Loop through bounding-polygon sides
                % calculating distances -- store minimum distance point.
                
                % Need to use spherical coordinates to do this properly.
                % Otherwise Earth's surface curvature cause problems when
                % using Euclidean geometry... Each bounding edge is
                % defined by 2 points, P1 and P2. Rotate axes such that P1
                % and P2 both lie on prime meridian (zero longitude) and
                % such that point Q (which is either westPoint or
                % eastPoint) lies on equator (zero latitude). Then
                % point Y (r=r, theta=0, phi=0) is the position on arc 
                % connecting P1 to P2 that is closest to Q. To find Y in
                % original coords use inverse transform...
                
                intersectPointWest = nan(size(coordsBound, 1)-1, 2);
                intersectDistWest = inf(size(intersectPointWest, 1), 1);
                intersectPointEast = nan(size(coordsBound, 1)-1, 2);
                intersectDistEast = inf(size(intersectPointEast, 1), 1);

                
                R = 6371; % approx Earth radius (km)

                for i = 1:size(coordsBound, 1) - 1
                    
                    % Convert lat-longs to Cartesian coords
                    coordsBound_Cart = R .* [ ...
                        cos(pi ./ 180 .* coordsBound(:,2)) .* ...
                        cos(pi ./ 180 .* coordsBound(:,1)), ... % x-coord
                        cos(pi ./ 180 .* coordsBound(:,2)) .* ...
                        sin(pi ./ 180 .* coordsBound(:,1)), ... % y-coord
                        sin(pi ./ 180 .* coordsBound(:,2))]; % z-coord
                    eastPoint_Cart = R .* [ ...
                        cos(pi ./ 180 .* eastPoint(2)) .* ...
                        cos(pi ./ 180 .* eastPoint(1)), ...
                        cos(pi ./ 180 .* eastPoint(2)) .* ...
                        sin(pi ./ 180 .* eastPoint(1)), ...
                        sin(pi ./ 180 .* eastPoint(2))];
                    westPoint_Cart = R .* [ ...
                        cos(pi ./ 180 .* westPoint(2)) .* ...
                        cos(pi ./ 180 .* westPoint(1)), ...
                        cos(pi ./ 180 .* westPoint(2)) .* ...
                        sin(pi ./ 180 .* westPoint(1)), ...
                        sin(pi ./ 180 .* westPoint(2))];

                    edgeVertices = coordsBound(i:i+1,:); % lat-lon coords of bounding edge vertices
                    % sort latitude into descending order
                    [~,o] = sort(edgeVertices(:,2), 'descend');
                    edgeVertices = edgeVertices(o,:);
                    edgeVertices_Cart = coordsBound_Cart(i:i+1,:); % xyz coords of bounding edge vertices
                    edgeVertices_Cart = edgeVertices_Cart(o,:);
                    
                    % Series of rotations...
                    % 1: about z-axis
                    ang = pi ./ 180 .* edgeVertices(1,1); % rotation angle (longitude of 1st vertex of bounding side, in radians)
                    K1 = [cos(ang), sin(ang), 0; ...
                        -sin(ang), cos(ang), 0; ...
                        0, 0, 1]; % rotation matrix
                    edgeVertices_Cart = (K1 * edgeVertices_Cart')';
                    westPoint_Cart =  (K1 * westPoint_Cart')';
                    eastPoint_Cart =  (K1 * eastPoint_Cart')';
                    % 2: about y-axis
                    ang = pi / 2 - pi ./ 180 .* edgeVertices(1,2); % rotation angle (colatitude of 1st vertex of bounding side, in radians)
                    K2 = [cos(ang), 0, sin(-ang); ...
                        0, 1, 0;
                        -sin(-ang), 0, cos(ang)];
                    edgeVertices_Cart = (K2 * edgeVertices_Cart')';
                    westPoint_Cart = (K2 * westPoint_Cart')';
                    eastPoint_Cart = (K2 * eastPoint_Cart')';
                    % 3: about z again, now that a point lies on north pole
                    xy = edgeVertices_Cart(2,1:2);
                    arg = xy(2) ./ abs(xy(1));
                    ang = atan(arg); % longitude  of 2nd vertex of bounding side, in radians
                    if xy(1) < 0
                        if xy(2) > 0
                            ang = ang + 90;
                        end
                        if xy(2) < 0
                            ang = -(90 + ang);
                        end
                    end
%                     ang = atan(edgeVertices_Cart(2,2) ./ edgeVertices_Cart(2,1)); % longitude  of 2nd vertex of bounding side, in radians
                    K3 = [cos(ang), sin(ang), 0; ...
                        -sin(ang), cos(ang), 0; ...
                        0, 0, 1];
                    edgeVertices_Cart = (K3 * edgeVertices_Cart')';
                    westPoint_Cart =  (K3 * westPoint_Cart')';
                    eastPoint_Cart =  (K3 * eastPoint_Cart')';
                    % Now the bounding vertices lie on x-z plane with the
                    % 1st vertex on z-axis.
                    % 4: a final rotation about y-axis places target point on (x,y,z) = (R,0,0)
                    
                    arg = westPoint_Cart(3) ./ sum(westPoint_Cart(1:2) .^ 2) .^ 0.5;
                    angWest = atan(arg);
                    if westPoint_Cart(1) < 0
                        angWest = pi - angWest;
                    end
                    arg = eastPoint_Cart(3) ./ sum(eastPoint_Cart(1:2) .^ 2) .^ 0.5;
                    angEast = atan(arg);
                    if eastPoint_Cart(1) < 0
                        angEast = pi - angEast;
                    end
                    

%                     angWest = atan(westPoint_Cart(3) ./ sum(westPoint_Cart(1:2) .^ 2) .^ 0.5);
%                     angEast = atan(eastPoint_Cart(3) ./ sum(eastPoint_Cart(1:2) .^ 2) .^ 0.5);
%                     if westPoint_Cart(1) < 0
%                         angWest = angWest +  (pi/2 - angWest);
%                     end
%                     if eastPoint_Cart(1) < 0
%                         angEast = angEast +  (pi/2 - angEast);
%                     end
                    KWest = [cos(angWest), 0, sin(angWest); ...
                        0, 1, 0;
                        -sin(angWest), 0, cos(angWest)];
                    KEast = [cos(angEast), 0, sin(angEast); ...
                        0, 1, 0;
                        -sin(angEast), 0, cos(angEast)];
                    
                    % Test validity of edge i for west and east point --
                    % fully transformed z-coords of bound should straddle 0
                    westPoint_Cart = (KWest * westPoint_Cart')';
                    eastPoint_Cart = (KEast * eastPoint_Cart')';
                    
                    
                    validWest = (KWest * edgeVertices_Cart')';
                    validWest = validWest(1,3) >= 0 && validWest(2,3) <= 0 && ...
                        validWest(2,3) < westPoint_Cart(3) && westPoint_Cart(3) < validWest(1,3);
                    validEast = (KEast * edgeVertices_Cart')';
                    validEast = validEast(1,3) >= 0 && validEast(2,3) <= 0 && ...
                        validEast(2,3) < eastPoint_Cart(3) && eastPoint_Cart(3) < validEast(1,3);

                    if ~validWest && ~validEast
                        continue
                    else
                        intersectWest = [R, 0, 0]';
                        intersectEast = [R, 0, 0]';
                        intersectWest_Cart = inv(K1) * (inv(K2) * (inv(K3) * (inv(KWest) * intersectWest)));
                        intersectEast_Cart = inv(K1) * (inv(K2) * (inv(K3) * (inv(KEast) * intersectEast)));
                        
                        % Convert to lat-longs
                        intersectEast = [sum(intersectEast_Cart .^ 2) .^ 0.5, ...
                            atan(intersectEast_Cart(2) ./ intersectEast_Cart(1)) .* 180 ./ pi, ...
                            atan(intersectEast_Cart(3) ./ sum(intersectEast_Cart(1:2) .^ 2) .^ 0.5) .* 180 ./ pi];
                        intersectWest = [sum(intersectWest_Cart .^ 2) .^ 0.5, ...
                            atan(intersectWest_Cart(2) ./ intersectWest_Cart(1)) .* 180 ./ pi, ...
                            atan(intersectWest_Cart(3) ./ sum(intersectWest_Cart(1:2) .^ 2) .^ 0.5) .* 180 ./ pi];
                        intersectWest = intersectWest(2:3); % remove 1st element -- which should equal R
                        intersectEast = intersectEast(2:3); % to leave only lat-longs
                        
                        intersectPointWest(i,:) = intersectWest;
                        intersectPointEast(i,:) = intersectEast;
                        
                        lon1 = westPoint(1) .* pi ./ 180;
                        lat1 = westPoint(2) .* pi ./ 180;
                        lon2 = intersectWest(1) .* pi ./ 180;
                        lat2 = intersectWest(2) .* pi ./ 180;
                        ang = acos(sin(lat1) .* sin(lat2) + cos(lat1) .* cos(lat2) .* cos(0.5 .* abs(lon2-lon1)));
                        intersectDistWest(i) = R .* ang;
                        
                        lon1 = eastPoint(1) .* pi ./ 180;
                        lat1 = eastPoint(2) .* pi ./ 180;
                        lon2 = intersectEast(1) .* pi ./ 180;
                        lat2 = intersectEast(2) .* pi ./ 180;
                        ang = acos(sin(lat1) .* sin(lat2) + cos(lat1) .* cos(lat2) .* cos(0.5 .* abs(lon2-lon1)));
                        intersectDistEast(i) = R .* ang;
                    end
                    
                end
                    
                wb = find(intersectDistWest == min(intersectDistWest)); % index of nearest bounding line
                eb = find(intersectDistEast == min(intersectDistEast));
                intersectPointWest = intersectPointWest(wb,:);
                intersectPointEast = intersectPointEast(eb,:);

                                % Include boundary values into vertices
                vertices = [intersectPointWest; vertices; intersectPointEast];
                [~,o] = sort(vertices(:,1)); % sort: west to east
                vertices = vertices(o,:);
                % Now all the info we need is stored in vertices and
                % coordsBound... define 2 polygons that neatly fit into the
                % bounding polygon... both polygons include 'vertices', but
                % use different bounding points.
                % Insert the intersecting points into the bounding
                % vertices, then trace out 2 routes...
                if eb < wb
                    coordsBound_ = [coordsBound(1:eb,:); ... 
                        intersectPointEast; coordsBound(eb+1:wb,:); ... 
                        intersectPointWest; coordsBound(wb+1:end,:)];
                else
                    coordsBound_ = [coordsBound(1:wb,:); ... 
                        intersectPointWest; coordsBound(wb+1:eb,:); ... 
                        intersectPointEast; coordsBound(eb+1:end,:)];
                end
                % split coordsBound_ into 2 parts, as routes from west to east
                wb = find(all(intersectPointWest == coordsBound_, 2));
                eb = find(all(intersectPointEast == coordsBound_, 2));
                nc = size(coordsBound_, 1);
                
                if eb < wb
                    i1 = [wb:nc, 1:eb];
                    i2 = wb:-1:eb;
                else
                    i1 = wb:eb;
                    i2 = [wb:-1:1, nc:-1:eb];
                end
                route1 = coordsBound_(i1,:);
                route2 = coordsBound_(i2,:);                
                if max(route1(:,2)) > max(route2(:,2))
                    routeArc = route1;
                    routeAtl = route2;
                else
                    routeArc = route2;
                    routeAtl = route1;
                end
                
                % Combine bounding values with water-mass intersection
                routeArc_ = [routeArc; flip(vertices(1:end-1,:))];
                routeAtl_ = [routeAtl; flip(vertices(1:end-1,:))];
                
                % These routes define 2 polygons that should fit neatly
                % into the bounding polygon that defines total sample area.
                % The 2 polygons are coloured according to water origin.
                col = colTable.col{strcmp(colTable.waterMass, 'Arctic')};
                m_patch(routeArc_(:,1), routeArc_(:,2), col);
                col = colTable.col{strcmp(colTable.waterMass, 'Atlantic')};
                m_patch(routeAtl_(:,1), routeAtl_(:,2), col);
        end
    end
end


%% Sample locations

switch colourByDataType

    case false

        ind = Data.scalar.Year == Year;
        dat = unique(table(Data.scalar.Event(ind), Data.scalar.Latitude(ind), ...
            Data.scalar.Longitude(ind)));
        dat = [dat, Data.scalar.waterMass];
        dat.Properties.VariableNames = {'Event', 'Latitude', 'Longitude', 'waterMass'};
        waterMasses = unique(dat.waterMass);
        for j = 1:length(waterMasses)
            waterMass = waterMasses{j};
            ind = strcmp(dat.waterMass, waterMass);
            lon = dat.Longitude(ind);
            lat = dat.Latitude(ind);
            col = colTable.col{strcmp(colTable.waterMass, waterMass)};
            m_plot(lon, lat, 'LineStyle', 'none', 'Marker', Marker, ...
                'MarkerSize', pointSize, 'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', col);
            if j == 1, hold on; end
            if j == length(waterMasses), hold off; end
        end
        
    case true

        ind = Data.scalar.Year == Year;
        dat = unique(table(Data.scalar.Event(ind), Data.scalar.Latitude(ind), ...
            Data.scalar.Longitude(ind)));
        dat = [dat, Data.scalar.waterMass];
        dat.Properties.VariableNames = {'Event', 'Latitude', 'Longitude', 'waterMass'};
        % find which data types are available from each sample event
        OM = {'PON','POC','chl_a'}; % index POM and chl as single group because these measurements derive from single data set
        vars = Data.scalar.Variable;
        vars(ismember(vars, OM)) = {'OM'};
        inCostFunction = Data.scalar.inCostFunction;
        vars = arrayfun(@(x) unique(vars(Data.scalar.Event == x & ...
            inCostFunction)), dat.Event, 'UniformOutput', false);
        dat.dataTypes = vars;
        sizeSamples = unique(Data.sizeFull.Event);
        for i = 1:length(sizeSamples)
            ev = sizeSamples(i);
            dat.dataTypes{dat.Event == ev} = [dat.dataTypes{dat.Event == ev}; 'size'];
        end
        % reorder dat by number of data types so that the busiest points
        % are plotted last and layered on top
        nd = cellfun(@(x) length(x), dat.dataTypes);
        [~,o] = sort(nd);
        dat = dat(o,:);
        coords = dat(:,{'Longitude','Latitude'});
        
        switch fullAreaPolygon, case true
            
            if ~stripedBorder
                
                coords_ = coords;
                coords_.Properties.VariableNames = {'lon','lat'};
                sampleArea = polygon_from_stations(coords_, 0);
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
                    
                    
%                     lat_inc = latDiff / sn;                    
%                     lat_new =  p1_rad(2):-lat_inc:p2_rad(2);                    
%                     if sn ~= round(sn)
%                         lat_new = [lat_new, p2_rad(2)];
%                     end
%                     lon_new = zeros(size(lat_new));
%                     
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
        
        % colours for each data type
        colOM = [0 1 0];
        colN = [1 1 0];
        colSize = [0 1 1];
        colData = table({'OM'; 'N'; 'size'}, [colOM; colN; colSize]);
        colData.Properties.VariableNames = {'data', 'col'};
        % for each sample event plot points as pies sliced according to which data
        % types were sampled
        tt = cell(height(dat), 1); % pie divisions
        for i = 1:height(dat)
            nd = length(dat.dataTypes{i});
            tt{i} = (1:nd) ./ nd .* 2 .* pi;
        end
        for i = 1:height(dat)
            nd = length(tt{i});
            for j = 1:nd
                dv = dat.dataTypes{i}(j);
                % start/end points of arcs for pie outer cirlces
                if j > 1
                    t(1) = tt{i}(j-1);
                else
                    t(1) = 0;
                end
                t(2) = tt{i}(j);
                % approx arc from t(1) to t(2)
                a = linspace(t(1), t(2), 100);
                % convert from lat-long coords to standardised coords
                % defined using plot dimensions -- this prevents distortion
                lon = coords.Longitude(i);
                lat = coords.Latitude(i);
                lon_ = (lon - min(mapLon)) ./ lonRange;
                lat_ = (lat - min(mapLat)) ./ latRange;
                % create slice
                if nd == 1
                    x_ = lon_ + pieSize .* [cos(t(2)); cos(t(1)); cos(a(:))];
                    y_ = lat_ + pieSize .* [sin(t(2)); sin(t(1)); sin(a(:))];
                else
                    x_ = lon_ + pieSize .* [cos(t(2)); 0; cos(t(1)); cos(a(:))];
                    y_ = lat_ + pieSize .* [sin(t(2)); 0; sin(t(1)); sin(a(:))];
                end
                x = min(mapLon) + x_ .* lonRange;
                y = min(mapLat) + y_ .* latRange;
                col = colData.col(ismember(colData.data, dv),:);
                p = m_patch(x, y, col);
            end
        end
                
        switch includeLegend, case true
            EdgeColour = p.EdgeColor;
            p1 = m_scatter(0,0, 'MarkerEdgeColor', EdgeColour, ...
                'MarkerFaceColor', colData.col(strcmp(colData.data, 'N'),:));
            p2 = m_scatter(0,0, 'MarkerEdgeColor', EdgeColour, ...
                'MarkerFaceColor', colData.col(strcmp(colData.data, 'OM'),:));
            p3 = m_scatter(0,0, 'MarkerEdgeColor', EdgeColour, ...
                'MarkerFaceColor', colData.col(strcmp(colData.data, 'size'),:));
            leg = legend([p1,p2,p3], 'DIN', 'OM', 'size', 'location', legendPosition);
            leg.Title.String = legendTitle;
            leg.Title.FontSize  = legendTitleSize;
            leg.FontSize = legendTextSize;
        end
end

% 
% mTitle = 'Hausgarten sample sites';
% 
% m_text(mean(mapLon), mapLat(2) + 0.1 * diff(mapLat), mTitle, ...
%     'FontSize', axisTextSize+4, 'FontWeight', 'bold', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom');
% 
% m_text(mapLon(1) + 0.2 * diff(mapLon), mapLat(2) - 0.2 * diff(mapLat), ...
%     '2016', 'Color', col2016, 'FontSize', axisTextSize + 4, 'FontWeight', 'bold')
% m_text(mapLon(1) + 0.21 * diff(mapLon), mapLat(2) - 0.3 * diff(mapLat), ...
%     '2017', 'Color', col2017, 'FontSize', axisTextSize + 4, 'FontWeight', 'bold')
% 
% m_text(19, 78.75, 'Spitsbergen', 'FontSize', axisTextSize - 2, 'FontWeight', 'bold', ... 
%     'HorizontalAlignment', 'right')
% 
% % save figure
% figName = 'sampleSites';
% outDir = '/home/aidan/Documents/work/microARC/slides';
% figFile = fullfile(outDir, figName);
% 
% % Adjust figure window dimensions
% fig1.Units = 'inches';
% fig1.Position = [0 0 6 3];
% fig1.PaperPositionMode = 'auto'; % save figure window as-is
% 
% switch savePlot, case true, print(fig1, figFile, '-r300', '-dpng'); end
% 

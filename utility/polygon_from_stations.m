function poly = polygon_from_stations(stations, radius)
    % poly = polygon_from_stations(stations, radius)
    % =====================================================================
    % FUNCTION: Determine polygon encompassing all 'stations' and extended
    %           in each direction by 'radius' (in km).
    % INPUT:
    %       - stations ... structure with 'lon' and 'lat' fields (in degree)
    %       - radius   ... distannce by which polygon around stations is
    %                      expanded
    % OUTPUT:
    %       - poly ... structure with 'lon' and 'lat' fields (in degree) of
    %         corner points of polygon encompassing all stations
    % =====================================================================
    if  nargin<nargin(@polygon_from_stations)
        error('Not enough input arguments.');
    end
    
    % remove duplicate locations
    if isrow(stations.lon), stations.lon = stations.lon'; end
    if isrow(stations.lat), stations.lat = stations.lat'; end
    ll = unique([stations.lon, stations.lat], 'rows', 'stable');
    
    % get polygon around all observational locations
    ixy = convhull(ll(:,1), ll(:,2), 'Simplify', true);
    % calculate centre point of polygon
    warning off
    [pLonC, pLatC] = centroid(polyshape(ll(ixy,1), ll(ixy,2)));
    warning on
    % widen polygon by calculating vectors describing outer points
    % from centre point and increasing their length by the defined radius
    v = [ll(ixy,1) - pLonC, ll(ixy,2) - pLatC];
    d = distance_on_earth(deg2rad(ll(ixy,:)), deg2rad([pLonC, pLatC]));
    fac = 1. + 1.e3*radius./d;
    v_new = v .* fac;
    poly.lon = pLonC + v_new(:,1)';
    poly.lat = pLatC + v_new(:,2)';
end
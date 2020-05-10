function distance = distance_on_earth(loc1, loc2)
% distance_on_earth calculates the distance between two (or more)
% geographical locations, loc1 and loc2
%
% input arguments:
% - loc1: M-by-2 array containing M decimal [lon,lat]-locations [rad]
% - loc2: N-by-2 array containing N decimal [lon,lat]-locations [rad]
%
% return value:
% - distance [m]: if M==N     => vector of length M containing the distance
%                                for each i-th pair of [loc1(i,:),loc2(i,:)]
%                 if M>1,N==1 => vector of length M containing the distance
%                                for each i-th pair of [loc1(i,:),loc2]
%                 if M==1,N>1 => vector of length N containing the distance
%                                for each i-th pair of [loc1,loc2(i,:)]
%                 else        => ERROR

if nargin<2, error('Function requires two input arguments'); end

if  (size(loc1,2)~=2) || (size(loc2,2)~=2)
    error('Location arrays must have exactly to columns: [lon, lat].');
end

% check if values are in [rad]
if  (sum(min([loc1(:,1);loc2(:,1)])<-2*pi)>0) || ...
    (sum(max([loc1(:,1);loc2(:,1)])> 2*pi)>0) || ...
    (sum(min([loc1(:,2);loc2(:,2)])<-2*pi)>0) || ...
    (sum(max([loc1(:,2);loc2(:,2)])> 2*pi)>0)    
    error('Values must be in radian.');
end

% check size of location arrays
if  sum(size(loc1))~=sum(size(loc2))
    if  (sum(size(loc1))~=3) && (sum(size(loc2))~=3)
        error('Location arrays must have same size, otherwise at least one array must be of size [1,2]');
    end
    if  sum(size(loc1))==3
        loc1 = repmat(loc1,size(loc2,1),1);
    else
        loc2 = repmat(loc2,size(loc1,1),1);
    end
end

earthRadius = 6371000; % [m]
dLon = abs(loc2(:,1) - loc1(:,1)); % difference in longitude
dLat = abs(loc2(:,2) - loc1(:,2)); % difference in latitude
f_haversine = ...         % haversine function
              sin(0.5*dLat).*sin(0.5*dLat) + ...
              cos(loc1(:,2)).*cos(loc2(:,2)) .* sin(0.5*dLon).*sin(0.5*dLon);

distance = 2*earthRadius*atan2(sqrt(f_haversine),sqrt(1-f_haversine));

return

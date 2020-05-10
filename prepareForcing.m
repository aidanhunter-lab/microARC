function [Forc, FixedParams] = prepareForcing(F,fixedParams)
% Extract forcing data used in size-structured model from F, and interpolate over
% the depth layers specified in FixedParams.

FixedParams = fixedParams;

for iy = 1:length(FixedParams.years)    
    y_index = ['y' num2str(FixedParams.years(iy))];
    forcing = F.(y_index);    
    % exclude trajectories that go on land
    iSurf = forcing.z == max(forcing.z);
    iKeep = sum(~isnan(forcing.T(:,iSurf,:)),1) == size(forcing.T,1);
    forcing = remove_landTrajs(forcing,iKeep);    
    % forcing data dimensions for each year
    FixedParams.(y_index).nt = length(forcing.t);    % number of forcing data time steps
    FixedParams.(y_index).nTraj = length(forcing.iTraj); % number of forcing trajectories            
    FixedParams.(y_index).lat = forcing.y; % latitude and longitude
    FixedParams.(y_index).lon = forcing.x;
    
    % interpolate forcing data over depth to match dimensions of biological model
    % dry-wet mask    
    H3d = repmat(reshape(forcing.H, [size(forcing.H,1) 1 size(forcing.H,2)]), ...
        [1 FixedParams.nz 1]);
    z3d = repmat(reshape(FixedParams.z, [1 FixedParams.nz 1]), ...
        [length(forcing.t) 1 length(forcing.iTraj)]);    
    wet = z3d > H3d;
    % if necessary extend 3D tracer fields to surface
    if  max(FixedParams.z)>max(forcing.z)
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.T(:,1,:), forcing.T);
    else
        z_ext = forcing.z;
        v_ext = forcing.T;
    end
    
    % temperature
    forcing.T = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % temperature at centre of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);    
    forcing.T = flip(gapFill_forc(flip(forcing.T,2), wet),2); % fill gaps
    
    % diffusivity - find at depth layer edges and centers
    if  max(FixedParams.z)>max(forcing.z)
        v_ext = cat(2, forcing.kv(:,1,:), forcing.kv);
    else
        v_ext = forcing.kv;
    end
    forcing.kv_center = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % vertical diffusivities at edges of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);        
    forcing.kv_center = flip(gapFill_forc(flip(forcing.kv_center,2), wet),2); % fill gaps        
    forcing.kv_center = forcing.kv_center * 24*60*60; % convert m2/s -> m2/day    
    forcing.kv_center(~wet) = 0.0; % set diffusivity on land/in sediment to zero to avoid mixing into sediment
    
    forcing.kv = interp_forc(forcing.kv, 1:length(forcing.iTraj) , ... % vertical diffusivities at edges of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, forcing.z, FixedParams.zw);    
    forcing.kv = forcing.kv(:,2:end-1,:); % only require diffusivities between depth layers, not at surface or bottom
    wet = wet(:,2:end,:);    
    forcing.kv = flip(gapFill_forc(flip(forcing.kv,2), wet),2); % fill gaps    
    forcing.kv = forcing.kv * 24*60*60; % convert m2/s -> m2/day    
%     forcing.kv = max(min(forcing.kv, 10^-2), 10^-5.5) * 86400; % convert m2/s -> m2/day
    forcing.kv(~wet) = 0.0; % set diffusivity on land/in sediment to zero to avoid mixing into sediment
    
    % calculate PAR at center of depth layers
    PARfrac = 0.43; % photosynthetically available fraction of incoming shortwave radiation at surface - PARAMETER COPIED FROM BIOMAS MODEL    
    PARsurf = repmat(reshape(PARfrac .* forcing.swrad, ... 
        [FixedParams.(y_index).nt 1 FixedParams.(y_index).nTraj]), ...
        [1 FixedParams.nz 1]);
    attSW = 0.04; % light attenuation by seawater (1 / m)
    att_z = attSW; % light is also attenuated by plankton, so this calcualtion may be better within the main model function - as BIOMAS does...        
    forcing.PAR = PARsurf .* repmat(exp(-att_z .* abs(FixedParams.z))', ... 
        [FixedParams.(y_index).nt 1 FixedParams.(y_index).nTraj]);    
    % Convert PAR from W/m2 -> muEin/s/m2    
    h=6.62607004e-34; % Planck constant
    c=2.99792458e8;   % light speed constant (should this be changed to speed in seawater?)
    a=6.02214086e23;  % Avogadro constant
    l=570e-9;         % assume average wavelength of 570 nm    
    e=h*c/l;          % energy of single photon (J)
    E=a*e;            % energy (J / mole) => 1 Einstein = E J => 1 J = 1/E Einstein => 1 J = 1e6/E micro Einsteins
    
    forcing.PAR = 1e6/E .* forcing.PAR;
    

    % Variables useful to set initial conditions, NO3, PS and PL
    wet = z3d > H3d;

    % NO3
    if  max(FixedParams.z)>max(forcing.z)
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.NO3ic(:,1,:), forcing.NO3ic);
    else
        z_ext = forcing.z;
        v_ext = forcing.NO3ic;
    end
    forcing.NO3ic = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % NO3 at centre of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);    
    forcing.NO3ic = flip(gapFill_forc(flip(forcing.NO3ic,2), wet),2); % fill gaps
    
    % PS
    if  max(FixedParams.z)>max(forcing.z)
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.PSic(:,1,:), forcing.PSic);
    else
        z_ext = forcing.z;
        v_ext = forcing.PSic;
    end
    forcing.PSic = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % PS at centre of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);    
    forcing.PSic = flip(gapFill_forc(flip(forcing.PSic,2), wet),2); % fill gaps
    
    % PL
    if  max(FixedParams.z)>max(forcing.z)
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.PLic(:,1,:), forcing.PLic);
    else
        z_ext = forcing.z;
        v_ext = forcing.PLic;
    end
    forcing.PLic = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % PL at centre of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);    
    forcing.PLic = flip(gapFill_forc(flip(forcing.PLic,2), wet),2); % fill gaps
    
    Forc.(y_index) = forcing;
end

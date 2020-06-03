function Forc = prepareForcing(F,fixedParams)

% Extracts required forcing data from F, omitting spurious trajectories,
% and interpolates over depth layers specified in FixedParams.

FixedParams = fixedParams;

for iy = 1:length(FixedParams.years)    
    y_index = ['y' num2str(FixedParams.years(iy))];
    forcing = F.(y_index);
    
    % exclude trajectories that go on land
    iSurf = forcing.z == max(forcing.z);    
    iKeep = sum(~isnan(forcing.T(:,iSurf,:)),1) == size(forcing.T,1);
    forcing = remove_landTrajs(forcing,iKeep);    
    
    nt = length(forcing.t);
    np = length(forcing.iTraj);
    forcing.iTraj = 1:length(forcing.iTraj);
    
    
    % interpolate forcing data over depth to match dimensions of biological model
    H3d = repmat(reshape(forcing.H, [size(forcing.H,1) 1 size(forcing.H,2)]), ...
        [1 FixedParams.nz+1 1]);
    z3d = repmat(reshape(FixedParams.zw, [1 FixedParams.nz+1 1]), ...
        [length(forcing.t) 1 length(forcing.iTraj)]);    
    % dry-wet mask    
    wet = z3d > H3d;
    wet = wet(:,2:end,:);
    forcing.wet = wet;
    
    % temperature
    if  max(FixedParams.z)>max(forcing.z) % if necessary extend 3D tracer fields to surface
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.T(:,1,:), forcing.T);
    else
        z_ext = forcing.z;
        v_ext = forcing.T;
    end    
    forcing.T = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % temperature at centre of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);    
    forcing.T = flip(gapFill_forc(flip(forcing.T,2), wet),2); % fill gaps

    % diffusivity - at depth layer centers
    if  max(FixedParams.z)>max(forcing.z)
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.kv(:,1,:), forcing.kv);
    else
        z_ext = forcing.z;
        v_ext = forcing.kv;
    end    
    forcing.K_center = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % vertical diffusivities at centers of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.z);    
    forcing.K_center = flip(gapFill_forc(flip(forcing.K_center,2), wet),2); % fill gaps
    forcing.K_center = forcing.K_center * (24*60*60); % convert m2/s -> m2/day
    forcing.K_center(~wet) = 0.0; % set diffusivity on land/in sediment to zero to avoid mixing into sediment

    % diffusivity - at depth layer edges
    if  max(FixedParams.zw(2:end))>max(forcing.z)
        z_ext = [0; forcing.z];
        v_ext = cat(2, forcing.kv(:,1,:), forcing.kv);
    else
        z_ext = forcing.z;
        v_ext = forcing.kv;
    end
    
    forcing.K = interp_forc(v_ext, 1:length(forcing.iTraj) , ... % vertical diffusivities at edges of depth layers
        1:length(forcing.iTraj), forcing.t, forcing.t, z_ext, FixedParams.zw);    
    forcing.K = forcing.K(:,2:end-1,:); % only require diffusivities between depth layers, not at surface or bottom
    wet_w = wet(:,2:end,:);    
    forcing.K = flip(gapFill_forc(flip(forcing.K,2), wet_w),2); % fill gaps
    forcing.K = forcing.K * (24*60*60); % convert m2/s -> m2/day
    forcing.K(~wet_w) = 0.0; % set diffusivity on land/in sediment to zero to avoid mixing into sediment
    
    forcing = rmfield(forcing, 'kv');
    
    % Convert surface light from swrad (W/m^2) to PAR (muEin/s/m^2)
    PARfrac = 0.43; % photosynthetically available fraction of incoming shortwave radiation at surface
    forcing.PARsurf = reshape(PARfrac .* forcing.swrad, [nt 1 np]); % convert swrad to PAR (W/m^2)    
    h=6.62607004e-34; % Planck constant
    c=2.99792458e8;   % light speed constant (should this be changed to speed in seawater?)
    a=6.02214086e23;  % Avogadro constant
    l=570e-9;         % assume average wavelength of 570 nm    
    e=h*c/l;          % energy of single photon (J)
    E=a*e;            % energy (J / mole) => 1 Einstein = E J => 1 J = 1/E Einstein => 1 J = 1e6/E micro Einsteins    
    forcing.PARsurf = 1e6/E .* forcing.PARsurf; % surface PAR (muEin/s/m^2)
    % Light attenuation is calculated within ODEs because it depends on
    % plankton concentrations 
    

    % Variables useful to set initial conditions, NO3, PS and PL

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
    
    if iy == 1, Forc.years = FixedParams.years; end
    
    forcing = rmfield(forcing, 'z');
    forcing.t = repmat(forcing.t, [1 np]);
 
    fields = fieldnames(forcing);
    for j = 1:length(fields)
        fs = size(forcing.(fields{j}));
        if length(fs) == 3 && (fs(1) == nt && fs(3) == np)
            forcing.(fields{j}) = permute(forcing.(fields{j}), [2 1 3]);
        end
        fl = length(fs);
        if iy == 1
            Forc.(fields{j}) = [];
        end
        Forc.(fields{j}) = cat(fl, Forc.(fields{j}), forcing.(fields{j}));
    end
        
end

nTraj = length(Forc.iTraj);

% Account for changes in modelled water column depth caused by
% shallow portions of particle trajectories.
wet = Forc.wet;
zChange = [zeros(1, 1, nTraj) diff(sum(wet))]; % trajectories & time steps where water column deepens/shallows
% If water column deepens then fill new depth layers assuming
% concentrations equal those of previous deepest layer
deepens = zChange > 0;
infillDepth = [false(FixedParams.nz, 1, nTraj) ... 
    (wet(:,2:end,:) & ~wet(:,1:end-1,:))]; % depths to infill from layer above
replicateConc = [infillDepth(2:end,:,:) & ~infillDepth(1:end-1,:,:); ... 
    false(1, FixedParams.nt, nTraj)]; % index concentration to replicate into deeper layers

Forc.deepens = deepens;
Forc.infillDepth = infillDepth;
Forc.replicateConc = replicateConc;


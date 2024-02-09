% plot modelled size spectra either as histograms or as spectra

function plt = plotModelSizeSpectra(Out,auxVars, fixedParams, data, forc, type, ...
    var, plankton, time, intAve, varargin)

% out: model output
% auxVars contains cell abundance and biovolume
% type: 'histogram' or 'spectrum'
% var: one of: 'abundance', 'biovolume', 'Cbiomass', 'Nbiomass', 'Chl'
% plankton: 'phyto' or 'zoo'
% time: time (doy of year)
% intAve: 'integrated' or 'averaged' over depth 

% optional

% trajectory selection
% per default all trajs are selected
% traj: select trajectories by their ID (one or more)
% event: mean of all trajectories associated with this event ID
% waterOrigin: all 'Arctic'or all 'Atlantic' trajectories

% depth selection
% per default all depths are used
% depthLayer: one depth layer ([m])
% maxDepth: all depth layers until maximum depth ([m])

% if type is 'spectrum':
% normalised: only if type is spectrum, should spectrum be normalised? true
%     or false


% for testing
% fixedParams = FixedParams; 
% data = Data; 
% forc = Forc; 
% type = 'histogram';
% var = 'Cbiomass';
% var = 'abundance';
% plankton = 'phyto';
% time = 200; 
% intAve = 'averaged'; 



% get size intervals
ESD_edges = fixedParams.PPdia_intervals;


% get size binned data of variable var
if strcmp(var, 'Cbiomass')
    if strcmp(plankton, 'phyto')
        mod = squeeze(Out.P(:,:,fixedParams.PP_C_index,:,:));
    elseif strcmp(plankton, 'zoo')
        mod = squeeze(Out.P(:,:,fixedParams.ZP_C_index,:,:));
    end
    unit = 'mmol C'; 
elseif strcmp(var, 'Nbiomass')
    if strcmp(plankton, 'phyto')
        mod = squeeze(Out.P(:,:,fixedParams.PP_N_index,:,:)); 
    elseif strcmp(plankton, 'zoo')
        mod = squeeze(Out.P(:,:,fixedParams.ZP_N_index,:,:)); 
    end
    unit = 'mmol N'; 
elseif strcmp(var, 'Chl')
    if strcmp(plankton, 'phyto')
        mod = squeeze(Out.P(:,:,fixedParams.PP_Chl_index,:,:)); 
    elseif strcmp(plankton, 'zoo')
        mod = squeeze(Out.P(:,:,fixedParams.ZP_Chl_index,:,:)); 
    end
    unit = 'mg Chl';
elseif strcmp(var, 'abundance') 
    if strcmp(plankton, 'phyto')
        mod = auxVars.cellDensity(fixedParams.phytoplankton,:,:,:);
    elseif strcmp(plankton, 'zoo')
        mod = auxVars.cellDensity(fixedParams.zooplankton,:,:,:);
    end
    unit = 'cells';
    
elseif strcmp(var, 'biovolume')
    if strcmp(plankton, 'phyto')
        mod = auxVars.biovolume(fixedParams.phytoplankton,:,:,:);
    elseif strcmp(plankton, 'zoo')
        mod = auxVars.biovolume(fixedParams.zooplankton,:,:,:);
    end
    unit = 'µm^3';
end
% mod(sizeclass, depth, doy, itraj)

% pick time
mod = squeeze(mod(:, :, time, :));
timeLabel = ['doy = ' num2str(time)];

% pick trajectories
% by traj (one or more), event, waterOrigin
% else use average of all trajectories (mod stays as is)

% per default set trajectory label to all trajectories
trajLabel = 'all trajectories';

if ~isempty(varargin)
    if any(strcmp(varargin, 'traj'))
        traj = varargin{find(strcmp(varargin, 'traj'))+1};
        if ~ismember(forc.iTraj, traj)
            warning('Select only trajectories present in Forc struct')
        end
        % get indices of trajs
        trajIndex = find(ismember(forc.iTraj, traj));
        
        % filter mod
        mod = mod(:, :, trajIndex); % tst = 
        % set trajectory label
        ntraj = length(traj); 
        if ntraj > 1
            trajLabel = ['trajectories ', regexprep(num2str(traj),'\s+',',')];
        else
            trajLabel = ['trajectrory ' num2str(traj)];
        end
        
    elseif any(strcmp(varargin, 'event'))
        event = varargin{find(strcmp(varargin, 'event'))+1};
        if data.scalar.nEvents < event
            warning('Select only events present in Data struct')
        end
        % get indices of all trajs associated to this event(s)
        trajIndex = data.scalar.evTraj(:, event);
        trajIndex = unique(trajIndex(:));
        
        % filter mod
        mod = mod(:, :, trajIndex);  
        % set trajectory label
        trajLabel = ['all trajectories associated with event ' num2str(event)];
        
        
    elseif any(strcmp(varargin, 'waterOrigin'))
        waterOrigin = varargin{find(strcmp(varargin, 'waterOrigin'))+1};
        if all([~strcmp('Atlantic', waterOrigin) ~strcmp('Arctic', waterOrigin)])
            warning('waterOrigin must be either "Arctic" or "Atlantic"')
        end
        % get indices of all trajs associated with this watermass
        trajIndex = strcmp(forc.waterMass, waterOrigin); 
        
        
        % filter mod
        mod = mod(:, :, trajIndex);  
        % set trajectory label
        trajLabel = ['all ' waterOrigin ' trajectories'];
    end
end

% get mean over all trajectories, result should be 9x9 matrix
mod = mean(mod, 3);


% deal with depth: pick one, integrate or average to max depth
% default: no depth is selected, all are used.
depthLabel = [intAve ' over the whole water column'];
if ~isempty(varargin)
    if any(strcmp(varargin, 'depthLayer'))
        depthLayer = varargin{find(strcmp(varargin, 'depthLayer'))+1};
        
        % depth layer is in [m], so find closest layer in model
        [minValue,closestIndex] = min(abs(fixedParams.z*-1-depthLayer));
        closestValue = fixedParams.z(closestIndex);
        
        % filter for size class
        mod = mod(:,closestIndex);
        % set depthLabel
        depthLabel = ['at ' num2str(depthLayer) 'm (' num2str(closestValue) 'm) depth'];
    elseif any(strcmp(varargin, 'maxDepth'))
        maxDepth = varargin{find(strcmp(varargin, 'maxDepth'))+1};
       
        
        % depth layer is in [m], so find closest layer in model
        [minValue,closestIndex] = min(abs(fixedParams.z*-1-maxDepth));
        closestValue = fixedParams.z(closestIndex);
        
        % keep only layers until maxDepth
        mod = mod(:, (1:closestIndex)); 
        % set depthLabel
        depthLabel = [intAve ' over 0m-' num2str(maxDepth) 'm (' num2str(closestValue) 'm) depth']; 
        
        
    end 
end

% now average or integrate over depths (if theres more than 1 depth)
numDepth = size(mod, 2);
if numDepth > 1
    if strcmp(intAve, 'integrated')
        % multiply mod with thickness of each depth class
        Y = mod .* reshape(fixedParams.zwidth(1:numDepth), [1 numDepth]);
        % get sum of all depths per size class
        mod = sum(Y, 2);

    elseif strcmp(intAve, 'averaged')

        % multiply mod with thickness of each depth class
        Y = mod .* reshape(fixedParams.zwidth(1:numDepth), [1 numDepth]);
        % get sum of all depths per size class
        Y2 = sum(Y, 2);
        % and divide over thickness of selected water column
        mod = Y2 ./ sum(fixedParams.zwidth(1:numDepth));
    end
end

% mod should now be 9x1 struct (one value per size class)
% begin plotting

plt = figure; 

if strcmp(type, 'histogram')
    
  histogram('BinEdges',ESD_edges','BinCounts', mod)
  ax = gca;
  ax.XScale = 'log'; 
  ax.YScale = 'log';
  xlabel('cell diameter (µm)')
  xticks(round(ESD_edges, 2))
  
  switch intAve
      case 'integrated'
          yunit = [unit '/ m^2'];
      case 'averaged'
          yunit = [unit '/ m^3'];
  end
  ylabel([plankton 'plankton ' var ' ( ' yunit ' )'])
  ylim([min(mod)*0.9 max(mod)*1.1])
  title({['Size binned ' plankton 'plankton ' var ', ' timeLabel]  , ...
    [trajLabel ' ' depthLabel]})

    
elseif strcmp(type, 'spectrum')
    % scale binned data to size vector, (calculate densities)
    % then plot
    
    % get cell concentration/biomass/biovolume densities from mod
    dlog10ESD = diff(log10(ESD_edges(1:2))); 
   
    modSpec = mod ./ dlog10ESD; 
    modNorm = modSpec ./ sum(mod); % normalise
    
    % default: plot spectrum (non-normalised)
    modPlot = modSpec; 
    % set unit for modSpec
    switch intAve
      case 'integrated'
          yunit = [unit ' / m^2 / log_{10}(ESD/1 µm)'];
      case 'averaged'
          yunit = [unit ' / m^3 / log_{10}(ESD/1 µm)'];
    end
    ylab = [ plankton 'plankton ' var ' density ( ' yunit ' )' ]; 
    
    
    % test: is sum of mod same as intrgral of modSpec?
    % sum(modSpec) * dlog10ESD
    % sum(mod)
    % sum(modNorm) * dlog10ESD 
    

    
    % in case of normalised spectrum plot:  
    if ~isempty(varargin)
        if any(strcmp(varargin, 'normalised'))
            normalised = varargin{find(strcmp(varargin, 'normalised'))+1};
            if normalised
                modPlot = modNorm
                ylab = [ 'normalised ' plankton 'plankton ' var ' density ( log_{10}(ESD/1 µm)^{-1} )' ];
            end
        end
    end
    
    
    % repeat each element of modPlot twice
    modPlot = repelem(modPlot, 2);
    % and get according size vector
    ESDvec = [ESD_edges(1), repelem(ESD_edges(2:end-1), 2)', ESD_edges(end)]'; 
    
    % plot
    loglog(ESDvec, modPlot)
    xlabel('cell diameter (µm)')
    xticks(round(ESD_edges, 2))
    ylabel(ylab)
    title({[plankton 'plankton ' var ' density spectrum, ' timeLabel]  , ...
    [trajLabel ' ' depthLabel]})

else 
    warning('Specify plot type (histogram or spectrum)')
end


end
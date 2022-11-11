function modData = matchModOutput2Data2(out, auxVars, Data, FixedParams, varargin)
% Returns model output that 'matches' the data used for optimisation.
% Output struct 'modData' has the same form as the 'Data' struct.

% Model output times are matched to sample times; outputs are interpolated
% over depths to match sampled depth layers.
% Then model outputs are transformed using the same functions that
% standardised the data.

fitToFullSizeSpectra = []; % determines how modelled output is matched to size data
extractVarargin(varargin)
%evTrajID can be imported manually
if isempty(fitToFullSizeSpectra) % || ~islogical(fitToFullSizeSpectra)
    fitToFullSizeSpectra = false; % By default assume that the binned/integrated size data is used to tune parameters
end


ntrajOut = size(out.N, ndims(out.N));
ntrajData = size(Data.scalar.EventTraj, 2);

if ntrajOut ~= ntrajData
    warning('Number of trajectories in model output "out" does not equal number of trajectories used in from the Data.')
    if ~exist('evTrajIX') 
        error('evTrajIX needs to be imported')  % this is summerTrajsIndex or autumnTrajsIndex, made when Forc is further subsetted for single trajecory runs.
    end
end



%% Scalar data
Vars = Data.scalar.obsInCostFunction;

evTraj = Data.scalar.evTraj; % trajectories used for each sample event
nsamples = size(evTraj, 1); % number of trajectories per sampling event

nEvent = Data.scalar.nEvents; % number of sample events
depths_mod = abs(FixedParams.z); % modelled depth layers

% Build modData to match the form of Data -- we don't need to include every
% field, only the fields required for the cost function.
modData.scalar.Yearday = nan(Data.scalar.nSamples,1);
modData.scalar.Depth = nan(Data.scalar.nSamples,1);
modData.scalar.waterMass = cell(Data.scalar.nSamples,1);
modData.scalar.Variable = cell(Data.scalar.nSamples,1);
modData.scalar.Value = nan(Data.scalar.nSamples, nsamples);
modData.scalar.scaled_Value = modData.scalar.Value;

% Standardise model output with respect to depth and event using linear mixed models
for i = 1:nEvent
    iEvent = Data.scalar.Event == i; % index event
    itraj = evTraj(:,i); % trajectories used for event i (their indeces in Forc.iTraj!; not their IDs)
    
    % skip the rest when no modelled traj belongs to this event 
    % but only, when ntrajOut ~= ntrajData
     if ntrajOut ~= ntrajData
        % look up if any modelled trajectory is a member of of itraj
        % evTrajIX: (original; before subsetting Forc) indices of modelled trajectories
        % evTraj: List of (original) trajectory indices per event

        if ~any(ismember(itraj, evTrajIX)) % in case none of the moddeled trajs matches an observation event
            continue % skip to next iteration without running the rest of the for loop
        end
        
        % reassign which trajs should be extracted from model output 
        itraj = find(ismember(evTrajIX, itraj));
    end
    
        vars = unique(Data.scalar.Variable(iEvent)); % variables measured in event i
        vars = vars(ismember(vars, Vars));
        iEvent = ismember(Data.scalar.Variable, vars) & iEvent; % omit unmodelled variables from the event index
        Yearday = Data.scalar.Yearday(find(iEvent, 1));
        Depth = Data.scalar.Depth(iEvent);
        clear waterMass
        [waterMass{1:sum(iEvent),1}] = deal(Data.scalar.waterMass{i});
        Variable = Data.scalar.Variable(iEvent);
        modData.scalar.Yearday(iEvent) = Yearday;
        modData.scalar.Depth(iEvent) = Depth;
        modData.scalar.waterMass(iEvent) = waterMass;
        modData.scalar.Variable(iEvent) = Variable;
        for j = 1:length(vars)
            % loop through all data types sampled during event i
            jvar = vars{j};
            ind = iEvent & strcmp(Data.scalar.Variable, jvar);
            depth = modData.scalar.Depth(ind);
            switch jvar
                % modelled values [depth,trajectory]
                case 'N'
                    ymod = squeeze(out.N(:,:,Yearday,itraj));
                case 'PON'
                    ymod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_N_index,Yearday,itraj));
                case 'POC'
                    ymod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_C_index,Yearday,itraj));
                case 'chl_a'
                    ymod = squeeze(sum(out.P(:,:,FixedParams.PP_Chl_index,Yearday,itraj)));
            end
            % Interpolate modelled output to match observation depths. The
            % 'extrap' option prevents NaNs at observation depths below the
            % midpoints of the deepest modelled layer. 'extrap' is not required
            % if samples from deeper than the midpoint of the deepest modelled
            % layer are discarded; it is required when all data within the
            % deepest modelled layer are retained -- see omitDeepSamples.m.
            ymod = interp1(depths_mod, ymod, depth, 'linear', 'extrap');
            ymod_scaled = Data.scalar.(['scaleFun_' jvar])(Data.scalar.scale_mu(ind), ...
                Data.scalar.scale_sig(ind), ymod); % scale model output using same functions that scaled the data
%             modData.scalar.Value(ind,:) = ymod;  
%             modData.scalar.scaled_Value(ind,:) = ymod_scaled; % change
%             back in case assignment does not work any more for full model
%             runs
            modData.scalar.Value(ind,1:size(ymod,2)) = ymod;
            modData.scalar.scaled_Value(ind,1:size(ymod_scaled,2)) = ymod_scaled;
        end
end

% Input names of measurements not included in model/cost function
modData.scalar.Variable(~Data.scalar.inCostFunction) = ...
    Data.scalar.Variable(~Data.scalar.inCostFunction);


%% Size spectra

switch fitToFullSizeSpectra
    
    case 'trueByWaterMass'
        
        % creates waterOriginSpectra (similar to S-Spectra in the Lampe et al 2021 paper) by aggregating
        % by trajectory origin (full size spectra)
        
        % Forc must be passed as varArgIn, to look up watermass of
        % simulated trajectories!
        if ~exist('Forc','var')
            error('For this method "Forc" must be passed as an optional argument' )
        end
        
       
        timespan = [Data.size.YeardayFirst(1) : Data.size.YeardayLast]; % get time range of sampling campaigns
        % get depth range of observtions
        minDepth = Data.size.DepthMin(1);
        maxDepth = Data.size.DepthMax(1);
        
        % FixedParams.zw % edges of depth layers
        [minVal, minIndex] = min(abs(-minDepth-FixedParams.zw));
        [minVal, maxIndex] = min(abs(-maxDepth-FixedParams.zw));
        maxIndex = maxIndex-1;  % only layers shallower than closest depth layer edge! 
          clear minVal
        depthIndexRange = [minIndex:maxIndex];
        
        % get indices of trajectories by watermass
        atlTrajs = strcmp(Forc.waterMass, 'Atlantic' );
        arcTrajs = strcmp(Forc.waterMass, 'Arctic' );
        
        ESD = unique(Data.sizeFull.ESD);
        nESD = length(ESD);
        
        % average all trajs by water mass 
        % whithin the time frame of the sampling campaigns of S-spectra
        % keep only upper depth layers x (50?)m 
        % result: to average spectra, one for arctic and one for atlantic
        % origin -> these are the eqiuvalents to the S spectra
        
       
        allVarsSize = {'cellDensity', 'BioVolDensity'}';
        
        % prepare modData struct, compatible to Data struct 
        modData.size.scenario = Data.size.scenario;
        
            regime = Data.size.regime; 
            wm = strrep(regime, 'warm', 'Atlantic'); 
            wm = strrep(wm, 'cold', 'Arctic');
        
        modData.size.waterMass = wm;
        modData.size.trophicLevel = Data.size.trophicLevel;
        modData.size.ESD = Data.size.ESD;
        modData.size.cellDensity = nan(size(Data.size.cellDensity));
        modData.size.BioVolDensity = nan(size(Data.size.BioVolDensity));
        
        
        watermasses = unique(modData.size.waterMass);
        trophicLevels = unique(modData.size.trophicLevel);
        for i =  1:length(watermasses)
            WM = watermasses{i};
            indWM = strcmp(modData.size.waterMass, WM);
            
            switch WM
                case 'Atlantic'
                    WMtrajs = atlTrajs; 
                case 'Arctic'
                    WMtrajs = arcTrajs;
            end
            
            
            for j = 1:length(trophicLevels)
                TL = trophicLevels{j};
                indTL = strcmp(modData.size.trophicLevel, TL);
                
                switch TL
                    case 'autotroph'
                        planktonclass = FixedParams.phytoplankton;
                    case 'heterotroph' 
                        planktonclass = FixedParams.zooplankton;
                end
                
                indY = indWM & indTL; % assign extrapolated model output here
                
                % find and average model output correspondenting to indY
                for k = 1:length(allVarsSize)
                    var = allVarsSize{k}; 
                    
                    switch var
                        case 'cellDensity'
                            ymod = auxVars.cellDensity;
                        case 'BioVolDensity'
                            ymod = auxVars.biovolume;
                    end
                    % filter for class, depth, time, watermass
                    ymod = ymod(planktonclass,depthIndexRange,timespan,WMtrajs);
                  
                    % average over size class (1st) dimension
                    ymod = mean(ymod, 2:4); % [cells or m^3 / m^3]
                    
                    % Transform and interpolate model output to match data.
                    % Units of ymod are either cells / m3 or m3 / m3 (within each size
                    % class interval). The full size spectra data have units
                    % cells / m3 / log10(ESD)
                    ymod_ = []; % [cells or m^3 / m^3 / log10(ESD)]
                    for sc = 1:FixedParams.nPP_size % only works if PP and ZP size sclasses are equal
                        y = ymod(sc); % [cells or m^3 / m^3]
                        interval = FixedParams.PPdia_intervals(sc:sc+1); % size class j interval boundaries [mu m]
                        ESDsc = ESD(interval(1) <= ESD & ESD <= interval(2));
                        nESDsc = length(ESDsc);
                        y_ = repmat(y ./ diff(log10(interval)), [nESDsc, 1]); % [cells or m^3 / m^3 / log10(ESD)]
                        ymod_ = vertcat(ymod_, y_); % [cells or m^3 / m^3 / log10(ESD)]  
                    end
                    % size(ymod_)
                    
                    % ymod_ now contains cell concentration or biovolume
                    % densities! 
                    % assign ymod_ to modData 
                    modData.size.(var)(indY) = ymod_; % [cells or m^3 / m^3 / log10(ESD)]

                end 
            end
        end

        
    
    case false

        allVarsSize = unique(Data.size.dataBinned.Variable);
        waterMasses = []; % water origin may not be specified if using data aggregated over all sampling events
        if isfield(Data.size.dataBinned, 'waterMass')
            waterMasses = unique(Data.size.dataBinned.waterMass);
        end
        default = isempty(waterMasses); % use fully aggregated size data and all trajectories by default
        
        % There are differences between the fully aggregated size data and the
        % sample-specific size data => use if-statement to construct modData.size
        
        if ~default
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Construction of modData when fitting to sample-specific size data
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Size-data vectors were derived from sample-specific data that was
            % grouped by water origin. There are data vectors associated to each
            % water origin group used to fit model. Each data vector is therefore
            % associated to specific particle trajectories => extract these as
            % trajectories from model output, then match times and depths to data.
            nwaterMasses = length(waterMasses);
            % Index trajectories, times and depths to extract from model output --
            % these must match the data
            if default
                traj = true(1, size(Data.scalar.EventTraj, 2));
            else
                traj = Data.size.evTraj;
                nsamples = size(traj, 1);
            end
            
            time = Data.size.dataBinned.Yearday(1); % the binned size data is averaged over sample events (time and depth) => all values are identical and we can use the 1st elements as indices
            depth = Data.size.dataBinned.Depth(1);
            
            % The size component of 'modData' is built to match the 'binnedData'
            % struct contained in 'Data.size'.
            n = size(Data.size.dataBinned.Year, 1);
            modData.size.waterMass = cell(n,1);
            modData.size.trophicLevel = cell(n,1);
            modData.size.Variable = cell(n,1);
            modData.size.Value = nan(n, numel(traj));
            % modData.size.Value_Atlantic = nan(n,length(trajAtl));
            % modData.size.Value_Arctic = nan(n,length(trajArc));
            
            for i = 1:length(allVarsSize)
                vs = allVarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, vs); % index variable
                modData.size.Variable(ind0,:) = Data.size.dataBinned.Variable(ind0);
                modData.size.trophicLevel(ind0,:) = Data.size.dataBinned.trophicLevel(ind0);
                modData.size.waterMass(ind0,:) = Data.size.dataBinned.waterMass(ind0);
                
                switch vs
                    case 'CellConc'
                        ymod = auxVars.cellDensity(:,:,time,traj(:));
                    case 'BioVol'
                        ymod = auxVars.biovolume(:,:,time,traj(:));
                end
                
                ymod = permute(ymod, [2 1 3 4]); % move depth dimension to front [depth, size, ntraj]
                ymodi = interp1(FixedParams.z , ymod, -abs(depth)); % interpolate model output to match observed depth
                ymodi = squeeze(ymodi); % dimension: [size, trajectories]
                
                % autotrophs
                ind1 = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                if ~isempty(waterMasses)
                    for j = 1:nwaterMasses
                        ind2 = ind1 & strcmp(Data.size.dataBinned.waterMass, waterMasses{j});
                        wi = strcmp(waterMasses{j}, Data.size.waterMass);
                        ind_ = (find(wi, 1) - 1)  * nsamples + 1:find(wi, 1, 'last') * nsamples;
                        modData.size.Value(ind2,ind_) = ymodi(FixedParams.phytoplankton,ind_);
                    end
                end
                
                % heterotrophs
                ind1 = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                if ~isempty(waterMasses)
                    for j = 1:nwaterMasses
                        ind2 = ind1 & strcmp(Data.size.dataBinned.waterMass, waterMasses{j});
                        wi = strcmp(waterMasses{j}, Data.size.waterMass);
                        ind_ = (find(wi, 1) - 1)  * nsamples + 1:find(wi, 1, 'last') * nsamples;
                        if FixedParams.nZP_size > 1
                            modData.size.Value(ind2,ind_) = ymodi(FixedParams.zooplankton,ind_);
                        else
                            modData.size.Value(ind2,ind_) = repmat(ymodi(FixedParams.zooplankton,ind_), [sum(ind2), 1]);
                        end
                    end
                end
                
            end
        else
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Construction of modData when fitting to fully aggregated size data
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % There is a single size-data vector (for each trophic level) assumed
            % to correspond to all sample events (fully aggregated).
            % For each sample, extract corresponding (trajectory, time, depth) model
            % output. All of these outputs are fitted to the same data vector.
            
            
            ind = ismember(Data.scalar.Year, unique(Data.size.Year)); % omit sample events from years without size data
            ev = unique(Data.scalar.Event(ind)); % all sample events
            et = evTraj(:,ev); % sampling events and associated trajectories
            nevent = size(et, 2);
            etime = nan(nevent, 1); % sample times of each event
            for i = 1:nevent
                ue = unique(Data.scalar.Yearday(Data.scalar.Event == ev(i)));
                etime(i) = ue(1); % this line looks like a quick-fix... check if there's an issue with double-sampling times in prepareFittingData.m
            end
            depths_obs = [min(Data.size.DepthMin) max(Data.size.DepthMax)]; % sample depth range (the aggregated size data don't specify precise depths)
            depth_ind = -FixedParams.zw(2:end) > depths_obs(1,1) & ...
                -FixedParams.zw(1:end-1) < depths_obs(1,2); % modelled depth layers corresponding to samples
            
            % VarsSize = Data.size.obsInCostFunction;
            allVarsSize = unique(Data.size.dataBinned.Variable);
            
            % The size component of 'modData' is built to match the 'binnedData'
            % sub-struct contained in 'Data'.
            n = size(Data.size.dataBinned.Year, 1);
            modData.size.Variable = cell(n,1);
            modData.size.trophicLevel = cell(n,1);
            modData.size.Value = nan(n,nsamples);
            modData.size.Value_allEvents = nan(n,nevent,nsamples);
            modData.size.scaled_Value = nan(n,nsamples);
            modData.size.scaled_Value_allEvents = nan(n,nevent,nsamples);
            
            for i = 1:length(allVarsSize)
                vs = allVarsSize{i};
                ind0 = strcmp(Data.size.dataBinned.Variable, vs);
                modData.size.Variable(ind0,:) = Data.size.dataBinned.Variable(ind0);
                modData.size.trophicLevel(ind0,:) = Data.size.dataBinned.trophicLevel(ind0);
                
                for j = 1:nsamples
                    itraj = et(j,:); % trajectory associated with each sampling event
                    [~, J] = max(sum(out.P(:,depth_ind,FixedParams.PP_Chl_index,etime,itraj))); % use modelled values from chl-max depth layer to compare to data
                    J = squeeze(J);
                    switch vs
                        case 'CellConc'
                            %                 ymod = auxVars.cellDensity(1:1:FixedParams.nPP_size,depth_ind,etime,itraj);
                            ymod = auxVars.cellDensity(:,depth_ind,etime,itraj);
                            ymod_ = nan(size(ymod,1), nevent);
                            for k = 1:nevent
                                ymod_(:,k) = ymod(:,J(k,k),k,k);
                            end
                            ymod = ymod_;
                            ymodMean = mean(ymod, 2); % average size-spectra over sampling events
                        case 'BioVol'
                            %                 ymod = auxVars.biovolume(1:FixedParams.nPP_size,depth_ind,etime,itraj);
                            ymod = auxVars.biovolume(:,depth_ind,etime,itraj);
                            ymod_ = nan(size(ymod,1), nevent);
                            for k = 1:nevent
                                ymod_(:,k) = ymod(:,J(k,k),k,k);
                            end
                            ymod = ymod_;
                            ymodMean = mean(ymod, 2); % average size-spectra over sampling events
                    end
                    
                    % autotrophs
                    ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                    
                    modData.size.Value(ind,j) = ymodMean(FixedParams.phytoplankton);
                    modData.size.Value_allEvents(ind,:,j) = ymod(FixedParams.phytoplankton,:);
                    
                    modData.size.scaled_Value(ind,j) = Data.size.(['scaleFun_' vs])( ...
                        Data.size.dataBinned.scale_mu(ind), ...
                        Data.size.dataBinned.scale_sig(ind), ymodMean(FixedParams.phytoplankton));
                    modData.size.scaled_Value_allEvents(ind,:,j) = Data.size.(['scaleFun_' vs])( ...
                        Data.size.dataBinned.scale_mu(ind), ...
                        Data.size.dataBinned.scale_sig(ind), ymod(FixedParams.phytoplankton,:));
                    
                    % heterotrophs
                    ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                    
                    if sum(FixedParams.zooplankton) ~= 1
                        % Number of modelled heterotrph size classes matches number
                        % of partitions in size data (there may only be one size class)
                        modData.size.Value(ind,j) = ymodMean(FixedParams.zooplankton);
                        modData.size.Value_allEvents(ind,:,j) = ymod(FixedParams.zooplankton,:);
                    else
                        modData.size.Value(ind,j) = repmat(ymodMean(FixedParams.zooplankton), [sum(ind) 1]);
                        modData.size.Value_allEvents(ind,:,j) = repmat(ymod(FixedParams.zooplankton,:), [sum(ind) 1]);
                    end
                    
                    modData.size.scaled_Value(ind,j) = Data.size.(['scaleFun_' vs])( ...
                        Data.size.dataBinned.scale_mu(ind), ...
                        Data.size.dataBinned.scale_sig(ind), ymodMean(FixedParams.zooplankton));
                    modData.size.scaled_Value_allEvents(ind,:,j) = Data.size.(['scaleFun_' vs])( ...
                        Data.size.dataBinned.scale_mu(ind), ...
                        Data.size.dataBinned.scale_sig(ind), ymod(FixedParams.zooplankton,:));
                end
            end
        end
        
        
    case true
        % this method for full size spectra (not binned), matched by
        % station
        % exchange all .size with .sizeFull
        
        allVarsSize = {'cellDensity'; 'BioVolDensity'};
        waterMasses = []; % water origin may not be specified if using data aggregated over all sampling events
        if isfield(Data.sizeFull, 'waterMass')
            waterMasses = unique(Data.sizeFull.waterMass);
        end
        % use fully aggregated size data and all trajectories by default (not a
        % good idea as size spectra seem to differ depending on water origin and
        % this influences numerically estimated parameter values).
       
        
        default = isempty(waterMasses);
        
        if ~default
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Construction of modData when fitting to sample-specific size data
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % Size-data vectors were derived from sample-specific data that was
            % grouped by water origin. There are data vectors associated to each
            % water origin group used to fit model. Each data vector is therefore
            % associated to specific particle trajectories => extract these as
            % trajectories from model output, then match times and depths to data.
            n = Data.sizeFull.nSamples;
            nwaterMasses = length(waterMasses);
            trophicLevels = unique(Data.sizeFull.trophicLevel);
            ntrophicLevels = length(trophicLevels);
            % Index trajectories, times and depths to extract from model output --
            % these must match the data
            traj = Data.sizeFull.evTraj; % particle trajectories used for each sample event
            nsamples = size(traj, 1); % number of partucle trajectories used for each sample event
            ESD = unique(Data.sizeFull.ESD);
            nESD = length(ESD);
            
            sampleEvents = unique(Data.sizeFull.Event, 'stable');
            eventTime = unique([Data.sizeFull.Event, Data.sizeFull.Yearday], 'rows');
            
            modData.sizeFull.Event = Data.sizeFull.Event;
            modData.sizeFull.Yearday = Data.sizeFull.Yearday;
            modData.sizeFull.Depth = Data.sizeFull.Depth;
            modData.sizeFull.trophicLevel = Data.sizeFull.trophicLevel;
            modData.sizeFull.ESD = Data.sizeFull.ESD;
            modData.sizeFull.cellVolume = Data.sizeFull.cellVolume;
            modData.sizeFull.evTraj = Data.sizeFull.evTraj;
            modData.sizeFull.waterMass = Data.sizeFull.waterMass;
            

            
            for i = 1:length(allVarsSize)
                modData.sizeFull.(allVarsSize{i}) = nan(n, 10);  %% nan(n, 1);
            end

            
            for i = 1:length(allVarsSize)
                vs = allVarsSize{i};
                switch vs
                    case 'cellDensity'
                        ymod = auxVars.cellDensity;
                    case 'BioVolDensity'
                        ymod = auxVars.biovolume;
                end
                
                % Transform and interpolate model output to match data.
                % Units of ymod are either cells / m3 or m3 / m3 (within each size
                % class interval). The full size spectra data have units
                % cells / m3 / log10(ESD)
                ymodSize = size(ymod);
                
                % autotrophs
                ymodP = reshape(ymod(FixedParams.phytoplankton,:), ...
                    [FixedParams.nPP_size, ymodSize(2:end)]);
                ymod_ = nan([nESD, prod(ymodSize(2:end))]);
                rowIndex = 0;
                for j = 1:FixedParams.nPP_size
                    y = ymodP(j,:);
                    interval = FixedParams.PPdia_intervals(j:j+1); % size class j interval boundaries
                    ESDj = ESD(interval(1) <= ESD & ESD <= interval(2));
                    nESDj = length(ESDj);
                    rowIndex = max(rowIndex) + (1:nESDj);
                    y_ = repmat(y ./ diff(log10(interval)), [nESDj, 1]);
                    ymod_(rowIndex,:) = y_;
                end
                ymodP = reshape(ymod_, [nESD, ymodSize(2:end)]); % model output units & number of size classes should now match the data
                
                % heterotrophs
                ymodZ = reshape(ymod(FixedParams.zooplankton,:), ...
                    [FixedParams.nZP_size, ymodSize(2:end)]);
                ymod_ = nan([nESD, prod(ymodSize(2:end))]);
                rowIndex = 0;
                for j = 1:FixedParams.nZP_size
                    y = ymodZ(j,:);
                    interval = FixedParams.ZPdia_intervals(j:j+1); % size class j interval boundaries
                    ESDj = ESD(interval(1) <= ESD & ESD <= interval(2));
                    nESDj = length(ESDj);
                    rowIndex = max(rowIndex) + (1:nESDj);
                    y_ = repmat(y ./ diff(log10(interval)), [nESDj, 1]);
                    ymod_(rowIndex,:) = y_;
                end
                ymodZ = reshape(ymod_, [nESD, ymodSize(2:end)]); % model output units & number of size classes should now match the data
                clear ymod
                ymod{1} = ymodP;
                ymod{2} = ymodZ;
                %         ymod = [ymodP; ymodZ];
                clear ymod_ ymodP ymodZ
                
                % For each sampling event where size data were measured, extract
                % the relevant trajectories and the times matching the data
                
                % BUT: this messes up the trajectory indices; the result is
                % trajectrory numbers in evTraj not matching with the
                % dimensions of ymod anymore. So we leave this out for now!
                %   %%
                % 
%                 for j = 1:ntrophicLevels
%                     ymod{j} = ymod{j}(:,:,:,modData.sizeFull.evTraj); % filtered out irrelevant sampling events
%                     ymodSize = size(ymod{j});
%                     ymod_ = nan([ymodSize(1:2), ymodSize(4)]); % dimension [size, depth, sample event]
%                     for jj = 1:size(eventTime, 1)
%                         ymod_(:,:,jj) = ymod{j}(:,:,eventTime(j,2),jj);
%                     end
%                     ymod{j} = ymod_;
%                 end
%                 
                % Interpolate modelled depths to derive values corresponding to
                % sampled depths.
                
                % ~~~~
                % here I think there is a mistake in aidans code: event ev is used as a
                % index for the trajectory dimension
                
                
                
                for j = 1:length(sampleEvents)
                    %%
                    ev = sampleEvents(j);
%                     ind0 = modData.sizeFull.Event == ev;
%                     for jj = 1:ntrophicLevels
%                         ymod_ = ymod{jj}(:,:,j); % dimension: [size, depth]
%                         ymod_ = permute(ymod_, [2, 1]); % move depth to first dimension
%                         ind1 = ind0 & strcmp(modData.sizeFull.trophicLevel, trophicLevels{jj});
%                         depths = unique(modData.sizeFull.Depth(ind1));
%                         ymod_ = interp1(FixedParams.z, ymod_, -abs(depths));
%                         ymod_ = permute(ymod_, [2, 1]);
%                         for ij = 1:length(depths)
%                             depth = depths(ij);
%                             ind2 = ind1 & modData.sizeFull.Depth == depth;
%                             modData.sizeFull.(vs)(ind2) = ymod_(:,ij);
%                         end
%                     end
                    
                    % fix: extract all trajs for each event 
                    
                    evTraj = modData.sizeFull.evTraj(:,j); % all 10 trajectories for event j
                    
                    % now, if not all trajs are modellled: 
                    if ntrajOut ~= ntrajData
                        % look up if any modelled trajectory is a member of of evTaj
                        % evTrajIX: (original; before subsetting Forc) indices of modelled trajectories
                        % evTraj: List of (original) trajectory indices per event
                        
                        if ~any(ismember(evTraj, evTrajIX))
                            continue % skip to the next iteration of the loop, because no trajs were modeled for this event
                        end
                        
                        evTraj = find(ismember(evTrajIX, evTraj));  
                    end
                    
                    ind0 = modData.sizeFull.Event == ev;
                    for jj = 1:ntrophicLevels
                        % ymod_ = ymod{jj}(:,:,j); % dimension: [size, depth]
                        ymod_2 = ymod{jj}(:,:,evTraj); % dimension: [size, depth, trajs] 
                         
                 
                        
                        % ymod_ = permute(ymod_, [2, 1]); % move depth to first dimension
                        ymod_2 = permute(ymod_2, [2, 1, 3]); % move depth to first dimension %%
                        
                        % I need a 450x10 matrix as ymod_ in
                        % modData.sizeFull.(vs)(ind2), not a
                        % 450x1 -> so that cost functoin can be
                        % calculatedfor every trajectory
                        
                        ind1 = ind0 & strcmp(modData.sizeFull.trophicLevel, trophicLevels{jj});
                        depths = unique(modData.sizeFull.Depth(ind1));
                        
                        % temp = nan(2,450,10);
                        temp = nan([length(depths), size(ymod_2, 2), 10]);  % dims: no of depths, ESD, trajs / not anymore size(ymod_2, 2, 3), since not always all trajs are modelled 
                       % interpolate to get depths corresponding to
                       % observations, per event trajectory
                       
                       for jjj = 1:length(evTraj) 
                            temp(:,:,jjj) = interp1(FixedParams.z, ymod_2(:,:,jjj), -abs(depths));
                       end 
                       % reassign temp to ymod_2
                       ymod_2 = temp; 
                       
                       ymod_2 = permute(ymod_2, [2, 1, 3]);
                        for ij = 1:length(depths)
                            depth = depths(ij);
                            ind2 = ind1 & modData.sizeFull.Depth == depth;
                            modData.sizeFull.(vs)(ind2,:) = squeeze(ymod_2(:,ij,:)); %% add ,: and squeeze
                        end
                    end
                    %%
                end
            end
        else
            
            
            % THE BELOW STILL NEEDS UPDATED
            
            
            %     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %     % Construction of modData when fitting to fully aggregated size data
            %     %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %     % There is a single size-data vector (for each trophic level) assumed
            %     % to correspond to all sample events (fully aggregated).
            %     % For each sample, extract corresponding (trajectory, time, depth) model
            %     % output. All of these outputs are fitted to the same data vector.
            %
            %
            %     ind = ismember(Data.scalar.Year, unique(Data.size.Year)); % omit sample events from years without size data
            %     ev = unique(Data.scalar.Event(ind)); % all sample events
            %     et = evTraj(:,ev); % sampling events and associated trajectories
            %     nevent = size(et, 2);
            %     etime = nan(nevent, 1); % sample times of each event
            %     for i = 1:nevent
            %         ue = unique(Data.scalar.Yearday(Data.scalar.Event == ev(i)));
            %         etime(i) = ue(1); % this line looks like a quick-fix... check if there's an issue with double-sampling times in prepareFittingData.m
            %     end
            %     depths_obs = [min(Data.size.DepthMin) max(Data.size.DepthMax)]; % sample depth range (the aggregated size data don't specify precise depths)
            %     depth_ind = -FixedParams.zw(2:end) > depths_obs(1,1) & ...
            %         -FixedParams.zw(1:end-1) < depths_obs(1,2); % modelled depth layers corresponding to samples
            %
            %     % VarsSize = Data.size.obsInCostFunction;
            %     allVarsSize = unique(Data.size.dataBinned.Variable);
            %
            %     % The size component of 'modData' is built to match the 'binnedData'
            %     % sub-struct contained in 'Data'.
            %     n = size(Data.size.dataBinned.Year, 1);
            %     modData.size.Variable = cell(n,1);
            %     modData.size.trophicLevel = cell(n,1);
            %     modData.size.Value = nan(n,nsamples);
            %     modData.size.Value_allEvents = nan(n,nevent,nsamples);
            %     modData.size.scaled_Value = nan(n,nsamples);
            %     modData.size.scaled_Value_allEvents = nan(n,nevent,nsamples);
            %
            %     for i = 1:length(allVarsSize)
            %         vs = allVarsSize{i};
            %         ind0 = strcmp(Data.size.dataBinned.Variable, vs);
            %         modData.size.Variable(ind0,:) = Data.size.dataBinned.Variable(ind0);
            %         modData.size.trophicLevel(ind0,:) = Data.size.dataBinned.trophicLevel(ind0);
            %
            %         for j = 1:nsamples
            %             itraj = et(j,:); % trajectory associated with each sampling event
            %             [~, J] = max(sum(out.P(:,depth_ind,FixedParams.PP_Chl_index,etime,itraj))); % use modelled values from chl-max depth layer to compare to data
            %             J = squeeze(J);
            %             switch vs
            %                 case 'CellConc'
            %                     %                 ymod = auxVars.cellDensity(1:1:FixedParams.nPP_size,depth_ind,etime,itraj);
            %                     ymod = auxVars.cellDensity(:,depth_ind,etime,itraj);
            %                     ymod_ = nan(size(ymod,1), nevent);
            %                     for k = 1:nevent
            %                         ymod_(:,k) = ymod(:,J(k,k),k,k);
            %                     end
            %                     ymod = ymod_;
            %                     ymodMean = mean(ymod, 2); % average size-spectra over sampling events
            %                 case 'BioVol'
            %                     %                 ymod = auxVars.biovolume(1:FixedParams.nPP_size,depth_ind,etime,itraj);
            %                     ymod = auxVars.biovolume(:,depth_ind,etime,itraj);
            %                     ymod_ = nan(size(ymod,1), nevent);
            %                     for k = 1:nevent
            %                         ymod_(:,k) = ymod(:,J(k,k),k,k);
            %                     end
            %                     ymod = ymod_;
            %                     ymodMean = mean(ymod, 2); % average size-spectra over sampling events
            %             end
            %
            %             % autotrophs
            %             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            %
            %             modData.size.Value(ind,j) = ymodMean(FixedParams.phytoplankton);
            %             modData.size.Value_allEvents(ind,:,j) = ymod(FixedParams.phytoplankton,:);
            %
            %             modData.size.scaled_Value(ind,j) = Data.size.(['scaleFun_' vs])( ...
            %                 Data.size.dataBinned.scale_mu(ind), ...
            %                 Data.size.dataBinned.scale_sig(ind), ymodMean(FixedParams.phytoplankton));
            %             modData.size.scaled_Value_allEvents(ind,:,j) = Data.size.(['scaleFun_' vs])( ...
            %                 Data.size.dataBinned.scale_mu(ind), ...
            %                 Data.size.dataBinned.scale_sig(ind), ymod(FixedParams.phytoplankton,:));
            %
            %             % heterotrophs
            %             ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
            %
            %             if sum(FixedParams.zooplankton) ~= 1
            %                 % Number of modelled heterotrph size classes matches number
            %                 % of partitions in size data (there may only be one size class)
            %                 modData.size.Value(ind,j) = ymodMean(FixedParams.zooplankton);
            %                 modData.size.Value_allEvents(ind,:,j) = ymod(FixedParams.zooplankton,:);
            %             else
            %                 modData.size.Value(ind,j) = repmat(ymodMean(FixedParams.zooplankton), [sum(ind) 1]);
            %                 modData.size.Value_allEvents(ind,:,j) = repmat(ymod(FixedParams.zooplankton,:), [sum(ind) 1]);
            %             end
            %
            %             modData.size.scaled_Value(ind,j) = Data.size.(['scaleFun_' vs])( ...
            %                 Data.size.dataBinned.scale_mu(ind), ...
            %                 Data.size.dataBinned.scale_sig(ind), ymodMean(FixedParams.zooplankton));
            %             modData.size.scaled_Value_allEvents(ind,:,j) = Data.size.(['scaleFun_' vs])( ...
            %                 Data.size.dataBinned.scale_mu(ind), ...
            %                 Data.size.dataBinned.scale_sig(ind), ymod(FixedParams.zooplankton,:));
            %         end
            %     end
        end
        
end


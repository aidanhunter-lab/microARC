function modData = matchModOutput2Data(out, auxVars, Data, FixedParams)
% Returns model output that 'matches' the data used for optimisation.
% Output struct 'modData' has the same form as the 'Data' struct.

% Model output times are matched to sample times; outputs are interpolated
% over depths to match sampled depth layers.
% Then model outputs are transformed using the same functions that
% standardised the data.

ntrajOut = size(out.N, ndims(out.N));
ntrajData = size(Data.scalar.EventTraj, 2);

if ntrajOut ~= ntrajData
    warning('Number of trajectories in model output "out" does not equal number of trajectories used in from the Data.')
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
    itraj = evTraj(:,i); % trajectories used for event i
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
        modData.scalar.Value(ind,:) = ymod;
        modData.scalar.scaled_Value(ind,:) = ymod_scaled;
    end
end

% Input names of measurements not included in model/cost function
modData.scalar.Variable(~Data.scalar.inCostFunction) = ... 
    Data.scalar.Variable(~Data.scalar.inCostFunction);


%% Size spectra

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




    
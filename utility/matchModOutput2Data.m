function modData = matchModOutput2Data(out, auxVars, Data, FixedParams)
% Returns model output that 'matches' the data used for optimisation.
% Output struct 'modData' has the same form as the 'Data' struct.

% Model output times are matched to sample times; outputs are interpolated
% over depths to match sampled depth layers.
% Then model outputs are transformed using the same functions that
% standardised the data.


%% Scalar data
Vars = Data.scalar.obsInCostFunction;

evTraj = Data.scalar.evTraj; % trajectories used for each sample event
nsamples = size(evTraj, 1); % number of trajectories per sampling event

nEvent = Data.scalar.nEvents; % number of sample events
depths_mod = abs(FixedParams.z); % modelled depth layers

% Build modData to match the form of Data -- we don't need to include every
% field, only the fields required for the cost function.
modData.scalar.Yearday = nan(Data.scalar.nSamples,1);
modData.scalar.Depth = modData.scalar.Yearday;
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
    Variable = Data.scalar.Variable(iEvent);
    modData.scalar.Yearday(iEvent) = Yearday;
    modData.scalar.Depth(iEvent) = Depth;
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
        ymod = interp1(depths_mod, ymod, depth); % interpolate modelled output to match observation depths
        ymod_scaled = Data.scalar.(['scaleFun_' jvar])(Data.scalar.scale_mu(ind), ...
            Data.scalar.scale_sig(ind), ymod); % scale model output using same functions that scaled the data
        modData.scalar.Value(ind,:) = ymod;
        modData.scalar.scaled_Value(ind,:) = ymod_scaled;
    end
end


%% Size spectra

ind = ismember(Data.scalar.Year, unique(Data.size.Year)); % index relavent sampling events
ev = unique(Data.scalar.Event(ind));
et = evTraj(:,ev); % sampling events and associated trajectories
nevent = size(et, 2);
% sample times of each event
etime = nan(nevent, 1);
for i = 1:nevent
    ue = unique(Data.scalar.Yearday(Data.scalar.Event == ev(i)));
    etime(i) = ue(1); % this line looks like a quick-fix... check if there's an issue with double-sampling times in prepareFittingData.m
end
depths_obs = [min(Data.size.DepthMin) max(Data.size.DepthMax)]; % sample depth range
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
        itraj = et(j,:); % trajectories associated with sampling event j
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


% Full size data -- grouped by origin then averaged over depths & events.
% Factors are Variable, trophicLeevel, and waterMass

waterMasses = unique(Data.sizeFull.dataBinned.groupedByOrigin.waterMass);

atl = strcmp(Data.sizeFull.waterMass, 'Atlantic');
arc = strcmp(Data.sizeFull.waterMass, 'Arctic');

trajAtl = Data.sizeFull.evTraj(:,atl);
trajArc = Data.sizeFull.evTraj(:,arc);
trajAtl = trajAtl(:); % All trajectories used for size spectra measurements
trajArc = trajArc(:);

timeAtl = unique(Data.sizeFull.dataBinned.groupedByOrigin.Yearday( ...
    strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, 'Atlantic')));
timeArc = unique(Data.sizeFull.dataBinned.groupedByOrigin.Yearday( ...
    strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, 'Arctic')));

depthAtl = unique(Data.sizeFull.dataBinned.groupedByOrigin.Depth( ...
    strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, 'Atlantic')));
depthArc = unique(Data.sizeFull.dataBinned.groupedByOrigin.Depth( ...
    strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, 'Arctic')));

% VarsSize = Data.size.obsInCostFunction;
allVarsSize = unique(Data.sizeFull.dataBinned.groupedByOrigin.Variable);

n = size(Data.sizeFull.dataBinned.groupedByOrigin.Year, 1);
modData.sizeFull.Variable = cell(n,1);
modData.sizeFull.trophicLevel = cell(n,1);
modData.sizeFull.waterMass = cell(n,1);
modData.sizeFull.Value_Atlantic = nan(n,length(trajAtl));
modData.sizeFull.Value_Arctic = nan(n,length(trajArc));

for i = 1:length(allVarsSize)
    vs = allVarsSize{i};
    ind0 = strcmp(Data.sizeFull.dataBinned.groupedByOrigin.Variable, vs); % index variable
    modData.sizeFull.Variable(ind0,:) = Data.sizeFull.dataBinned.groupedByOrigin.Variable(ind0);
    modData.sizeFull.trophicLevel(ind0,:) = Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel(ind0);
    modData.sizeFull.waterMass(ind0,:) = Data.sizeFull.dataBinned.groupedByOrigin.waterMass(ind0);

    for w = 1:length(waterMasses)
        wm = waterMasses{w};
        traj = eval(['traj' wm(1:3)]);
        time = eval(['time' wm(1:3)]);
        ind1 = ind0 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, wm); % index variable and waterMass
        switch vs
            case 'CellConc'
                ymod = auxVars.cellDensity(:,:,time,traj);
            case 'BioVol'
                ymod = auxVars.biovolume(:,:,time,traj);
        end
        ymod = permute(ymod, [2 1 3 4]); % move depth dimension to front [depth, size, ntraj]
        ymodi = interp1(FixedParams.z , ymod, -eval(['depth' wm(1:3)])); % interpolate model output to match observed depth
        ymodi = squeeze(ymodi); % dimension: [size, trajectories]

        % autotrophs
        ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'autotroph');
        modData.sizeFull.(['Value_' wm])(ind,:) = ymodi(FixedParams.phytoplankton,:);

        % heterotrophs
        ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'heterotroph');
        if FixedParams.nZP_size > 1
            modData.sizeFull.(['Value_' wm])(ind,:) = ymodi(FixedParams.zooplankton,:);
        else
            modData.sizeFull.(['Value_' wm])(ind,:) = repmat(ymodi(FixedParams.zooplankton,:), [sum(ind), 1]);
        end
    end
end


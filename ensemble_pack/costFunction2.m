function [cost, costComponents, costFunctionChoices, costLabel_dataType] = costFunction2(varargin)


% this version is edited, for Vanessas ensemble runs. 


% Labels for cost function options -- code & comments for each below.
% To use a different cost function simply write its label here then include
% another case below -- it must be fully implementable from Data and
% modData. Note below the requirement of cost function label when fitting
% to full (unbinned) size spectra.
costFunctionChoices = { ...
    'syntheticLikelihood_scalar';...
    'IQD_vectorFullSpectrum'; ...               % IQD for FullSpectrum vector data
    'IQD_S_vectorFullSpectrum';...              % IQD for FullSpectrum vector data, S Spectra
    'RMS_scalar'; ...                           % RMS, but only for scalar data (binned)
    'Hellinger_vector'; ...                     % Hellinger, but only for vector data (binned)
    'Hellinger_groupWaterOrigin_vector'; ...    % Hellinger, but only for vector data (binned)
    'Hellinger_vectorFullSpectrum';... 
    } ;      

% Specify type of size data used for each cost function choice: either
% 'binnedSizeData' or 'fullSizeSpectra'.
% To do this automatically we require the cost function label to contain
% the string 'FullSpectrum' iff it's designed to fit to the unbinned size
% spectra data.
dataTypes = {'binnedSizeData', 'fullSizeSpectra'};
di = 1 + cellfun(@(z) contains(z, 'FullSpectrum'), costFunctionChoices);
costLabel_dataType = cell(length(costFunctionChoices), 2);
costLabel_dataType(:,1) = costFunctionChoices;
costLabel_dataType(:,2) = arrayfun(@(z) dataTypes(z), di);
costLabel_dataType = cell2table(costLabel_dataType, ... 
    'VariableNames', {'label','dataType'});

% Evaluating cost requires passing name-value pairs for 'label', 'Data' and
% 'modData' as optional varargin values. The 'label' specifies which cost
% function to use; 'Data' are measurements; 'modData' are the modelled
% equivalents.
extractVarargin(varargin)

calculateCost = exist('Data', 'var') && exist('modData', 'var') && exist('label', 'var');
if ~calculateCost
    cost = nan;
    costComponents = nan;
    return
end

label = eval('label');
if ~ismember(label, costFunctionChoices)
    error('Cost function name must match one of the available options in costFunction.m.')
end

% Data types used to fit the model
Vars = Data.scalar.obsInCostFunction;
VarsSize = Data.size.obsInCostFunction;

%% Evaluate cost function with name given by 'label'
switch label

        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        case 'syntheticLikelihood_scalar' % fit every data point...
            
            
         % THIS IS ADJUSTED FROM AIDANS 'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet'
         % which got lost (moved into comments) in the commit 'Option for fitting to subsets of trajectories -- Atlantic or Arctic '  on 17 May 2021 
         
         % removed the dirichlet part for vector data
         
         
       % Model misfit to data described using a 'synthetic likelihood',
        % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
        % Standardised scalar data are approximately normally distributed.
        % A variety of trajectory combinations are used to generate model
        % outputs, which are transformed identically to the data. Model
        % outputs define normal distributions used as likelihood terms --
        % what is likelihood of observing data given model expectations?
        % Running the model over multiple forcing data trajectories
        % generates the output variability (process error which is probably
        % underestimated).
        
        % Scalar data
        % Each standardised data point is assigned an independent normal 
        % distribution parameterised using model outputs over multiple 
        % trajectories. The likelihood is the product of probabilities of 
        % observing all data points given these model-estimated 
        % distributions.
%        log2pi = log(2*pi);
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             n = length(yobs);
            % n = size(ymod, 2);
            ymod_nonans = ymod(:,all(~isnan(ymod))); % fix is needed for when ntrajs_mod~=ntrajs_obs
            n = size(ymod_nonans,2);
            mu_obs = mean(yobs);
            v_obs = var(yobs);
            Mu = mean(ymod);
            V = var(ymod); 
            Sig = std(ymod); 
            
            % likelihood of mu_obs
            mu = mean(Mu, 'omitnan');
            v = var(Mu, 'omitnan');
            
            Lik_var = prod( (1 ./ sqrt(2*pi*Sig.^2)).* exp(-1 .* ((mu_obs - Mu).^2 ./ (2.*Sig.^2))), 'omitnan' )^(1/n); 
        %    Lik_var = nthroot(prod( (1 ./ sqrt(2*pi*Sig.^2)).* exp(-1 .* ((mu_obs - Mu).^2 ./ (2.*Sig.^2))), 'omitnan' ), n); % avoid imaginary numbers when using the nth root
        
           
%             negLogLik_mu = log2pi + log(v) + (mu_obs - mu) .^ 2 ./ v;
%             % likelihood of v_obs
%             mu = mean(V);
%             v = var(V);
%             negLogLik_v = log2pi + log(v) + (v_obs - mu) .^ 2 ./ v;
%             % combine
%             negLogLik = 0.5 * n * (negLogLik_mu + negLogLik_v);


            L.(varLabel) = Lik_var;
        end
        

        
        costComponents = L;

        
        Lik = double.empty; 
        for i = 1:length(Vars)
            Lik = [Lik costComponents.(Vars{i})];
        end
        Lik = prod(Lik, 'omitnan'); 
        cost = -1 * log(Lik); 
        
        

            
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
        case 'IQD_vectorFullSpectrum'
        
        % cost function for full spectrum size data per station/watermass -> Data.sizeFull and
        % modData.sizeFull are used here
        
        % Integrated Quadratic Distance
        
        
        % Copied then adjusted from Aidans RMS_HellingerFullSpectrum
        trophicLevels = unique(Data.sizeFull.trophicLevel);
        groupedByWaterOrigin = isfield(Data.sizeFull, 'waterMass');
        if groupedByWaterOrigin
            waterMasses = unique(Data.sizeFull.waterMass);
        else
            waterMasses = {[]};
        end
        % VarsSize labels match binned data => redefine for full size spectra
        switch Data.sizeFull.obsInCostFunction{1}
            case 'BioVol'; VarsSize = {'BioVolDensity'};
            case 'CellConc'; VarsSize = {'cellDensity'};
        end
        allEvents = unique(Data.sizeFull.Event, 'stable');
        nVars = length(VarsSize);
        nWaterMasses = length(waterMasses);
        nTrophicLevels = length(trophicLevels);
        % HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
        IQDs = nan(nVars, nWaterMasses, nTrophicLevels); % -> size(IQDs) is now 1,3,2
       
        % at this point we don't need an abundance metric, we just need
        % metric for normalised spectra
%         totAbnMetric = nan(nVars, nWaterMasses, nTrophicLevels);
%         a = log(3) / log(2); % steepness of cost metric for total abundance
        
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind0 = ismember(Data.sizeFull.Event, allEvents(strcmp(Data.sizeFull.waterMass, wm)));
                Events = unique(Data.sizeFull.Event(ind0));
                counter = 0;
                for ie = 1:length(Events)
                    event = Events(ie);
                    ind1 = ind0 & Data.sizeFull.Event == event;
                    depths = unique(Data.sizeFull.Depth(ind1));
                    
                    for id = 1:length(depths)
                        counter = counter + 1;
                        depth = depths(id);
                        ind2 = ind1 & Data.sizeFull.Depth == depth;
                        
                        for it = 1:length(trophicLevels)
                            trophicLevel = trophicLevels{it};
                            ind3 = ind2 & strcmp(Data.sizeFull.trophicLevel, trophicLevel);
                            
                            % filter for observed and modelled spectra
                            yobs = Data.sizeFull.(varLabel)(ind3);
                         %   ymod = modData.sizeFull.(varLabel)(ind3);
                            ymod = modData.sizeFull.(varLabel)(ind3,:); % all 10 trajectories 
%                             % Hellinger distance for relative abundance-at-size                            yobsRel = yobs ./ sum(yobs);
%                             yobsRel = yobs ./ sum(yobs);
%                             ymodRel = ymod ./ sum(ymod);
%                             HellingerDistance = (1 - sum((yobsRel .* ymodRel) .^ 0.5)) .^ 0.5;
%                             HellingerDistances(i,w,it,counter) = HellingerDistance;
                            % IQD for relative abundance or biovolume at
                            % size
                            
                            % normalise observed and modelled spectra
                            yobsRel = yobs ./ sum(yobs);
%                             ymodRel = ymod ./ sum(ymod);
                            ymodRel = nan(size(ymod)); % for all trajectories
                            for iy = 1:size(ymod, 2)
                                ymodRel(:,iy) = ymod(:,iy) ./ sum(ymod(:,iy));
                            end
                            
                       
                            % plot spectra und normalised spectra,
                            % check of everything looks ok                        
%                             figure
%                                 loglog(unique(Data.sizeFull.ESD), yobs) % obs
%                                 hold on
%                                 for iy = 1:size(ymod, 2)
%                                     loglog(unique(Data.sizeFull.ESD), ymod(:,iy), 'r')
%                                 end
%                                 hold off
%                             
%                             figure
%                                 loglog(unique(Data.sizeFull.ESD), yobsRel) % obs
%                                 hold on
%                                 for iy = 1:size(ymod, 2)
%                                     loglog(unique(Data.sizeFull.ESD), ymodRel(:,iy), 'r')
%                                 end
%                                 hold off
                            
                            % calculate IQD
                            ds = diff(log10(Data.sizeFull.ESD(1:2))); 
                            
                            for iy = 1:size(ymodRel, 2)
                                IQD(iy) = sum( (cumsum(ymodRel(:,iy)) - cumsum(yobsRel)).^2 ) * ds ;
                            end
%                             figure
%                                 loglog(unique(Data.sizeFull.ESD), cumsum(yobsRel)) % obs
%                                 hold on
%                                 for iy = 1:size(ymod, 2)
%                                     loglog(unique(Data.sizeFull.ESD), cumsum(ymodRel(:,iy)), 'r')
%                                 end
%                                 hold off

                            %  assign mean IQD to IQDs
                            IQDs(i,w,it,counter) = mean(IQD, 'omitnan');  % here there is now a new dimension to IQDs, when counter exceeds one. but: empty elements are filled with 0, not NaN, this affects the mean calculated from IQDs
                            
                            % skip total abundance metric
%                             % Metric for total abundance
%                             yobsTot = trapz(yobs);
%                             ymodTot = trapz(ymod);
%                             u = abs(log(ymodTot / yobsTot));
%                             u = exp(-a .* u);
%                             z = (1 - u) ./ (1 + u);
%                             totAbnMetric(i,w,it,counter) = z;
                        end
                    end
                end
            end
        end
        
        % testing:
        % when counter > 1 (more than 1 depth in a station), a new 
        % dimension is added to IQDs and empty spots are filled with
        % zero; this affects the means in IQDs_...
        % therefore replace 0 values with NaN, as realistically IQD will
        % never be 0
        IQDs(IQDs==0) = NaN;

%         % Average Hellinger distances over sample events & depths
%         HellingerDistances_ = mean(HellingerDistances, ndims(HellingerDistances));
%         totAbnMetric_ = mean(totAbnMetric, ndims(totAbnMetric));
        
        % Average Integrated Quadratic Distances over sample events & depths
        
        % TESTING: CHANGE BACK FOR FULL TRAJ RUNS! THIS IS ONLY FOR
        % PARTIAL MODEL RUNS (fewer trajectories)
       % IQDs_ = mean(IQDs, ndims(IQDs)); % full runs
       % IQDs_ = IQDs; % partial runs
        IQDs_ = mean(IQDs, ndims(IQDs), 'omitnan'); % should work for both

        % Store in costComponents
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for it = 1:length(trophicLevels)
                trophicLevel = trophicLevels{it};
                if groupedByWaterOrigin
                    for w = 1:length(waterMasses)
                        waterMass = waterMasses{w};
                        label = [varLabel '_' waterMass '_' trophicLevel];
                        
                        % a slash, as used in 'Actic/Atlantic' can not be
                        % used in a fieldname. fix this:
                        
                        label = matlab.lang.makeValidName(label); 
                        
                        
%                         costComponents.([label '_Rel']) = HellingerDistances_(i,w,it);
%                         costComponents.([label '_Tot']) = totAbnMetric_(i,w,it);
                        costComponents.([label '_Rel']) = IQDs_(i,w,it);
                    end
                else
                    label = [varLabel '_' trophicLevel];
%                     costComponents.([label '_Rel']) = HellingerDistances_(i,1,it);
%                     costComponents.([label '_Tot']) = totAbnMetric_(i,1,it);
                    costComponents.([label '_Rel']) = IQDs_(i,1,it);
                end
            end
        end
        
        %%% no weighting, just averanging!
%         % Within each size data group, weight relative abundance-at-size
%         % relative to total abundance.
%         weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' matlab.lang.makeValidName(wm) '_'];
                else
                    cLabel = [varLabel '_'];
                end
%                 cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
%                 cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
%                 costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
                costSize(j,1,i) = cra;
                costSize(j,2,i) = crh;
            end
        end
        
   
%         
%         cost = mean(cost); % finally, average cost over nutrient and size data components
         cost = mean(costSize(:), 'omitnan');
        
         
                 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
    case 'IQD_S_vectorFullSpectrum'
        % cost function for S spectrum size data per station -> Data.size and
        % modData.size are used here
        
        % throw error if optional agrument seasonConfig is not declared
        % this is needed because we need to know whether to use S1/S2 or
        % S3/S4
        
        
%         % BUT: throw error if there is no List with watermasses per
%         % stations somewhere
%         
% %         if ~exist('seasonConfig')
%             error("Set an option for argument 'seasonConfig' for the cost function 'IQD_S_vectorFullSpectrum' to work.")
%         end

        
        % Integrated Quadratic Distance
        
        
        % Copied then adjusted from RMS_HellingerFullSpectrum
        trophicLevels = unique(Data.size.trophicLevel);
        
        scenarios = unique(Data.size.scenario);
        
        % VarsSize labels match binned data => redefine for full size spectra
        switch Data.sizeFull.obsInCostFunction{1}
            case 'BioVol'; VarsSize = {'BioVolDensity'};
            case 'CellConc'; VarsSize = {'cellDensity'};
        end
        
        nVars = length(VarsSize);
        nScenarios = length(scenarios);
        nTrophicLevels = length(trophicLevels);
        % HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
        IQDs = nan(nVars, nScenarios, nTrophicLevels);
               
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for iSc = 1:length(scenarios)
                
                sc = scenarios{iSc};
                ind0 = ismember(Data.size.scenario, sc);
                
                for it = 1:length(trophicLevels)
                    trophicLevel = trophicLevels{it};
                    ind1 = ind0 & strcmp(Data.size.trophicLevel, trophicLevel);

                    % filter for observed and modelled spectra
                    yobs = Data.size.(varLabel)(ind1);           %% filter yobs for scenarios -> S1/S2 or S3/S4, depending on season config
                    
                    ymod = modData.size.(varLabel)(ind1);

 
                            
                    % IQD for relative abundance or biovolume at
                    % size
                            
                    % normalise observed and modelled spectra
                    yobsRel = yobs ./ sum(yobs);
                    ymodRel = ymod ./ sum(ymod);
                          
                       
%                     % erstmal spectra und normalised spectra plotten,
%                     % gucken ob alles ok aussieht                        
%                     figure
%                         loglog(unique(Data.sizeFull.ESD), yobs) % obs
%                         hold on
%                         for iy = 1:size(ymod, 2)
%                             loglog(unique(Data.sizeFull.ESD), ymod(:,iy), 'r')
%                         end
%                         hold off
% 
%                     figure
%                         loglog(unique(Data.sizeFull.ESD), yobsRel) % obs
%                         hold on
%                         for iy = 1:size(ymod, 2)
%                             loglog(unique(Data.sizeFull.ESD), ymodRel(:,iy), 'r')
%                         end
%                         hold off
                            
                    % calculate IQD
                    ds = diff(log10(Data.size.ESD(1:2))); 

                    
                    IQD = sum( (cumsum(ymodRel) - cumsum(yobsRel)).^2 ) * ds ;
                    
%                     figure
%                         loglog(unique(Data.sizeFull.ESD), cumsum(yobsRel)) % obs
%                         hold on
%                         for iy = 1:size(ymod, 2)
%                             loglog(unique(Data.sizeFull.ESD), cumsum(ymodRel(:,iy)), 'r')
%                         end
%                         hold off

                    %  assign mean IQD to IQDs
                    IQDs(i,iSc,it) = IQD;


                end
            end
        end

        

        % Store in costComponents
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for it = 1:length(trophicLevels)
                trophicLevel = trophicLevels{it}; 
                for iSc = 1:length(scenarios)
                    scenario = scenarios{iSc};
                    label = [varLabel '_' scenario '_' trophicLevel];

                    % a slash, as used in 'Actic/Atlantic' can not be
                    % used in a fieldname. fix this:

                    label = matlab.lang.makeValidName(label); 


%                         costComponents.([label '_Rel']) = HellingerDistances_(i,w,it);
%                         costComponents.([label '_Tot']) = totAbnMetric_(i,w,it);
                    costComponents.([label '_Rel']) = IQDs(i,iSc,it);
                end

            end
        end
        
%         
%         cost = mean(cost); % finally, average cost over size data components
         cost = mean(struct2array(costComponents));
        
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        % --- not sure if cost function options below here still work % ---
            
        case 'Hellinger_vector'
            
            % adjusted from 'RMS_Hellinger' to only include vector data
        
       
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
            end
        end
        
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
%         costNutrient = zeros(1,length(Vars));
%         for i = 1:length(Vars)
%             costNutrient(i) = costComponents.(Vars{i});
%         end
%         
%         % Treat PON & POC data as a single type, POM, thereby
%         % downweighting their combined cost contribution
%         POMi = ismember(Vars, {'PON','POC'});
%         costPOM = mean(costNutrient(POMi));
%         costNutrient(POMi) = [];
%         costNutrient = [costNutrient, costPOM];
%         
%         % Average the cost across data-types (nutrient & size)
%         cost = [mean(costNutrient), mean(costSize(:))];
%         % Assign size vs nutrient weighting
%         weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
%         weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
%         cost = weights .* cost;
        
%         cost = mean(cost); % finally, average cost over nutrient and size data components
%         %             cost = cost(2); % try fitting only to the size data
        cost = mean(costSize);
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        case 'Hellinger_groupWaterOrigin_vector'
            
        % unsure if this code section still works
        % only used this for testing a while ago
        % __________________________________
        
        
        % To improve the compatability of cost components between different
        % data types use Hellinger distances for scalar (nutrient) and vector
        % (size) data. The scalar data is already transformed to standard
        % normal distribution so using a distribution-based cost metric
        % ought to be workable...
        % In this function a different, but comparable, metric is used for
        % the total abundance info contained in the size data. This is
        % necessary because Hellinger distance is not applicable to these
        % single data points. The metric also returns values in the
        % interval [0,1], with 0 represetning a perfect fit and 1 an very
        % poor fit.        
        % THERE ARE ISSUES WITH THIS FUNCTION: PDFs FOR THE SCALAR DATA ARE
        % CALCULATED NUMERICALLY, WHICH REQUIRES CHOSING THE RESOLUTION OF
        % THE GRID OVER WHICH PDFs ARE GENERATED. THE GRID RESOLUTION 
        % (ARBITRARY CHOICE) AFFECTS THE COST... NOT ACCEPTABLE!
        % A POSSIBLE SOLUTION, KEEPING THE 'SPIRIT' OF THIS METHOD, WOULD
        % BE TO USE THE CDFs OF THE SCALAR DATA RATHER THAN PDFs.
        % INTERPOLATION WOULD BE REQUIRED TO ALIGN THE MODELLED CDFs WITH
        % OBSERVED CDFs, AND THE METRIC WOULD BE SOME AVERAGE BETWEEN-CDF
        % DISTANCE, RATHER THAN HELLINGER DISTANCE, BUT IT WOULD BE A
        % CLEANER (LESS BIASED BY ARBITRARY CHOICE) METHOD...
        
        
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
%         costNutrient = zeros(1,length(Vars));
%         for i = 1:length(Vars)
%             costNutrient(i) = costComponents.(Vars{i});
%         end
%         
%         % Average the cost across data-types (nutrient & size)
%         cost = [mean(costNutrient), mean(costSize(:))];
%         % Assign size vs nutrient weighting
%         weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
%         weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
%         cost = weights .* cost;
%         
%         cost = mean(cost); % finally, average cost over nutrient and size data components
        cost = mean(costSize); % finally, average cost over nutrient and size data components
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    case 'IQD_vector' 
        
        % unsure if this code section still works
        % __________________________________
        
        % IQD for binned vector data (per sample)
        % adjusted from IQD_vectorFullSpectrum
        
        % again: Data.sizeFull contains spectra per sample (.size only for
        % scenario 3 (autumn-arctic samples...)) and also contains a binned
        % subset which will be used here. modData must also have a .sizeFull
        % substruct, but then misses dataBinned. we will therefore
        % construct it first in matchModOutput2Data2.m 
      
        
        trophicLevels = unique(Data.sizeFull.trophicLevel);
        groupedByWaterOrigin = isfield(Data.sizeFull, 'waterMass');
        if groupedByWaterOrigin
            waterMasses = unique(Data.sizeFull.waterMass);
        else
            waterMasses = {[]};
        end
        % VarsSize labels match binned data => redefine for full size spectra
        switch Data.sizeFull.obsInCostFunction{1}
            case 'BioVol'; VarsSize = {'BioVolDensity'};
            case 'CellConc'; VarsSize = {'cellDensity'};
        end
        allEvents = unique(Data.sizeFull.Event, 'stable');
        nVars = length(VarsSize);
        nWaterMasses = length(waterMasses);
        nTrophicLevels = length(trophicLevels);
        % HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
        IQDs = nan(nVars, nWaterMasses, nTrophicLevels);
       
        % at this point we don't need an abundance metric, we just need
        % metric for normalised spectra
%         totAbnMetric = nan(nVars, nWaterMasses, nTrophicLevels);
%         a = log(3) / log(2); % steepness of cost metric for total abundance
        
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind0 = ismember(Data.sizeFull.Event, allEvents(strcmp(Data.sizeFull.waterMass, wm)));
                Events = unique(Data.sizeFull.Event(ind0));
                counter = 0;
                for ie = 1:length(Events)
                    event = Events(ie);
                    ind1 = ind0 & Data.sizeFull.Event == event;
                    depths = unique(Data.sizeFull.Depth(ind1));
                    for id = 1:length(depths)
                        counter = counter + 1;
                        depth = depths(id);
                        ind2 = ind1 & Data.sizeFull.Depth == depth;
                        for it = 1:length(trophicLevels)
                            trophicLevel = trophicLevels{it};
                            ind3 = ind2 & strcmp(Data.sizeFull.trophicLevel, trophicLevel);
                            
                            % filter for observed and modelled spectra
                            yobs = Data.sizeFull.(varLabel)(ind3);
                         %   ymod = modData.sizeFull.(varLabel)(ind3);
                            ymod = modData.sizeFull.(varLabel)(ind3,:); % all 10 trajectories 
%                             % Hellinger distance for relative abundance-at-size                            yobsRel = yobs ./ sum(yobs);
%                             yobsRel = yobs ./ sum(yobs);
%                             ymodRel = ymod ./ sum(ymod);
%                             HellingerDistance = (1 - sum((yobsRel .* ymodRel) .^ 0.5)) .^ 0.5;
%                             HellingerDistances(i,w,it,counter) = HellingerDistance;
                            % IQD for relative abundance or biovolume at
                            % size
                            
                            % normalise observed and modelled spectra
                            yobsRel = yobs ./ sum(yobs);
%                             ymodRel = ymod ./ sum(ymod);
                            ymodRel = nan(size(ymod)); % for all trajectories
                            for iy = 1:size(ymod, 2)
                                ymodRel(:,iy) = ymod(:,iy) ./ sum(ymod(:,iy));
                            end
                           
                            
                            % calculate IQD
                            ds = diff(log10(Data.sizeFull.ESD(1:2))); 
                            
                            for iy = 1:size(ymodRel, 2)
                                IQD(iy) = sum( (cumsum(ymodRel(:,iy)) - cumsum(yobsRel)).^2 ) * ds ;
                            end
%                             figure
%                                 loglog(unique(Data.sizeFull.ESD), cumsum(yobsRel)) % obs
%                                 hold on
%                                 for iy = 1:size(ymod, 2)
%                                     loglog(unique(Data.sizeFull.ESD), cumsum(ymodRel(:,iy)), 'r')
%                                 end
%                                 hold off

                            %  assign mean IQD to IQDs
                            IQDs(i,w,it,counter) = mean(IQD);
                            
                            % skip total abundance metric
%                             % Metric for total abundance
%                             yobsTot = trapz(yobs);
%                             ymodTot = trapz(ymod);
%                             u = abs(log(ymodTot / yobsTot));
%                             u = exp(-a .* u);
%                             z = (1 - u) ./ (1 + u);
%                             totAbnMetric(i,w,it,counter) = z;
                        end
                    end
                end
            end
        end

%         % Average Hellinger distances over sample events & depths
%         HellingerDistances_ = mean(HellingerDistances, ndims(HellingerDistances));
%         totAbnMetric_ = mean(totAbnMetric, ndims(totAbnMetric));
        % Average Integrated Quadratic Distances over sample events & depths
        IQDs_ = mean(IQDs, ndims(IQDs));
        

        % Store in costComponents
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for it = 1:length(trophicLevels)
                trophicLevel = trophicLevels{it};
                if groupedByWaterOrigin
                    for w = 1:length(waterMasses)
                        waterMass = waterMasses{w};
                        label = [varLabel '_' waterMass '_' trophicLevel];
                        
                        % a slash, as used in 'Actic/Atlantic' can not be
                        % used in a fieldname. fix this:
                        
                        label = matlab.lang.makeValidName(label); 
                        
                        
%                         costComponents.([label '_Rel']) = HellingerDistances_(i,w,it);
%                         costComponents.([label '_Tot']) = totAbnMetric_(i,w,it);
                        costComponents.([label '_Rel']) = IQDs_(i,w,it);
                    end
                else
                    label = [varLabel '_' trophicLevel];
%                     costComponents.([label '_Rel']) = HellingerDistances_(i,1,it);
%                     costComponents.([label '_Tot']) = totAbnMetric_(i,1,it);
                    costComponents.([label '_Rel']) = IQDs_(i,1,it);
                end
            end
        end
        
        %%% no weihgting, just averanging!
%         % Within each size data group, weight relative abundance-at-size
%         % relative to total abundance.
%         weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' matlab.lang.makeValidName(wm) '_'];
                else
                    cLabel = [varLabel '_'];
                end
%                 cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
%                 cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
%                 costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
                costSize(j,1,i) = cra;
                costSize(j,2,i) = crh;
            end
        end
        
        % skip scalar data for this method         
%         cost = mean(cost); % finally, average cost over nutrient and size data components
         cost = mean(costSize(:));
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
    
        
            
            
    case 'Hellinger_vectorFullSpectrum'
       % manipulated from 'RMS_HellingerFullSpectrum'
        
        % Same as case RMS_Hellinger but use model output representing the
        % full size spectra rather than binned (integrated) data.
        % Build the size data cost section to separately fit all size
        % spectra -- keep the returned cost components as autotrophs &
        % heterotrophs in Arctic and/or Atlantic waters => average the
        % costs over sampling events and depths within this function.
        

        
        % Vector (size) data
        trophicLevels = unique(Data.sizeFull.trophicLevel);
        groupedByWaterOrigin = isfield(Data.sizeFull, 'waterMass');
        if groupedByWaterOrigin
            waterMasses = unique(Data.sizeFull.waterMass);
        else
            waterMasses = {[]};
        end
        % VarsSize labels match binned data => redefine for full size spectra
        switch Data.sizeFull.obsInCostFunction{1}
            case 'BioVol'; VarsSize = {'BioVolDensity'};
            case 'CellConc'; VarsSize = {'cellDensity'};
        end
        allEvents = unique(Data.sizeFull.Event, 'stable');
        nVars = length(VarsSize);
        nWaterMasses = length(waterMasses);
        nTrophicLevels = length(trophicLevels);
        HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
%         totAbnMetric = nan(nVars, nWaterMasses, nTrophicLevels);
%         a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind0 = ismember(Data.sizeFull.Event, allEvents(strcmp(Data.sizeFull.waterMass, wm)));
                Events = unique(Data.sizeFull.Event(ind0));
                counter = 0;
                for ie = 1:length(Events)
                    event = Events(ie);
                    ind1 = ind0 & Data.sizeFull.Event == event;
                    depths = unique(Data.sizeFull.Depth(ind1));
                    for id = 1:length(depths)
                        counter = counter + 1;
                        depth = depths(id);
                        ind2 = ind1 & Data.sizeFull.Depth == depth;
                        for it = 1:length(trophicLevels)
                            trophicLevel = trophicLevels{it};
                            ind3 = ind2 & strcmp(Data.sizeFull.trophicLevel, trophicLevel);
                            yobs = Data.sizeFull.(varLabel)(ind3);
                            ymod = modData.sizeFull.(varLabel)(ind3);
                            % Hellinger distance for relative abundance-at-size                            yobsRel = yobs ./ sum(yobs);
                            yobsRel = yobs ./ sum(yobs);
                            ymodRel = ymod ./ sum(ymod);
                            HellingerDistance = (1 - sum((yobsRel .* ymodRel) .^ 0.5)) .^ 0.5;
                            HellingerDistances(i,w,it,counter) = HellingerDistance;
%                             % Metric for total abundance
%                             yobsTot = trapz(yobs);
%                             ymodTot = trapz(ymod);
%                             u = abs(log(ymodTot / yobsTot));
%                             u = exp(-a .* u);
%                             z = (1 - u) ./ (1 + u);
%                             totAbnMetric(i,w,it,counter) = z;
                        end
                    end
                end
            end
        end
        
        % Average Hellinger distances over sample events & depths
        HellingerDistances_ = mean(HellingerDistances, ndims(HellingerDistances));
%         totAbnMetric_ = mean(totAbnMetric, ndims(totAbnMetric));
        % Store in costComponents
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for it = 1:length(trophicLevels)
                trophicLevel = trophicLevels{it};
                if groupedByWaterOrigin
                    for w = 1:length(waterMasses)
                        waterMass = waterMasses{w};
                        label = [varLabel '_' waterMass '_' trophicLevel];
                        
                        % a slash, as used in 'Actic/Atlantic' can not be
                        % used in a fieldname. fix this:
                        
                        label = matlab.lang.makeValidName(label); 
                        
                        
                        costComponents.([label '_Rel']) = HellingerDistances_(i,w,it);
%                         costComponents.([label '_Tot']) = totAbnMetric_(i,w,it);
                    end
                else
                    label = [varLabel '_' trophicLevel];
                    costComponents.([label '_Rel']) = HellingerDistances_(i,1,it);
%                     costComponents.([label '_Tot']) = totAbnMetric_(i,1,it);
                end
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
%         weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' matlab.lang.makeValidName(wm) '_'];
                else
                    cLabel = [varLabel '_'];
                end
%                 cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
%                 cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
%                 costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
                costSize(j,1,i) = cra; 
                costSize(j,2,i) = crh;
            end
        end
        
%         costNutrient = zeros(1,length(Vars));
%         for i = 1:length(Vars)
%             costNutrient(i) = costComponents.(Vars{i});
%         end
        
%         % Treat PON & POC data as a single type, POM, thereby
%         % downweighting their combined cost contribution
%         POMi = ismember(Vars, {'PON','POC'});
%         costPOM = mean(costNutrient(POMi));
%         costNutrient(POMi) = [];
%         costNutrient = [costNutrient, costPOM];
%         
%         % Average the cost across data-types (nutrient & size)
%         cost = [mean(costNutrient), mean(costSize(:))];
%         % Assign size vs nutrient weighting
%         weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
%         weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
%         cost = weights .* cost;
%         
%         cost = mean(cost); % finally, average cost over nutrient and size data components
        cost = mean(costSize(:)); 
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end


%%

if ~exist('cost', 'var')
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if ~exist('costComponents', 'var')
    costComponents = nan;
end

function Data = standardiseFittingData(Data, varargin)
% Use linear mixed models to standardise fitting data.
% Returns linear-model coefficients used to rescale data as a function of
% depth and (if significant) sampling event.

extractVarargin(varargin)

obsInCostFunction_scalar = Data.scalar.obsInCostFunction;
obsInCostFunction_size = Data.size.obsInCostFunction;

if isfield(Data, 'scalar')
    dat = Data.scalar;
else
    dat = [];
end

if isfield(Data, 'size')
    sizeDatBinned = Data.size.dataBinned;
else
    sizeDatBinned = [];
end

if isfield(Data, 'sizeFull')
    sizeDatFullBinned = Data.sizeFull.dataBinned;
else
    sizeDatFullBinned = [];
end

dat.waterMass = dat.waterMass(dat.Event);
fields = fieldnames(dat);
nSamples = dat.nSamples;
for i = 1:length(fields)
    if size(dat.(fields{i}),1) ~= nSamples || contains(fields{i}, 'scaleFun')
        dat = rmfield(dat, fields{i});
    end
end

dat = struct2table(dat);

%% Find scaling factors for all fitting data

% Linearise POM data w.r.t. depth by log tranform of dependent variable
% (PON or POC), equivalent to an exponential decay with depth.

dat.log_Value = log(dat.Value);

dPOM = {'PON','POC'};
logy = dPOM;

for jj = 1:length(dPOM)
    
    dv = dPOM{jj};
    d_ind = strcmp(dat.Variable, dv);
    d = dat(d_ind,:);
    
    % Fit linear mixed models
    lme_mods{1} = fitlme(d, 'log_Value~1');
    lme_mods{2} = fitlme(d, 'log_Value~Depth');
    lme_mods{3} = fitlme(d, 'log_Value~Depth+(1|Event)');
    lme_mods{4} = fitlme(d, 'log_Value~Depth+(Depth-1|Event)');
    lme_mods{5} = fitlme(d, 'log_Value~Depth+(Depth|Event)');
    
    % Use AIC score to choose most appropriate model
    test_aic = cellfun(@(x) x.ModelCriterion.AIC, lme_mods);
    whichMod = test_aic==min(test_aic);
    mod = lme_mods{whichMod};    
    clear lme_mods
    
    fm = fitted(mod);
    dat.scale_mu(d_ind) = fm;
    rsd = mod.residuals;
    dat.scale_sig(d_ind) = std(rsd);
    
    Data.scalar.(['scaleFun_' dv]) = @(mu,sig,x) (log(x) - mu) ./ sig; % scaling function
    dat.scaled_Value(d_ind) = Data.scalar.(['scaleFun_' dv]) ...
        (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data

    
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Chlorophyll is treated identically to POM: linearise data with respecct
% to depth by log tranform of dependent variable (chl)

dv = 'chl_a';
logy = [logy 'chl_a'];

d_ind = strcmp(dat.Variable, dv);
d = dat(d_ind,:);

% Fit linear mixed models
lme_mods{1} = fitlme(d, 'log_Value~1');
lme_mods{2} = fitlme(d, 'log_Value~Depth');
lme_mods{3} = fitlme(d, 'log_Value~Depth+(1|Event)');
lme_mods{4} = fitlme(d, 'log_Value~Depth+(Depth-1|Event)');
lme_mods{5} = fitlme(d, 'log_Value~Depth+(Depth|Event)');

% Use AIC score to choose most appropriate model
test_aic = cellfun(@(x)x.ModelCriterion.AIC, lme_mods);
whichMod = test_aic==min(test_aic);
mod = lme_mods{whichMod};
clear lme_mods

fm = fitted(mod);
dat.scale_mu(d_ind) = fm;
rsd = mod.residuals;
dat.scale_sig(d_ind) = std(rsd);

Data.scalar.(['scaleFun_' dv]) = @(mu,sig,x) (log(x) - mu) ./ sig; % scaling function
dat.scaled_Value(d_ind) = Data.scalar.(['scaleFun_' dv]) ...
    (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Nitrogen data is treated differently because it has a different
% functional relationship to depth: linearise nitrogen data with respect to
% depth by log tranform of independent variable (depth).
% Also, as nitrogen measurements were sampled from either side of the
% N-cline there are 2 distinct patterns in the data. Therefore the raw data 
% cannot be linearised by applying a single LMM to the whole data set.
% Instead, group data above and below the N-cline, then fit two separate
% LMMs. This will produce residuals that are approximately normally
% distributed throughout the water column.

dat.log_Depth = log(dat.Depth);

dv = 'N';
logx = {dv};

d_ind = strcmp(dat.Variable, dv);
d = dat(d_ind,:);

if ~exist('NclineDepth', 'var') || isempty(NclineDepth)
    
    % Fit linear mixed models
    lme_mods{1} = fitlme(d, 'Value~1');
    lme_mods{2} = fitlme(d, 'Value~log_Depth');
    lme_mods{3} = fitlme(d, 'Value~log_Depth+(1|Event)');
    lme_mods{4} = fitlme(d, 'Value~log_Depth+(log_Depth-1|Event)');
    lme_mods{5} = fitlme(d, 'Value~log_Depth+(log_Depth|Event)');
    
    % Use AIC score to choose most appropriate model
    test_aic = cellfun(@(x)x.ModelCriterion.AIC, lme_mods);
    whichMod = test_aic==min(test_aic);
    mod = lme_mods{whichMod};
    clear lme_mods
    
    fm = fitted(mod);
    dat.scale_mu(d_ind) = fm;
    rsd = mod.residuals;
    dat.scale_sig(d_ind) = std(rsd);
    
    Data.scalar.(['scaleFun_' dv]) = @(mu,sig,x) (x - mu) ./ sig; % scaling function
    dat.scaled_Value(d_ind) = Data.scalar.(['scaleFun_' dv]) ...
        (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data

elseif isnumeric(NclineDepth) && length(NclineDepth) == 1
    % Account for N-cline by using two separate models
    shallow = abs(d.Depth) <= abs(NclineDepth);
    deep = ~shallow;
    d1 = d(shallow,:);
    d2 = d(deep,:);
    
    % Fit linear mixed models
    lme_mods1{1} = fitlme(d1, 'Value~1');
    lme_mods1{2} = fitlme(d1, 'Value~log_Depth');
    lme_mods1{3} = fitlme(d1, 'Value~log_Depth+(1|Event)');
    lme_mods1{4} = fitlme(d1, 'Value~log_Depth+(log_Depth-1|Event)');
    lme_mods1{5} = fitlme(d1, 'Value~log_Depth+(log_Depth|Event)');

    lme_mods2{1} = fitlme(d2, 'Value~1');
    lme_mods2{2} = fitlme(d2, 'Value~log_Depth');
    lme_mods2{3} = fitlme(d2, 'Value~log_Depth+(1|Event)');
    lme_mods2{4} = fitlme(d2, 'Value~log_Depth+(log_Depth-1|Event)');
    lme_mods2{5} = fitlme(d2, 'Value~log_Depth+(log_Depth|Event)');

    % Use AIC score to choose most appropriate model
    test_aic1 = cellfun(@(x)x.ModelCriterion.AIC, lme_mods1);
    test_aic2 = cellfun(@(x)x.ModelCriterion.AIC, lme_mods2);
    whichMod1 = test_aic1==min(test_aic1);
    whichMod2 = test_aic2==min(test_aic2);
    mod1 = lme_mods1{whichMod1};
    mod2 = lme_mods2{whichMod2};
    clear lme_mods1 lme_mods2
    
    fm1 = fitted(mod1);
    fm2 = fitted(mod2);
    scale_mu = nan(size(d, 1), 1);
    scale_mu(shallow) = fm1;
    scale_mu(deep) = fm2;
    dat.scale_mu(d_ind) = scale_mu;
    
    rsd1 = mod1.residuals;
    rsd2 = mod2.residuals;
    scale_sig = nan(size(d, 1), 1);
    scale_sig(shallow) = std(rsd1);
    scale_sig(deep) = std(rsd2);
    dat.scale_sig(d_ind) = scale_sig;
    
    Data.scalar.(['scaleFun_' dv]) = @(mu,sig,x) (x - mu) ./ sig; % scaling function
    dat.scaled_Value(d_ind) = Data.scalar.(['scaleFun_' dv]) ...
        (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data
    
end


%% Store outputs in the Data struct
Data.scalar.scale_mu = dat.scale_mu;
Data.scalar.scale_sig = dat.scale_sig;
Data.scalar.scaled_Value = dat.scaled_Value;

Data.scalar.scale_mu(~Data.scalar.inCostFunction) = nan;
Data.scalar.scale_sig(~Data.scalar.inCostFunction) = nan;
Data.scalar.scaled_Value(~Data.scalar.inCostFunction) = nan;


%% Standardise size spectra measurements

% There are no covariates available in these data, so standardise by log
% transforming, then substracting the mean and dividing by the standard 
% deviation across size classes.
% Box-Cox transforms would be better than log transforms as the resulting
% data would be closer to standard normal, but this will not influence the
% most effective (I think!) cost functions which use size data on natural
% scale.

% Aggregated size spectra -- no sampling event covariate

vars = unique(sizeDatBinned.Variable, 'stable'); % size spectra measurement types (cell concentration, bio-volume, or nitrogen concentration)
if isfield(sizeDatBinned, 'scenario')
    scenarios = unique(sizeDatBinned.scenario);
else
    scenarios = nan;
end
groups = unique(sizeDatBinned.trophicLevel);

for i = 1:length(vars)
    varLabel = vars{i};
    indi = strcmp(sizeDatBinned.Variable, vars{i});    
    switch varLabel
        % We may use different standardising functions for each data type -- see commented case below where a logit transform was used for biovolume (when the units were appropriate to constrain 0<biovol<1)
        case {'CellConc', 'NConc', 'BioVol'}
            Data.size.(['scaleFun_' varLabel]) = @(mu,sig,x) (log(x)-mu) ./ sig;
%         case 'BioVol'
%             Data.size.(['scaleFun_' varLabel]) = @(mu,sig,x) (log(x ./ (1-x)) - mu) ./ sig;
    end
    for j = 1:length(scenarios)
        if iscell(scenarios)
            indj = indi & strcmp(sizeDatBinned.scenario, scenarios{j});
        else
            indj = indi;
        end
        for k = 1:length(groups)
            indk = indj & strcmp(sizeDatBinned.trophicLevel, groups{k});
            y = sizeDatBinned.Value(indk);
            switch varLabel
                % See above comment about using different functions
                case {'CellConc', 'NConc', 'BioVol'}
                    yl = log(y);
%                 case 'BioVol'
%                     yl = log(y ./ (1-y));
            end
            mu = mean(yl);
            sig = std(yl);
            sizeDatBinned.scale_mu(indk,:) = mu;
            sizeDatBinned.scale_sig(indk,:) = sig;
            sizeDatBinned.scaled_Value(indk,:) = Data.size.(['scaleFun_' varLabel])(mu, sig, y);
        end
    end    
end

Data.size.dataBinned = sizeDatBinned;


% Full size spectra -- here we have a sampling event covariate, and some
% events have measures over multiple depths... but just repeat the process
% used for the aggregated data.

if ~isempty(sizeDatFullBinned)
    
    vars = unique(sizeDatFullBinned.Variable, 'stable'); % size spectra measurement types (cell concentration, bio-volume, or nitrogen concentration)
    cruises = unique(sizeDatFullBinned.Cruise);
    events = unique(sizeDatFullBinned.Event);
    groups = unique(sizeDatFullBinned.trophicLevel);
    
    
    for i = 1:length(vars)
        varLabel = vars{i};
        indi = strcmp(sizeDatFullBinned.Variable, vars{i});
        switch varLabel
            case {'CellConc', 'NConc', 'BioVol'}
                Data.sizeFull.(['scaleFun_' varLabel]) = @(mu,sig,x) (log(x)-mu) ./ sig;
%             case 'BioVol'
%                 Data.sizeFull.(['scaleFun_' varLabel]) = @(mu,sig,x) (log(x ./ (1-x)) - mu) ./ sig;
        end
        for j = 1:length(cruises)
            indj = indi & strcmp(sizeDatFullBinned.Cruise, cruises{j});
            for k = 1:length(groups)
                indk = indj & strcmp(sizeDatFullBinned.trophicLevel, groups{k});
                for je = 1:length(events)
                    inde = indk & sizeDatFullBinned.Event == events(je);
                    depths = unique(sizeDatFullBinned.Depth(inde));
                    for jd = 1:length(depths)
                        indd = inde & sizeDatFullBinned.Depth == depths(jd);
                        y = sizeDatFullBinned.Value(indd);
                        switch varLabel
                            case {'CellConc', 'NConc', 'BioVol'}
                                yl = log(y);
%                             case 'BioVol'
%                                 yl = log(y ./ (1-y));
                        end
                        mu = mean(yl);
                        sig = std(yl);
                        sizeDatFullBinned.scale_mu(indd,:) = mu;
                        sizeDatFullBinned.scale_sig(indd,:) = sig;
                        sizeDatFullBinned.scaled_Value(indd,:) = Data.sizeFull.(['scaleFun_' varLabel])(mu, sig, y);
                    end
                end
            end
        end
    end
    
    Data.sizeFull.dataBinned = sizeDatFullBinned;
    
end


%% Plots of standardised data

% Scalar data

if exist('plotScalarData', 'var') && plotScalarData
    for jj = 1:length(obsInCostFunction_scalar)
        vv = obsInCostFunction_scalar{jj};
        if ismember(vv, logy)
            dv = logy{cellfun(@(x)~isempty(x),regexp(vv,logy))};
            d_ind = strcmp(dat.Variable, dv);
            d = dat(d_ind,:);
            uev = unique(d.Event);
            nev = length(uev);
            % Plotting variables
            xraw = d.Depth; xmod = d.Depth;
            yraw = d.Value; ymod = d.log_Value;
            %             x = d.Depth;
            %             y = d.log_Value;
            % Plot labels
            xlabraw = [dv ' (mmol N m^{-3})']; xlabmod = ['ln(' dv ')']; xlabstd = ['standardised ' dv];
            ylabraw = 'depth (m)'; ylab = 'depth (m)'; ylab2 = 'sampling event';
            %             x_lab = ['ln(' dv ')']; x_lab2 = ['standardised ' dv];
            %             y_lab = 'depth (m)'; y_lab2 = 'sampling event';
        end
        if ismember(vv,logx)
            dv = 'N';
            d_ind = strcmp(dat.Variable, dv);
            d = dat(d_ind,:);
            uev = unique(d.Event);
            nev = length(uev);
            % Plotting model variables
            xraw = d.Depth; xmod = d.log_Depth;
            yraw = d.Value; ymod = d.Value;
            %             x = d.log_Depth;
            %             y = d.Value;
            % Plot labels
            xlabraw = [dv ' (mmol N m^{-3})']; xlabmod = [dv ' (mmol N m^{-3})']; xlabstd = ['standardised ' dv];
            ylabraw = 'depth (m)'; ylab = 'ln(depth)'; ylab2 = 'sampling event';
            %             x_lab = [dv ' (mmol N m^{-3})']; x_lab2 = ['standardised ' dv];
            %             y_lab = 'ln(depth)'; y_lab2 = 'sampling event';
        end
        
        ys = d.scaled_Value; % scaled data
        
        % Compare raw data to standardised data
        figure
        subplot(3,2,1)
        scatter(yraw, xraw, 'k')
        hold on
        ca = gca; xl = ca.XLim; yl = ca.YLim;
%         xl_ = [-0.06 0.06] .* diff(xl) + xl;
        yl_ = [-0.06 0.06] .* diff(yl) + yl;
        mu = mean(yraw);
        s = std(yraw);
        pgon = polyshape([mu-s mu+s mu+s mu-s mu-s], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
        pg = plot(pgon);
        pg.FaceColor = [1 1 1] ./ 255;
        pg.FaceAlpha = 0.2;
        plot(repmat(mu, [1 2]), yl, 'k')
        xlim(xl); ylim(yl);
        xlabel(xlabraw); ylabel(ylabraw);
        hold off
        
        subplot(3,2,2)
        G = cell(nev,1);
        for i = 1:nev
            G{i} = yraw(d.Event==uev(i));
        end
        grp = cell2mat(arrayfun(@(i){i*ones(numel(G{i}),1)},(1:numel(G))'));
        boxplot(vertcat(G{:}),grp,'orientation','horizontal')
        hold on
        xlabel(xlabraw); ylabel(ylab2);
        yticklabels(num2str(uev))
        mu = mean(yraw);
        ca = gca;
        plot(repmat(mu,[1 2]), ca.YLim,'k')
        hold off
        
        subplot(3,2,3)
        scatter(ymod, xmod, 'k')
        hold on
        ca = gca; xl = ca.XLim; yl = ca.YLim;
%         xl_ = [-0.06 0.06] .* diff(xl) + xl;
        yl_ = [-0.06 0.06] .* diff(yl) + yl;
        mu = mean(ymod);
        s = std(ymod);
        pgon = polyshape([mu-s mu+s mu+s mu-s mu-s], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
        pg = plot(pgon);
        pg.FaceColor = [1 1 1] ./ 255;
        pg.FaceAlpha = 0.2;
        plot(repmat(mu, [1 2]), yl, 'k')
        xlim(xl); ylim(yl);
        xlabel(xlabmod); ylabel(ylab);
        hold off
        
        subplot(3,2,4)
        G = cell(nev,1);
        for i = 1:nev
            G{i} = ymod(d.Event==uev(i));
        end
        grp = cell2mat(arrayfun(@(i){i*ones(numel(G{i}),1)},(1:numel(G))'));
        boxplot(vertcat(G{:}),grp,'orientation','horizontal')
        hold on
        xlabel(xlabmod); ylabel(ylab2);
        yticklabels(num2str(uev))
        mu = mean(ymod);
        ca = gca;
        plot(repmat(mu,[1 2]), ca.YLim,'k')
        hold off
        
        subplot(3,2,5)
        scatter(ys, xmod, 'k')
        hold on
        ca = gca; xl = ca.XLim; yl = ca.YLim;
%         xl_ = [-0.06 0.06] .* diff(xl) + xl;
        yl_ = [-0.06 0.06] .* diff(yl) + yl;
        mu = mean(ys);
        s = std(ys);
        pgon = polyshape([mu-s mu+s mu+s mu-s mu-s], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
        pg = plot(pgon);
        pg.FaceColor = [1 1 1] ./ 255;
        pg.FaceAlpha = 0.2;
        plot(repmat(mu, [1 2]), yl, 'k')
        xlim(xl); ylim(yl);
        xlabel(xlabstd); ylabel(ylab);
        hold off
        
        subplot(3,2,6)
        G = cell(nev,1);
        for i = 1:nev
            G{i} = ys(d.Event==uev(i));
        end
        grp = cell2mat(arrayfun(@(i){i*ones(numel(G{i}),1)},(1:numel(G))'));
        boxplot(vertcat(G{:}),grp,'orientation','horizontal')
        hold on
        xlabel(xlabstd); ylabel(ylab2);
        yticklabels(num2str(uev))
        mu = mean(ys);
        ca = gca;
        plot(repmat(mu,[1 2]), ca.YLim,'k')
        hold off
        
        sgtitle(['Compare raw and standardised ' dv ' measurements'])
        
        pause(0.1)
    end
end



% Size data
if exist('plotSizeData', 'var') && plotSizeData
    for jj = 1:length(obsInCostFunction_size)
        vv = obsInCostFunction_size{jj};
        
        d_ind = strcmp(sizeDatBinned.Variable, vv);
        d = structfun(@(x) x(d_ind), sizeDatBinned, 'UniformOutput', false);
        
        y = d.Value; % raw data
        ys = d.scaled_Value; % scaled data
        x = d.size;
        
        trophicLevels = unique(d.trophicLevel);
        ntrophicLevels = length(trophicLevels);
        
        nsizes = length(unique(d.sizeClass));
        nsizesP = length(unique(d.sizeClass(strcmp(d.trophicLevel, 'autotroph'))));
        nsizesZ = length(unique(d.sizeClass(strcmp(d.trophicLevel, 'heterotroph'))));
        y_ = nan(nsizes, ntrophicLevels);
        y_(1:nsizesP, strcmp(trophicLevels, 'autotroph')) = y(strcmp(d.trophicLevel, 'autotroph'));
        y_(1:nsizesZ, strcmp(trophicLevels, 'heterotroph')) = y(strcmp(d.trophicLevel, 'heterotroph'));
        y = y_;
        
        y_ = nan(nsizes, ntrophicLevels);
        y_(1:nsizesP, strcmp(trophicLevels, 'autotroph')) = ys(strcmp(d.trophicLevel, 'autotroph'));
        y_(1:nsizesZ, strcmp(trophicLevels, 'heterotroph')) = ys(strcmp(d.trophicLevel, 'heterotroph'));
        ys = y_;
        %             y = reshape(y, [nsizes, ntrophicLevels]);
        %             ys = reshape(ys, [nsizes, ntrophicLevels]);
        y_ = nan(nsizes, ntrophicLevels);
        y_(1:nsizesP, strcmp(trophicLevels, 'autotroph')) = x(strcmp(d.trophicLevel, 'autotroph'));
        y_(1:nsizesZ, strcmp(trophicLevels, 'heterotroph')) = x(strcmp(d.trophicLevel, 'heterotroph'));
        x = y_;
        %             x = reshape(x, [nsizes, ntrophicLevels]);
        
        Cols = [0 1 0; 1 0 0];
        
        % Compare raw data to standardised data
        figure
        subplot(1,2,1)
        
        lplt = loglog(x, y, 'o-');
        
        set(lplt, {'color'}, num2cell(Cols, 2))
        
        legend(trophicLevels)
        
        hold on
        
%         ca = gca;
%         xl = ca.XLim;
%         yl = ca.YLim;
        
        xlabel('ESD (\mum)');
        
        switch vv
            case 'CellConc'
                ylab = 'raw cell conc. data (cells m^{-3})';
            case 'BioVol'
                ylab = 'raw biovolume data (m^3 m^{-3})';
            case 'NConc'
                ylab = 'raw N conc. data (mmol N m^{-3})';
        end
        
        ylabel(ylab);
        hold off
        
        
        subplot(1,2,2)
        lplt = semilogx(x, ys, 'o-');
        set(lplt, {'color'}, num2cell(Cols, 2))
        
        hold on
        
%         ca = gca;
%         xl = ca.XLim;
%         yl = ca.YLim;
        
        xlabel('ESD (\mum)');
        
        switch vv
            case 'CellConc'
                ylab = 'standardised cell conc. data';
                tl = 'Compare raw and standardised cell conc. measurements';
            case 'BioVol'
                ylab = 'standardised biovolume data';
                tl = 'Compare raw and standardised biovolume measurements';
            case 'NConc'
                ylab = 'standardised N conc. data';
                tl = 'Compare raw and standardised N conc. measurements';
        end
        
        ylabel(ylab);
        hold off
        sgtitle(tl)
        
        pause(0.1)
    end
end


if exist('plotAllData', 'var') && plotAllData
    figure
    alph = 0.5;
    col_PON = [1 0 0];
    col_POC = [1 0.5 0];
    col_Chl = [0 1 0];
    col_N = [0 0 1];
    
    subplot(3,1,1)
    ind = strcmp(dat.Variable, 'PON');
    scatter(dat.scaled_Value(ind),dat.Depth(ind), 'MarkerEdgeColor', col_PON, 'MarkerFaceColor', col_PON, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph);
    hold on
    ylabel('depth (m)')
    ind = strcmp(dat.Variable, 'POC');
    scatter(dat.scaled_Value(ind),dat.Depth(ind), 'MarkerEdgeColor', col_POC, 'MarkerFaceColor', col_POC, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph);
    ind = strcmp(dat.Variable, 'chl_a');
    scatter(dat.scaled_Value(ind),dat.Depth(ind), 'MarkerEdgeColor', col_Chl, 'MarkerFaceColor', col_Chl, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph);
    ind = strcmp(dat.Variable, 'N');
    scatter(dat.scaled_Value(ind),dat.Depth(ind), 'MarkerEdgeColor', col_N, 'MarkerFaceColor', col_N, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph);
    
    ca = gca; xl = ca.XLim; yl = ca.YLim;
    % legend
    text(xl(1)+0.1*diff(xl),yl(1)+0.95*diff(yl),'PON')
    xl_ = xl(1) + [0.025 0.075] * diff(xl);
    yl_ = yl(1) + [0.975 0.925] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = [col_PON alph];
    
    text(xl(1)+0.1*diff(xl),yl(1)+0.85*diff(yl),'POC')
    yl_ = yl(1) + [0.875 0.825] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = [col_POC alph];
    
    text(xl(1)+0.1*diff(xl),yl(1)+0.75*diff(yl),'chl_a')
    yl_ = yl(1) + [0.775 0.725] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = [col_Chl alph];
    
    text(xl(1)+0.1*diff(xl),yl(1)+0.65*diff(yl),'N')
    yl_ = yl(1) + [0.675 0.625] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = [col_N alph];
    
    hold off
    
    subplot(3,1,2)
    ind = strcmp(dat.Variable, 'PON');
    scatter(dat.scaled_Value(ind), dat.Event(ind), 'MarkerEdgeColor', col_PON, 'MarkerFaceColor', col_PON, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph)
    hold on
    ylabel('sampling event')
    ind = strcmp(dat.Variable, 'POC');
    scatter(dat.scaled_Value(ind), dat.Event(ind), 'MarkerEdgeColor', col_POC, 'MarkerFaceColor', col_POC, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph)
    ind = strcmp(dat.Variable, 'chl_a');
    scatter(dat.scaled_Value(ind), dat.Event(ind), 'MarkerEdgeColor', col_Chl, 'MarkerFaceColor', col_Chl, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph)
    ind = strcmp(dat.Variable, 'N');
    scatter(dat.scaled_Value(ind), dat.Event(ind), 'MarkerEdgeColor', col_N, 'MarkerFaceColor', col_N, 'MarkerEdgeAlpha', alph, 'MarkerFaceAlpha', alph)
    hold off
    
    subplot(3,1,3)
    bw = 0.5; % smoothing kernel bandwidth to approximate distributions of scaled data
    
    ind = strcmp(dat.Variable, 'PON');
    y = dat.scaled_Value(ind);
    x_ = [min(y) max(y)];
    yf_PON = fitdist(y,'kernel','BandWidth',bw);
    
    ind = strcmp(dat.Variable, 'POC');
    y = dat.scaled_Value(ind);
    x_(1) = min([min(y), x_]); x_(2) = max([max(y), x_]);
    yf_POC = fitdist(y,'kernel','BandWidth',bw);
    
    ind = strcmp(dat.Variable, 'chl_a');
    y = dat.scaled_Value(ind);
    x_(1) = min([min(y), x_]); x_(2) = max([max(y), x_]);
    yf_Chl = fitdist(y,'kernel','BandWidth',bw);
    
    ind = strcmp(dat.Variable, 'N');
    y = dat.scaled_Value(ind);
    x_(1) = min([min(y), x_]); x_(2) = max([max(y), x_]);
    yf_N = fitdist(y,'kernel','BandWidth',bw);
    
    x_ = linspace(min(x_),max(x_),100);
    
    plot(x_,pdf(yf_PON,x_), 'LineWidth', 2, 'Color', [col_PON alph]);
    hold on
    plot(x_,pdf(yf_POC,x_), 'LineWidth', 2, 'Color', [col_POC alph])
    plot(x_,pdf(yf_Chl,x_), 'LineWidth', 2, 'Color', [col_Chl alph])
    plot(x_,pdf(yf_N,x_), 'LineWidth', 2, 'Color', [col_N alph])
    nd = (2*pi)^(-0.5) * exp(-0.5 .* x_ .^ 2);
    plot(x_,nd, 'LineWidth', 2, 'Color', [0 0 0 alph])
    ca = gca; xl = ca.XLim; yl = ca.YLim;
    xl_ = xl(1) + [0.025 0.075] * diff(xl);
    yl_ = yl(1)+0.85*diff(yl);
    text(xl(1)+0.1*diff(xl),yl_,'standard normal')
    plot(xl_, repmat(yl_, [1 2]), 'k')
    xlabel('standardised data')
    ylabel('sample density')
    hold off
    
end





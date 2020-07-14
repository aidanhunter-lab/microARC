function Data = standardiseFittingData(Data, varargin)

% Use linear mixed models to standardise fitting data.
% Returns linear-model coefficients used to rescale data as a function of
% depth and (if significant) sampling event.

fields = fieldnames(Data);
dat = Data;
dat.waterMass = Data.waterMass(Data.Event);
dat = rmfield(dat, ... 
    [{'nSamples';'nEvents';'EventTraj'}; fields(contains(fields,'scaleFun'))]);
dat = struct2table(dat);

% Extra arguments control whether plots are returned. Name-value pairs with
% options: 'plotScaledPON', 'plotScaledPOC', 'plotScaledN', 'plotAllData'
v = reshape(varargin, [2 0.5*length(varargin)]);

%% Find scaling factors for all fitting data

% Linearise POM data w.r.t. depth by log tranform of dependent variable
% (PON or POC), equivalent to an exponential decay with depth.

dat.log_Value = log(dat.Value);

dPOM = {'PON','POC'};

for jj = 1:length(dPOM)
    
    dv = dPOM{jj};
    d_ind = strcmp(dat.Variable, dv);
    d = dat(d_ind,:);
    
    uev = unique(d.Event);
    nev = length(uev);
    
    % Fit linear mixed models
    lme_mods{1} = fitlme(d, 'log_Value~1');
    lme_mods{2} = fitlme(d, 'log_Value~Depth');
    lme_mods{3} = fitlme(d, 'log_Value~Depth+(1|Event)');
    lme_mods{4} = fitlme(d, 'log_Value~Depth+(Depth-1|Event)');
    lme_mods{5} = fitlme(d, 'log_Value~Depth+(Depth|Event)');
    
    % Use AIC score to choose most appropriate model
    test_aic = cellfun(@(x)x.ModelCriterion.AIC, lme_mods);
    whichMod = find(test_aic==min(test_aic));
    mod = lme_mods{whichMod};    
    clear lme_mods
    
    % Extract fitted coefficients of linear models
    cf = zeros(2,1); % fixed-effect (depth) intercept & gradient
    cr = zeros(nev,2); % random-effect (event) intercept & gradient deviations    
    cf_ = mod.fixedEffects;
    cr_ = mod.randomEffects;    
    if length(cf_) == 2, cf = cf_; else, cf(1) = cf_; end        
    if length(cr_) == nev
        if whichMod ~= 4, cr(:,1) = cr_; else, cr(:,2) = cr_; end
    end
    if length(cr_) == 2*nev, cr = reshape(cr_, [2 nev])'; end

    rsd = mod.residuals;
    stdDev = std(rsd);
    % sqrt(mod.MSE)
    
    % Calculate depth- and event-dependent expected values and store
    % functions used to scale the data and model output
    intercept = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
    gradient = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
    dat.scale_mu(d_ind) = intercept + gradient .* dat.Depth(d_ind);
    dat.scale_sig(d_ind) = stdDev;
    
    Data.(['scaleFun_' dv]) = @(mu,sig,x) (log(x) - mu) ./ sig; % scaling function
    
    dat.scaled_Value(d_ind) = Data.(['scaleFun_' dv]) ... 
        (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data
    
    
%     % Store coefficients and standard deviation used to scale data.
%     dat.scale_intercept(d_ind) = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
%     dat.scale_gradient(d_ind) = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
%     dat.scale_stdDev(d_ind) = repmat(stdDev, [sum(d_ind) 1]);
%     % and the scaled data
%     mu = dat.scale_intercept(d_ind) + dat.scale_gradient(d_ind) .* dat.Depth(d_ind);
%     dat.scaled_Value(d_ind) = (dat.log_Value(d_ind) - mu) ./ dat.scale_stdDev(d_ind);
    
    
    
end

logy = dPOM;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Linearise nutrient data w.r.t. depth by log tranform of independent
% variable (depth)

dat.log_Depth = log(dat.Depth);

dv = 'N';

d_ind = strcmp(dat.Variable, dv);
d = dat(d_ind,:);

uev = unique(d.Event);
nev = length(uev);

% Fit linear mixed models
lme_mods{1} = fitlme(d, 'Value~1');
lme_mods{2} = fitlme(d, 'Value~log_Depth');
lme_mods{3} = fitlme(d, 'Value~log_Depth+(1|Event)');
lme_mods{4} = fitlme(d, 'Value~log_Depth+(log_Depth-1|Event)');
lme_mods{5} = fitlme(d, 'Value~log_Depth+(log_Depth|Event)');

% Use AIC score to choose most appropriate model
test_aic = cellfun(@(x)x.ModelCriterion.AIC, lme_mods);
whichMod = find(test_aic==min(test_aic));
mod = lme_mods{whichMod};
clear lme_mods

% Extract fitted coefficients of linear models
cf = zeros(2,1); % fixed-effect (depth) intercept & gradient
cr = zeros(nev,2); % random-effect (event) intercept & gradient deviations
cf_ = mod.fixedEffects;
cr_ = mod.randomEffects;
if length(cf_) == 2, cf = cf_; else, cf(1) = cf_; end
if length(cr_) == nev
    if whichMod ~= 4, cr(:,1) = cr_; else, cr(:,2) = cr_; end
end
if length(cr_) == 2*nev, cr = reshape(cr_, [2 nev])'; end

rsd = mod.residuals;
stdDev = std(rsd);

% Calculate depth- and event-dependent expected values and store
% functions used to scale the data and model output
intercept = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
gradient = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
dat.scale_mu(d_ind) = intercept + gradient .* dat.log_Depth(d_ind);
dat.scale_sig(d_ind) = stdDev;

Data.(['scaleFun_' dv]) = @(mu,sig,x) (x - mu) ./ sig; % scaling function

dat.scaled_Value(d_ind) = Data.(['scaleFun_' dv]) ...
    (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data


% % Store coefficients and standard deviation used to scale data.
% dat.scale_intercept(d_ind) = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
% dat.scale_gradient(d_ind) = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
% dat.scale_stdDev(d_ind) = repmat(stdDev, [sum(d_ind) 1]);
% % and the scaled data
% mu = dat.scale_intercept(d_ind) + dat.scale_gradient(d_ind) .* dat.log_Depth(d_ind);
% dat.scaled_Value(d_ind) = (dat.Value(d_ind) - mu) ./ dat.scale_stdDev(d_ind);

logx = {dv};


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Linearise chlorophll data w.r.t. depth by log tranform of dependent
% variable (chl)

dv = 'chl_a';

d_ind = strcmp(dat.Variable, dv);
d = dat(d_ind,:);

uev = unique(d.Event);
nev = length(uev);

% Fit linear mixed models
lme_mods{1} = fitlme(d, 'log_Value~1');
lme_mods{2} = fitlme(d, 'log_Value~Depth');
lme_mods{3} = fitlme(d, 'log_Value~Depth+(1|Event)');
lme_mods{4} = fitlme(d, 'log_Value~Depth+(Depth-1|Event)');
lme_mods{5} = fitlme(d, 'log_Value~Depth+(Depth|Event)');

% Use AIC score to choose most appropriate model
test_aic = cellfun(@(x)x.ModelCriterion.AIC, lme_mods);
whichMod = find(test_aic==min(test_aic));
mod = lme_mods{whichMod};
clear lme_mods

% Extract fitted coefficients of linear models
cf = zeros(2,1); % fixed-effect (depth) intercept & gradient
cr = zeros(nev,2); % random-effect (event) intercept & gradient deviations
cf_ = mod.fixedEffects;
cr_ = mod.randomEffects;
if length(cf_) == 2, cf = cf_; else, cf(1) = cf_; end
if length(cr_) == nev
    if whichMod ~= 4, cr(:,1) = cr_; else, cr(:,2) = cr_; end
end
if length(cr_) == 2*nev, cr = reshape(cr_, [2 nev])'; end

rsd = mod.residuals;
stdDev = std(rsd);

% Calculate depth- and event-dependent expected values and store
% functions used to scale the data and model output
intercept = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
gradient = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
dat.scale_mu(d_ind) = intercept + gradient .* dat.Depth(d_ind);
dat.scale_sig(d_ind) = stdDev;

Data.(['scaleFun_' dv]) = @(mu,sig,x) (log(x) - mu) ./ sig; % scaling function

dat.scaled_Value(d_ind) = Data.(['scaleFun_' dv]) ...
    (dat.scale_mu(d_ind), dat.scale_sig(d_ind), dat.Value(d_ind)); % scaled data


logy = [logy 'chl_a'];


%% Store outputs in the Data struct
% Data.log_Depth = dat.log_Depth;
% Data.log_Value = dat.log_Value;
Data.scale_mu = dat.scale_mu;
Data.scale_sig = dat.scale_sig;
% Data.scale_intercept = dat.scale_intercept;
% Data.scale_gradient = dat.scale_gradient;
% Data.scale_stdDev = dat.scale_stdDev;
Data.scaled_Value = dat.scaled_Value;



%% Plots of standardised data

if ~isempty(v) && any(contains(v(1,:),'Scaled'))
    v_ind = find(contains(v(1,:),'Scaled'));
    for jj = 1:length(v_ind)
        vv = v(:,v_ind(jj));
        
        if endsWith(vv{1},logy)
            dv = logy{cellfun(@(x)~isempty(x),regexp(vv{1},logy))};
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
        if endsWith(vv{1},{'ScaledN'})
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
        xl_ = [-0.06 0.06] .* diff(xl) + xl;
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
        xl_ = [-0.06 0.06] .* diff(xl) + xl;
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
        xl_ = [-0.06 0.06] .* diff(xl) + xl;
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
        
        suptitle(['Compare raw and standardised ' dv ' measurements'])

    end
end



if ~isempty(v) && any(contains(v(1,:),'All'))
    % Plot all standardised data sources together
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

    
    
%     nbins = 9;
%     ind = strcmp(dat.Variable, 'PON');
%     ph1 = histogram(dat.scaled_Value(ind),'Normalization','pdf');
%     ph1.FaceColor = 'green'; ph1.FaceAlpha = 0.2;
%     ph1.NumBins = nbins;
%     hold on
%     xlabel('standardised data')
%     ylabel('sample density')
%     ind = strcmp(dat.Variable, 'POC');
%     ph2 = histogram(dat.scaled_Value(ind),'Normalization','pdf');
%     ph2.FaceColor = 'red'; ph2.FaceAlpha = 0.2;
%     ph2.NumBins = nbins;
%     ind = strcmp(dat.Variable, 'N');
%     ph3 = histogram(dat.scaled_Value(ind),'Normalization','pdf');
%     ph3.FaceColor = 'blue'; ph3.FaceAlpha = 0.2;
%     ph3.NumBins = nbins;
%     % resize
%     be = [ph1.BinEdges ph2.BinEdges ph3.BinEdges];
%     be = linspace(min(be),max(be),nbins+1);
%     ph1.BinEdges = be; ph2.BinEdges = be; ph3.BinEdges = be;
%     % overlay a standard normal distibution
%     nx = linspace(min(be), max(be), 100);
%     nd = (2*pi)^(-0.5) * exp(-0.5 .* nx .^ 2);
%     plot(nx,nd,'k')
%     % legend
%     ca = gca; xl = ca.XLim; yl = ca.YLim;
%     xl_ = xl(1) + [0.025 0.075] * diff(xl);
%     yl_ = yl(1)+0.85*diff(yl);
%     text(xl(1)+0.1*diff(xl),yl_,'standard normal')
%     plot(xl_, repmat(yl_, [1 2]), 'k')
%     hold off

    
    
end



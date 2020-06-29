function Data = standardiseFittingData(Data, varargin)

% Use linear mixed models to standardise fitting data.
% Returns linear-model coefficients used to rescale data as a function of
% depth and (if significant) sampling event.

dat = Data;
dat.waterMass = Data.waterMass(Data.Event);
dat = rmfield(dat, {'nSamples','nEvents','EventTraj'});
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
    
    % Store coefficients and standard deviation used to scale data.
    dat.scale_intercept(d_ind) = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
    dat.scale_gradient(d_ind) = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
    dat.scale_stdDev(d_ind) = repmat(stdDev, [sum(d_ind) 1]);
    % and the scaled data
    mu = dat.scale_intercept(d_ind) + dat.scale_gradient(d_ind) .* dat.Depth(d_ind);
    dat.scaled_Value(d_ind) = (dat.log_Value(d_ind) - mu) ./ dat.scale_stdDev(d_ind);
end


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

% Store coefficients and standard deviation used to scale data.
dat.scale_intercept(d_ind) = cf(1) + cr(dat.Event(d_ind)-min(uev)+1,1);
dat.scale_gradient(d_ind) = cf(2) + cr(dat.Event(d_ind)-min(uev)+1,2);
dat.scale_stdDev(d_ind) = repmat(stdDev, [sum(d_ind) 1]);
% and the scaled data
mu = dat.scale_intercept(d_ind) + dat.scale_gradient(d_ind) .* dat.log_Depth(d_ind);
dat.scaled_Value(d_ind) = (dat.Value(d_ind) - mu) ./ dat.scale_stdDev(d_ind);



%% Store outputs in the Data struct
Data.log_Depth = dat.log_Depth;
Data.log_Value = dat.log_Value;
Data.scale_intercept = dat.scale_intercept;
Data.scale_gradient = dat.scale_gradient;
Data.scale_stdDev = dat.scale_stdDev;
Data.scaled_Value = dat.scale_intercept;



%% Plots of standardised data

if ~isempty(v) && any(contains(v(1,:),'Scaled'))
    v_ind = find(contains(v(1,:),'Scaled'));
    for jj = 1:length(v_ind)
        vv = v(:,v_ind(jj));
        
        if endsWith(vv{1},{'POC','PON'})
            dv = dPOM{cellfun(@(x)~isempty(x),regexp(vv{1},dPOM))};
            d_ind = strcmp(dat.Variable, dv);
            d = dat(d_ind,:);
            uev = unique(d.Event);
            nev = length(uev);
            % Linear model variables
            x = d.Depth;
            y = d.log_Value;
            % Plot labels
            x_lab = ['ln(' dv ')']; x_lab2 = ['standardised ' dv];
            y_lab = 'depth (m)'; y_lab2 = 'sampling event';
        end
        if endsWith(vv{1},{'ScaledN'})
            dv = 'N';
            d_ind = strcmp(dat.Variable, dv);
            d = dat(d_ind,:);
            uev = unique(d.Event);
            nev = length(uev);
            % Linear model variables
            x = d.log_Depth;
            y = d.Value;
            % Plot labels
            x_lab = [dv ' (mmol N m^{-3})']; x_lab2 = ['standardised ' dv];
            y_lab = 'ln(depth)'; y_lab2 = 'sampling event';
        end
        
        ys = d.scaled_Value; % scaled data

        % Compare raw data to standardised data
        figure
        subplot(2,2,1)
        scatter(y, x, 'k')
        hold on
        ca = gca; xl = ca.XLim; yl = ca.YLim;
        xl_ = [-0.06 0.06] .* diff(xl) + xl;
        yl_ = [-0.06 0.06] .* diff(yl) + yl;
        mu = mean(y);
        s = std(y);
        pgon = polyshape([mu-s mu+s mu+s mu-s mu-s], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
        pg = plot(pgon);
        pg.FaceColor = [1 1 1] ./ 255;
        pg.FaceAlpha = 0.2;
        plot(repmat(mu, [1 2]), yl, 'k')
        xlim(xl); ylim(yl);
        xlabel(x_lab); ylabel(y_lab);
        hold off
        
        subplot(2,2,2)
        G = cell(nev,1);
        for i = 1:nev
            G{i} = y(d.Event==uev(i));
        end
        grp = cell2mat(arrayfun(@(i){i*ones(numel(G{i}),1)},(1:numel(G))'));
        boxplot(vertcat(G{:}),grp,'orientation','horizontal')
        hold on
        xlabel(x_lab); ylabel(y_lab2);
        yticklabels(num2str(uev))
        mu = mean(y);
        ca = gca;
        plot(repmat(mu,[1 2]), ca.YLim,'k')
        hold off
        
        subplot(2,2,3)
        scatter(ys, x, 'k')
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
        xlabel(x_lab2); ylabel(y_lab);
        hold off
        
        subplot(2,2,4)
        G = cell(nev,1);
        for i = 1:nev
            G{i} = ys(d.Event==uev(i));
        end
        grp = cell2mat(arrayfun(@(i){i*ones(numel(G{i}),1)},(1:numel(G))'));
        boxplot(vertcat(G{:}),grp,'orientation','horizontal')
        hold on
        xlabel(x_lab2); ylabel(y_lab2);
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
    subplot(3,1,1)
    ind = strcmp(dat.Variable, 'PON');
    scatter(dat.scaled_Value(ind),dat.Depth(ind),'g')
    hold on
    ylabel('depth (m)')
    ind = strcmp(dat.Variable, 'POC');
    scatter(dat.scaled_Value(ind),dat.Depth(ind),'r')
    ind = strcmp(dat.Variable, 'N');
    scatter(dat.scaled_Value(ind),dat.Depth(ind),'b')
    ca = gca; xl = ca.XLim; yl = ca.YLim;
    % legend
    text(xl(1)+0.1*diff(xl),yl(1)+0.95*diff(yl),'PON')
    xl_ = xl(1) + [0.025 0.075] * diff(xl);
    yl_ = yl(1) + [0.975 0.925] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = 'green';
    text(xl(1)+0.1*diff(xl),yl(1)+0.85*diff(yl),'POC')
    yl_ = yl(1) + [0.875 0.825] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = 'red';
    text(xl(1)+0.1*diff(xl),yl(1)+0.75*diff(yl),'N')
    yl_ = yl(1) + [0.775 0.725] * diff(yl);
    pgon = polyshape([xl_(1) xl_(2) xl_(2) xl_(1) xl_(1)], [yl_(1) yl_(1) yl_(2) yl_(2) yl_(1)]);
    pg = plot(pgon); pg.FaceColor = 'blue';
    hold off
    
    subplot(3,1,2)
    ind = strcmp(dat.Variable, 'PON');
    scatter(dat.scaled_Value(ind),dat.Event(ind),'g')
    hold on
    ylabel('sampling event')
    ind = strcmp(dat.Variable, 'POC');
    scatter(dat.scaled_Value(ind),dat.Event(ind),'r')
    ind = strcmp(dat.Variable, 'N');
    scatter(dat.scaled_Value(ind),dat.Event(ind),'b')
    hold off
    
    subplot(3,1,3)
    nbins = 9;
    ind = strcmp(dat.Variable, 'PON');
    ph1 = histogram(dat.scaled_Value(ind),'Normalization','pdf');
    ph1.FaceColor = 'green'; ph1.FaceAlpha = 0.2;
    ph1.NumBins = nbins;
    hold on
    xlabel('standardised data')
    ylabel('sample density')
    ind = strcmp(dat.Variable, 'POC');
    ph2 = histogram(dat.scaled_Value(ind),'Normalization','pdf');
    ph2.FaceColor = 'red'; ph2.FaceAlpha = 0.2;
    ph2.NumBins = nbins;
    ind = strcmp(dat.Variable, 'N');
    ph3 = histogram(dat.scaled_Value(ind),'Normalization','pdf');
    ph3.FaceColor = 'blue'; ph3.FaceAlpha = 0.2;
    ph3.NumBins = nbins;
    % resize
    be = [ph1.BinEdges ph2.BinEdges ph3.BinEdges];
    be = linspace(min(be),max(be),nbins+1);
    ph1.BinEdges = be; ph2.BinEdges = be; ph3.BinEdges = be;
    % overlay a standard normal distibution
    nx = linspace(min(be), max(be), 100);
    nd = (2*pi)^(-0.5) * exp(-0.5 .* nx .^ 2);
    plot(nx,nd,'k')
    % legend
    ca = gca; xl = ca.XLim; yl = ca.YLim;
    xl_ = xl(1) + [0.025 0.075] * diff(xl);
    yl_ = yl(1)+0.85*diff(yl);
    text(xl(1)+0.1*diff(xl),yl_,'standard normal')
    plot(xl_, repmat(yl_, [1 2]), 'k')
    hold off
    
end



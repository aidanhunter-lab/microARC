function plt = plot_fitToData_depthBoxplot(v, Data, modData, logPlot, varargin)
% Plot to compare model output to data used to tune the parameters.

% plt = figure;

% Default colours
coldat = [0 0 0];
colmod = [0 1 0];

%%
% Different selection of plots for scalar data than for size spectra data
switch v
    case unique(Data.scalar.Variable)
        % Dimensions
%         plt.Units = 'inches';
%         plt.Position = [0 0 16 12];
        
        % Plot labels
        switch v
            case 'N'
                Var = 'DIN';
                xlab = 'DIN (mmol N m^{-3})';
                xlab_log = 'DIN (log_{10} mmol N m^{-3})';
                xlab_std = 'standardised DIN';
            case 'PON'
                Var = 'PON';
                xlab = 'PON (mmol N m^{-3})';
                xlab_log = 'PON (log_{10} mmol N m^{-3})';
                xlab_std = 'standardised PON';
            case 'POC'
                Var = 'POC';
                xlab = 'POC (mmol C m^{-3})';
                xlab_log = 'POC (log_{10} mmol C m^{-3})';
                xlab_std = 'standardised POC';
            case 'chl_a'
                Var = 'chl-a';
                xlab = 'chl-a (mg chl-a m^{-3})';
                xlab_log = 'chl-a (log_{10} mg chl-a m^{-3})';
                xlab_std = 'standardised chl-a';
        end
        
        %% All standardised data -- ungrouped and sorted
        ind = strcmp(Data.scalar.Variable, v);
        n = sum(ind); % number of data points
        x = 1:n;
        o = Data.scalar.(['sortOrder_' v]);
        ydat = Data.scalar.scaled_Value(ind);      % observed data
        ymod = modData.scalar.scaled_Value(ind,:); % modelled equivalent
        ydat = ydat(o); % reorder
        ymod = ymod(o,:);
        
        event = Data.scalar.Event(ind);
        depth = Data.scalar.Depth(ind);
        value = Data.scalar.Value(ind); % observed data
        value_log = log10(value);
        value_scaled = Data.scalar.scaled_Value(ind);
        uevent = unique(event);
        udepth = unique(depth);
        nevent = length(uevent);
        ndepth = length(udepth);
        valueMod = modData.scalar.Value(ind,:); % modelled equivalent of data
        valueMod_log = log10(valueMod);
        valueMod_scaled = modData.scalar.scaled_Value(ind,:);
        clear allValues allValues_log allValues_scaled
        % Create matrix structure for boxplot
        for i = ndepth:-1:1
            di = depth == udepth(i);
            x = value(di);
            allValues(1:numel(x),2*ndepth-2*i+1) = x(:);
            x = valueMod(di,:);
            allValues(1:numel(x),2*ndepth-2*i+2) = x(:);
            x = value_log(di);
            allValues_log(1:numel(x),2*ndepth-2*i+1) = x(:);
            x = valueMod_log(di,:);
            allValues_log(1:numel(x),2*ndepth-2*i+2) = x(:);
            x = value_scaled(di);
            allValues_scaled(1:numel(x),2*ndepth-2*i+1) = x(:);
            x = valueMod_scaled(di,:);
            allValues_scaled(1:numel(x),2*ndepth-2*i+2) = x(:);
        end
        allValues(~isnan(allValues) & allValues == 0) = nan;
        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
        allValues_log(~isnan(allValues_log) & allValues_log == 0) = nan;
        
        switch logPlot
            case true
                bp = boxplot2(allValues_log, 'orientation', 'horizontal', 'barwidth', 0.7);
            case false
                bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7);
        end
        for i = 1:size(allValues,2)
            bp.out(i).Marker = '.';
        end
        for i = 1:size(allValues,2) / 2
            bp.out(2*i-1).MarkerFaceColor = coldat;
            bp.out(2*i-1).MarkerEdgeColor = coldat;
            bp.out(2*i).MarkerFaceColor = colmod;
            bp.out(2*i).MarkerEdgeColor = colmod;
            bp.box(2*i-1).Color = coldat;
            bp.box(2*i).Color = colmod;
            bp.med(2*i-1).Color = [0 0 0];
            bp.med(2*i).Color = [0 0 0];
        end
        gc = gca;
        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
        gc.YTick = mean(tt);
        gc.YTickLabel = num2str(flip(udepth));
        gc.TickLength = 0.5 * gc.TickLength;
        xl = gc.XLim; yl = gc.YLim;
        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
        hold on
        for i = 1:length(tt)
            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        end
        hold off
        switch logPlot
            case true
                xlabel(xlab_log)
            case false
                xlabel(xlab)
        end
        ylabel('depth (m)')
        % legend
        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
        
        
end

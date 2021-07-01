function plt = plot_fitToData(v, Data, modData, logPlot, varargin)
% Plot to compare model output to data used to tune the parameters.

extractVarargin(varargin)

if ~exist('errorBars', 'var')
    errorBars = false;
end

plt = figure;

% Default colours
coldat = [0 0 0];
colmod = [0 1 0];

%%
% Different selection of plots for scalar data than for size spectra data
switch v
    case unique(Data.scalar.Variable)
        % Dimensions
        plt.Units = 'inches';
        plt.Position = [0 0 16 12];
        
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
        subplot(2,3,1)        
        ind = strcmp(Data.scalar.Variable, v);
        n = sum(ind); % number of data points
        x = 1:n;
        o = Data.scalar.(['sortOrder_' v]);
        ydat = Data.scalar.scaled_Value(ind);      % observed data
        ymod = modData.scalar.scaled_Value(ind,:); % modelled equivalent
        ydat = ydat(o); % reorder
        ymod = ymod(o,:);
        % scatter-plot
        for i = 1:size(ymod, 2)
            scatter(x, ymod(:,i), 'Marker', '.', 'MarkerEdgeColor', colmod);
            if i == 1, hold on; end
        end
        scatter(x, ydat, 'MarkerEdgeColor', coldat)
        % asthetics
        gc = gca;
        gc.XLim = [0.5, n+0.5];        
        % legend
        xl = gc.XLim;
        yl = gc.YLim;
        
        xleg = [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
        yleg = [yl(1)+0.125*diff(yl), yl(1)+0.025*diff(yl)];
        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
            [1 1 1], 'EdgeColor', [1 1 1])
        text(xl(2)-0.15*diff(xl), yl(1)+0.1*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(1)+0.05*diff(yl), 'model')        
        scatter(xl(2)-0.17*diff(xl), yl(1)+0.1*diff(yl), 'MarkerEdgeColor', coldat)
        scatter(xl(2)-0.17*diff(xl), yl(1)+0.05*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)        
        ylabel(['sorted ' xlab_std])
        title([Var  ' data'])
        hold off

        %% All data on natural or log scale -- ungrouped, but sorted
        subplot(2,3,4)
        ydat = Data.scalar.Value(ind);      % observed data
        ymod = modData.scalar.Value(ind,:); % modelled equivalent        
        [~,o] = sort(ydat); % reorder
        ydat = ydat(o);
        ymod = ymod(o,:);
        if logPlot
            ydat = log10(ydat);
            ymod = log10(ymod);
        end
        % scatter plot
        for i = 1:size(ymod, 2)
            scatter(x, ymod(:,i), 'Marker', '.', 'MarkerEdgeColor', colmod);
            if i == 1, hold on; end
        end
        scatter(x, ydat, 'MarkerEdgeColor', coldat)
        % asthetics
        gc = gca;
        gc.XLim = [0.5, n+0.5];        
        % legend
        xl = gc.XLim;
        yl = gc.YLim;
        xleg = [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
        yleg = [yl(1)+0.125*diff(yl), yl(1)+0.025*diff(yl)];
        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
            [1 1 1], 'EdgeColor', [1 1 1])
        text(xl(2)-0.15*diff(xl), yl(1)+0.1*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(1)+0.05*diff(yl), 'model')
        scatter(xl(2)-0.17*diff(xl), yl(1)+0.1*diff(yl), 'MarkerEdgeColor', coldat)
        scatter(xl(2)-0.17*diff(xl), yl(1)+0.05*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)        
        switch logPlot
            case true
                ylabel(['sorted ' xlab_log])
            case false
                ylabel(['sorted ' xlab])
        end
        
        %% Standardised data grouped by depth -- boxplot
        subplot(2,3,2)
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
        % boxplot
        bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7);
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
        % asthetics
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
        xlabel(xlab_std)
        ylabel('depth (m)')
        title([Var  ' data grouped by depth'])
        gc.XLim = xl;
        % legend
        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
        
        %% data on natural or log scale grouped by depth -- boxplot                        subplot(2,3,5)
        subplot(2,3,5)
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
        
        %% standardised data grouped by event -- boxplot
        subplot(2,3,3)
        clear allValues allValues_log allValues_scaled
        for i = 1:nevent
            ei = event == uevent(i);
            x = value(ei);
            allValues(1:numel(x),2*i-1) = x(:);
            x = valueMod(ei,:);
            allValues(1:numel(x),2*i) = x(:);
            x = value_log(ei);
            allValues_log(1:numel(x),2*i-1) = x(:);
            x = valueMod_log(ei,:);
            allValues_log(1:numel(x),2*i) = x(:);
            x = value_scaled(ei);
            allValues_scaled(1:numel(x),2*i-1) = x(:);
            x = valueMod_scaled(ei,:);
            allValues_scaled(1:numel(x),2*i) = x(:);
        end
        allValues(~isnan(allValues) & allValues == 0) = nan;
        allValues_log(~isnan(allValues_log) & allValues_log == 0) = nan;
        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
        bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7);
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
        gc.YTickLabel = num2str(uevent);
        gc.TickLength = 0.5 * gc.TickLength;
        xl = gc.XLim; yl = gc.YLim;
        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
        hold on
        for i = 1:length(tt)
            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        end
        hold off
        xlabel(xlab_std)
        ylabel('sampling event')
        title([Var  ' data grouped by event'])
        % legend
        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
        
        %% data on natural or log scale grouped by event -- boxplot
        subplot(2,3,6)        
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
        gc.YTickLabel = num2str(uevent);
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
        ylabel('sampling event')        
        % legend
        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
    
    case unique(Data.size.dataBinned.Variable)
        
        if ~exist('waterOrigin', 'var')
            % By default do not group by water origin, just plot the fully
            % averaged size data
            waterOrigin = false;
        end
        
        switch waterOrigin
            
            case false
                
                plt.Units = 'inches';
                plt.Position = [0 0 16 6];
                
                switch v
                    case 'CellConc'
                        Var = 'cell concentration';
                        ylab = 'cell conc. (cells m^{-3})';
                        ylab_rel = 'relative cell conc.';
                        ylab_tot = 'total cell conc. (cells m^{-3})';
                    case 'BioVol'
                        Var = 'bio-volume';
                        ylab = 'bio-volume (m^3 m^{-3})';
                        ylab_rel = 'relative bio-volume';
                        ylab_tot = 'total bio-volume (m^3 m^{-3})';
                    case 'NConc'
                        Var = 'nitrogen concentration';
                        ylab = 'N conc. (mmol N m^{-3})';
                        ylab_rel = 'relative N conc.';
                        ylab_tot = 'total N conc. (mmol N m^{-3})';
                end
                
                %% Relative abundances
                subplot(1,3,1)
                ind = strcmp(Data.size.dataBinned.Variable, v);
                
                if ~isempty(varargin) && any(strcmp(varargin, 'trophicGroup'))
                    trophicGroup = varargin{find(strcmp(varargin, 'trophicGroup')) + 1};
                    ind = ind & strcmp(Data.size.dataBinned.trophicLevel, trophicGroup);
                end
                
                x = unique(Data.size.dataBinned.size(ind));
                ydat = Data.size.dataBinned.Value(ind);  % observed spectra
                ymod = modData.size.Value(ind,:);        % modelled equivalent
                ydat_tot = sum(ydat);
                ymod_tot = sum(ymod);
                ydat_rel = ydat ./ ydat_tot;
                ymod_rel = ymod ./ ymod_tot;
                switch logPlot
                    case 'loglog'
                        plotFun = @loglog;
                    case 'semilogx'
                        plotFun = @semilogx;
                end
                for i = 1:size(ymod_rel, 2)
                    plotFun(x, ymod_rel(:,i), 'Marker', '.', ...
                        'MarkerEdgeColor', colmod, 'Color', colmod);
                    if i == 1, hold on; end
                end
                plotFun(x, ydat_rel, '-o', 'Color', coldat)
                gc = gca;
                xlabel('ESD (\mum)')
                ylabel(ylab_rel)
                % legend
                xl = gc.XLim;
                yl = gc.YLim;
                switch logPlot
                    case 'loglog'
                        xl = log10(xl);
                        yl = log10(yl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = 10 .^ [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                    case 'semilogx'
                        xl = log10(xl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.05*diff(yl), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.1*diff(yl), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.1*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                end
                hold off
                
                %% total abundance
                subplot(1,3,2)
                scatter(repmat(1/3, [1 length(ymod_tot)]), ymod_tot, 'MarkerFaceColor', colmod, ...
                    'MarkerEdgeColor', colmod);
                hold on
                scatter(2/3, ydat_tot, 'MarkerEdgeColor', coldat);
                gc = gca;
                yl = gc.YLim;
                yl(1) = 0;
                yl(2) = 1.5 .* yl(2);
                gc.YLim = yl;
                gc.XLim = [0 1];
                gc.XTick = [];
                ylabel(ylab_tot)
                % legend
                xl = gc.XLim;
                yl = gc.YLim;
                xleg = [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                yleg = [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                    [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                    [1 1 1], 'EdgeColor', [1 1 1])
                text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                scatter(xl(2)-0.17*diff(xl), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', coldat)
                scatter(xl(2)-0.17*diff(xl), yl(2)-0.1*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                hold off
                
                %% Abundance spectra
                subplot(1,3,3)
                for i = 1:size(ymod, 2)
                    plotFun(x, ymod(:,i), 'Marker', '.', ...
                        'MarkerEdgeColor', colmod, 'Color', colmod);
                    if i == 1, hold on; end
                end
                plotFun(x, ydat, '-o', 'Color', coldat)
                gc = gca;
                xlabel('ESD (\mum)')
                ylabel(ylab)
                % legend
                xl = gc.XLim;
                yl = gc.YLim;
                switch logPlot
                    case 'loglog'
                        xl = log10(xl);
                        yl = log10(yl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = 10 .^ [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                    case 'semilogx'
                        xl = log10(xl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.05*diff(yl), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.1*diff(yl), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.1*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                end
                hold off
                
                if exist('trophicGroup', 'var')
                    sgtitle(['Model fit to ' trophicGroup ' ' Var ' data'])
                else
                    sgtitle(['Model fit to ' Var ' data'])
                end
                
                
            case {'Arctic','Atlantic'}
                
                
                plt.Units = 'inches';
                plt.Position = [0 0 16 6];
                
                switch v
                    case 'CellConc'
                        Var = 'cell concentration';
                        ylab = 'cell conc. (cells m^{-3})';
                        ylab_rel = 'relative cell conc.';
                        ylab_tot = 'total cell conc. (cells m^{-3})';
                    case 'BioVol'
                        Var = 'bio-volume';
                        ylab = 'bio-volume (m^3 m^{-3})';
                        ylab_rel = 'relative bio-volume';
                        ylab_tot = 'total bio-volume (m^3 m^{-3})';
                    case 'NConc'
                        Var = 'nitrogen concentration';
                        ylab = 'N conc. (mmol N m^{-3})';
                        ylab_rel = 'relative N conc.';
                        ylab_tot = 'total N conc. (mmol N m^{-3})';
                end
                
                %% Relative abundances
                subplot(1,3,1)

                dat = Data.size.dataBinned;
%                 dat = Data.sizeFull.dataBinned.groupedByOrigin;
                modDat = modData.size;
%                 modDat = modData.sizeFull;
                
                ind = strcmp(dat.Variable, v) & strcmp(dat.waterMass, waterOrigin);
                
                if exist('trophicGroup', 'var')
                    ind = ind & strcmp(dat.trophicLevel, eval('trophicGroup'));
                end
                
                x = unique(dat.size(ind));
                ydat = dat.Value(ind); % observation
                ymod = modDat.Value(ind,:); % modelled equivalents
%                 ymod = modDat.(['Value_' waterOrigin])(ind,:); % modelled equivalents

                ydat_tot = sum(ydat);
                ymod_tot = sum(ymod);
                ydat_rel = ydat ./ ydat_tot;
                ymod_rel = ymod ./ ymod_tot;
                switch logPlot
                    case 'loglog'
                        plotFun = @loglog;
                    case 'semilogx'
                        plotFun = @semilogx;
                end
                for i = 1:size(ymod_rel, 2)
                    plotFun(x, ymod_rel(:,i), 'Marker', '.', ...
                        'MarkerEdgeColor', colmod, 'Color', colmod);
                    if i == 1, hold on; end
                end
                plotFun(x, ydat_rel, '-o', 'Color', coldat)
                
                % GET BACK TO THIS... 
%                 switch errorBars, case true
%                     ydatSE = dat.ValueSE(ind); % observation
%                     
%                     
%                 end

                
                gc = gca;
                xlabel('ESD (\mum)')
                ylabel(ylab_rel)
                % legend
                xl = gc.XLim;
                yl = gc.YLim;
                switch logPlot
                    case 'loglog'
                        xl = log10(xl);
                        yl = log10(yl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = 10 .^ [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                    case 'semilogx'
                        xl = log10(xl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.05*diff(yl), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.1*diff(yl), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.1*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                end
                hold off
                
                %% total abundance
                subplot(1,3,2)
                scatter(repmat(1/3, [1 length(ymod_tot)]), ymod_tot, 'MarkerFaceColor', colmod, ...
                    'MarkerEdgeColor', colmod);
                hold on
                scatter(2/3, ydat_tot, 'MarkerEdgeColor', coldat);
                gc = gca;
                yl = gc.YLim;
                yl(1) = 0;
                yl(2) = 1.5 .* yl(2);
                gc.YLim = yl;
                gc.XLim = [0 1];
                gc.XTick = [];
                ylabel(ylab_tot)
                % legend
                xl = gc.XLim;
                yl = gc.YLim;
                xleg = [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                yleg = [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                    [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                    [1 1 1], 'EdgeColor', [1 1 1])
                text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                scatter(xl(2)-0.17*diff(xl), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', coldat)
                scatter(xl(2)-0.17*diff(xl), yl(2)-0.1*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                hold off
                
                %% Abundance spectra
                subplot(1,3,3)
                for i = 1:size(ymod, 2)
                    plotFun(x, ymod(:,i), 'Marker', '.', ...
                        'MarkerEdgeColor', colmod, 'Color', colmod);
                    if i == 1, hold on; end
                end
                plotFun(x, ydat, '-o', 'Color', coldat)
                gc = gca;
                xlabel('ESD (\mum)')
                ylabel(ylab)
                % legend
                xl = gc.XLim;
                yl = gc.YLim;
                switch logPlot
                    case 'loglog'
                        xl = log10(xl);
                        yl = log10(yl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = 10 .^ [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.05*diff(yl)), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), 10 .^ (yl(2)-0.1*diff(yl)), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                    case 'semilogx'
                        xl = log10(xl);
                        xleg = 10 .^ [xl(2)-0.2*diff(xl), xl(2)-0.0125*diff(xl)];
                        yleg = [yl(2)-0.125*diff(yl), yl(2)-0.025*diff(yl)];
                        fill([xleg(1), xleg(2), xleg(2), xleg(1), xleg(1)], ...
                            [yleg(1), yleg(1), yleg(2), yleg(2), yleg(1)], ...
                            [1 1 1], 'EdgeColor', [1 1 1])
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.05*diff(yl), 'data')
                        text(10 .^ (xl(2)-0.15*diff(xl)), yl(2)-0.1*diff(yl), 'model')
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', coldat)
                        scatter(10 .^ (xl(2)-0.17*diff(xl)), yl(2)-0.1*diff(yl), 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                end
                hold off
                
                if exist('trophicGroup', 'var')
                    sgtitle({['Model fit to ' eval('trophicGroup') ' ' Var ' data'], [waterOrigin ' waters']})
                else
                    sgtitle({['Model fit to ' Var ' data'], [waterOrigin ' waters']})
                end
                
                
        end
        
        
        
        
        
        
end

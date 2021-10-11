function plot_standardisedData2(Type, group, Data, varargin)

pause(0.1) % this pause prevents misnumbering the y-axis... I don't know why...

extractVarargin(varargin)

% Colours for data sampled from different water masses
if ~exist('colArctic', 'var')
    colArctic = [0 0 1];
end
if ~exist('colAtlantic', 'var')
    colAtlantic = [1 0 0];
end
if ~exist('colMix', 'var')
    colMix = [1 0 1];
end
waterMass = unique(Data.scalar.waterMass);
nWaterMass = length(waterMass);
col = cell(nWaterMass, 1);
for i = 1:length(waterMass)
    switch waterMass{i}
        case 'Atlantic'
            col{i,1} = colAtlantic;
        case 'Arctic'
            col{i,1} = colArctic;
        case 'Arctic/Atlantic'
            col{i,1} = colMix;
    end
end
colTable = table(waterMass, col);

if ~exist('pointAlpha', 'var')
    pointAlpha = 0.25;
end
if ~exist('lineAlpha', 'var')
    lineAlpha = 1;
end

if ~exist('includeLegend', 'var')
    includeLegend = false;
end
if ~exist('legendPosition', 'var')
    legendPosition = 'southwest';
end

if ~exist('covariate', 'var')
    covariate = [];
end
if isempty(covariate) && strcmp(Type, 'scalar')
    error('When Type="scalar" the "covariate" option must be set as "Depth" or "Event".')
end

if ~exist('densityCurve', 'var')
    densityCurve = false;
end

if ~exist('centrePlot', 'var')
    centrePlot = true;
end

if ~exist('Bandwidth', 'var')
    Bandwidth = 0.5; % smoothing parameter for kernel density approximations
end

if ~exist('groupAtlantic', 'var')
    groupAtlantic = false;
end
if ~exist('groupArctic', 'var')
    groupArctic = false;
end

if ~exist('logScale', 'var')
    logScale = false;
end

if ~exist('legendType', 'var')
    legendType = 'none';
end


switch Type
    case 'scalar'
        dat = Data.scalar;
        ind = ismember(dat.Variable, group);
        waterMass = dat.waterMass;
        event = unique(dat.Event);
        tt = table(event, waterMass);
        tt = innerjoin(tt, colTable);
        tt = sortrows(tt, 1); % sort by event to match dat
        event = dat.Event;
        tt_ = table(event);
        tt_ = innerjoin(tt_, tt);
        waterMass = unique(waterMass);
        
        % different plotting symbol for each member of group
        pltSymb = {'o','s','^','v','d','x'};
        if ischar(group)
            ngroup = 1;
        else
            ngroup = length(group);
        end
        pltSymb = pltSymb(1:ngroup);
        
        xg = cell(1,ngroup);
        
        if ngroup == 1, group = {group}; end

        for k = 1:ngroup
            for j = 1:length(waterMass)
                
                ind_ = ind & strcmp(tt_.waterMass, waterMass{j}) & ...
                    strcmp(dat.Variable, group{k}) ;
                x = dat.scaled_Value(ind_);
                y = dat.(covariate)(ind_);
                
                xg{:,k} = [xg{:,k}; x];
                
                switch covariate, case 'Depth'
                    y = -y;
                end
                col = unique(cell2mat(tt_.col(ind_)), 'rows');
                
                scatter(x, y, pltSymb{k}, ...
                    'MarkerEdgeColor', col, 'MarkerEdgeAlpha', pointAlpha);
                if j == 1 && k == 1, hold on, end
%                 if j == length(waterMass) && k == ngroup, hold off, end
            end
        end
        
        if ngroup == 1, group = group{1}; end
        
        xl = get(gca, 'XLim');
        yl = get(gca, 'YLim');
        
        switch legendType
            % Create legend from dummy plots
            case 'variables'
                pltLeg = cell(1, ngroup);
                for ii = 1:ngroup
                    pltLeg{ii} = scatter(10000, 10000, pltSymb{ii}, ...
                        'MarkerEdgeColor', [0 0 0]);
                end
                set(gca, {'XLim', 'YLim'}, {xl,yl})
            case 'distributions'
                pltLeg = cell(1, 2);
                lty = {'--','-'};
                if ~exist('lineWidth', 'var'), lineWidth = 0.5; end
                for ii = 1:2
                    pltLeg{ii} = plot([10000 10001], [10000 10001], lty{ii}, ...
                        'Color', [0 0 0 lineAlpha], 'LineWidth', lineWidth);
                end
                set(gca, {'XLim', 'YLim'}, {xl,yl})
        end
        
        switch logScale, case true
            set(gca, 'YScale', 'log')
        end
        
        if exist('YLim', 'var')
            set(gca, 'YLim', YLim)
        end
        
        if ischar(group)
            switch group
                case 'N'
                    xLab = 'DIN';
                case 'POC'
                    xLab = 'POC';
                case 'PON'
                    xLab = 'PON';
                case 'chl_a'
                    xLab = 'chl a';
            end
        else
            xLab = 'standarised data';
        end
        switch covariate
            case 'Depth'
                yLab = 'depth (m)';
            case 'Event'
                yLab = 'sample event';
        end
        
        % Adjust y-axis limits so points do not lay on axes
        adj = 0.06;
        miny = min(abs(dat.(covariate)(ind)));
        maxy = max(abs(dat.(covariate)(ind)));
        switch covariate, case 'Depth', miny = - miny; maxy = -maxy; end        
        
        
        switch covariate
            case 'Depth'
                YTickLabel = get(gca, 'YTickLabel');
                for j = 1:length(YTickLabel)
                    YTickLabel{j} = strrep(YTickLabel{j},'-','');
                end
                set(gca, 'YTickLabel', YTickLabel)
            case 'Event'
                uev = unique(dat.Event);
                nev = length(uev);
                if nev >= 15
                    YTick = 5:5:max(nev);
                    YTickLabel = arrayfun(@(z) {num2str(z)}, YTick);
                    YTickLabel = YTickLabel(:);
                    set(gca, {'YTick','YTickLabel'}, {YTick, YTickLabel})
                else
                    YTick = get(gca, 'YTick');
                    if YTick(1) == 0
                        YTickLabel = get(gca, 'YTickLabel');
                        set(gca, {'YTick','YTickLabel'}, {YTick(2:end), YTickLabel(2:end)})
                    end
                end
        end
        
        switch logScale, case true
            ll = get(gca, 'YTickLabels');
            for j = 1:length(ll)
                ll{j} = strrep(ll{j}, '{', '');
                ll{j} = strrep(ll{j}, '}', '');
                ll{j} = num2str(str2num(ll{j}));
            end
            set(gca, 'YTickLabels', ll)
        end
        
        yl = get(gca, 'YLim');
        
        yld = abs(diff(yl));
        if yl(1) == maxy
            yl(1) = yl(1) - adj * yld;
        end
        if yl(2) == miny
            yl(2) = yl(2) + adj * yld;
        end
        
        set(gca, 'YLim', yl)        

        xl = get(gca, 'XLim');
        switch centrePlot, case true
            xlm = max(abs(xl));
            xl = [-xlm, xlm];
            set(gca, 'XLim', xl)
        end
        
        xlabel(xLab)
        ylabel(yLab)
        
        
        switch densityCurve, case true
            hold on
            if ~exist('lineWidth', 'var'), lineWidth = 0.5; end
            x_ = linspace(xl(1), xl(2), 100);
%             lty
            for j = 1:ngroup
                yDist = fitdist(xg{j}, 'kernel', 'Bandwidth', Bandwidth);
                y_ = pdf(yDist, x_); % sample density of standardised data
                y_ = (y_ - min(y_)) ./ (max(y_) - min(y_));
                if exist('lineLim', 'var')
                    y_ = lineLim(1) + y_ .* diff(lineLim);
                else
                    y_ = yl(1) + y_ .* diff(yl);
                end
                plot(x_, y_, 'Color', [0 0 0 lineAlpha], 'LineWidth', lineWidth)
            end
            nd = (2*pi)^(-0.5) * exp(-0.5 .* x_ .^ 2); % standard normal density
            nd = (nd - min(nd)) ./ (max(nd) - min(nd));
            if exist('lineLim', 'var')
                nd = lineLim(1) + nd .* diff(lineLim);
            else
                nd = yl(1) + nd .* diff(yl);
            end
            plot(x_, nd, '--', 'Color', [0 0 0 lineAlpha], 'LineWidth', 2*lineWidth)
            hold off
        end

        
        
        switch includeLegend, case true
            if ~exist('legFontWeight', 'var')
                legFontWeight = 'bold';
            end
            switch legendType
                case 'variables'
                    if ~exist('legLab', 'var')
                        legLab = group;
                    end
                    if ~exist('legPosition', 'var')
                        legPosition = 'east';
                    end
                    if ~exist('legOrientation', 'var')
                        legOrientation = 'vertical';
                    end
                    if ~exist('legColumns', 'var')
                        legColumns = [];
                    end
                    pltLeg_ = nan(size(pltLeg));
                    for ii = 1:ngroup
                        pltLeg_(ii) = pltLeg{ii};
                    end
                    leg = legend(pltLeg_, legLab, 'Location', legPosition, ... 
                        'Orientation', legOrientation, 'NumColumns', legColumns);
                    if exist('legTitle', 'var')
                        leg.Title.String = eval('legTitle');
                        leg.Title.FontWeight = eval('legFontWeight');
                    end
                case 'distributions'
                    if ~exist('legLab', 'var')
                        legLab = [];
                    end
                    if ~exist('legPosition', 'var')
                        legPosition = 'east';
                    end
                    if ~exist('legOrientation', 'var')
                        legOrientation = 'vertical';
                    end
                    if ~exist('legColumns', 'var')
                        if strcmp(legOrientation, 'vertical')
                            legColumns = 1;
                        else
                            legColumns = 2;
                        end
                    end
                    pltLeg_ = nan(size(pltLeg));
                    for ii = 1:length(pltLeg_)
                        pltLeg_(ii) = pltLeg{ii};
                    end
                    leg = legend(pltLeg_, legLab, 'Location', legPosition, ...
                        'Orientation', legOrientation, 'NumColumns', legColumns);
                    if exist('legTitle', 'var')
                        leg.Title.String = eval('legTitle');
                        leg.Title.FontWeight = eval('legFontWeight');
                    end

                    
            end
        end
                
end


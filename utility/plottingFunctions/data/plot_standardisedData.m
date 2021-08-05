function plt = plot_standardisedData(Type, group, Data, varargin)

extractVarargin(varargin)

% if ~exist('Ylim', 'var')
%     Ylim = [];
% end

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

% if ~exist('groupAutotroph', 'var')
%     groupAutotroph = false;
% end
% if ~exist('groupHeterotroph', 'var')
%     groupHeterotroph = false;
% end
if ~exist('groupAtlantic', 'var')
    groupAtlantic = false;
end
if ~exist('groupArctic', 'var')
    groupArctic = false;
end


switch Type
    case 'scalar'
        dat = Data.scalar;
        ind = strcmp(dat.Variable, group);
        waterMass = dat.waterMass;
        event = unique(dat.Event);
        tt = table(event, waterMass);
        tt = innerjoin(tt, colTable);
        tt = sortrows(tt, 1); % sort by event to match dat
        event = dat.Event;
        tt_ = table(event);
        tt_ = innerjoin(tt_, tt);
        waterMass = unique(waterMass);
        for j = 1:length(waterMass)
            ind_ = ind & strcmp(tt_.waterMass, waterMass{j});
            x = dat.scaled_Value(ind_);
            y = dat.(covariate)(ind_);
            switch covariate, case 'Depth'
                y = -y;
            end
%             y = dat.Depth(ind_);
            col = unique(cell2mat(tt_.col(ind_)), 'rows');
            if j == 1
                plt = scatter(x, y, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
                    'MarkerFaceAlpha', pointAlpha, 'MarkerEdgeAlpha', pointAlpha);
                hold on
            else
                scatter(x, y, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
                    'MarkerFaceAlpha', pointAlpha, 'MarkerEdgeAlpha', pointAlpha)
            end
            if j == length(waterMass), hold off; end
        end
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
        xlabel(xLab)
        switch covariate
            case 'Depth'
                yLab = 'depth (m)';
            case 'Event'
                yLab = 'sample event';
        end
        ylabel(yLab)
        gc = gca;
%         if ~isempty(Ylim)
%             gc.YLim = flip(-Ylim);
%         end

        switch covariate, case 'Depth'
            for j = 1:length(gc.YTickLabel)
                gc.YTickLabel{j} = strrep(gc.YTickLabel{j},'-','');
            end
        end
        % Adjust y-axis limits so points do not lay on axes
        adj = 0.06;
        miny = min(abs(dat.(covariate)(ind)));
        maxy = max(abs(dat.(covariate)(ind)));
        switch covariate, case 'Depth', miny = - miny; maxy = -maxy; end        
        yl = gc.YLim;
        yld = abs(diff(yl));
        if yl(1) == maxy
            yl(1) = yl(1) - adj * yld;
        end
        if yl(2) == miny
            yl(2) = yl(2) + adj * yld;
        end
        gc.YLim = yl;
        
        switch centrePlot, case true
            xl = gc.XLim;
            xl = max(abs(xl));
            gc.XLim = [-xl, xl];
        end
        
        
        switch densityCurve, case true
            hold on
            yDist = fitdist(x, 'kernel', 'Bandwidth', Bandwidth);
            x_ = linspace(gc.XLim(1), gc.XLim(2), 100);
            y_ = pdf(yDist, x_); % sample density of standardised data
            y_ = (y_ - min(y_)) ./ (max(y_) - min(y_));
            y_ = gc.YLim(1) + y_ .* diff(gc.YLim);
            plot(x_, y_, 'Color', [0 0 0])
            nd = (2*pi)^(-0.5) * exp(-0.5 .* x_ .^ 2); % standard normal density
            nd = (nd - min(nd)) ./ (max(nd) - min(nd));
            nd = gc.YLim(1) + nd .* diff(gc.YLim);
            plot(x_, nd, '--', 'Color', [0 0 0])
            hold off
        end
        
        
        
        switch includeLegend, case true
            legend(waterMass, 'Location', legendPosition)
        end
        
        
        
        
end


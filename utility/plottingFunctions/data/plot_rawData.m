function plt = plot_rawData(Type, group, Data, varargin)

extractVarargin(varargin)

if ~exist('Ylim', 'var')
    Ylim = [];
end

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
if ~exist('legendTextSize', 'var')
    legendTextSize = 10;
end
if ~exist('legendTitle', 'var')
    legendTitle = [];
end
if ~exist('legendTitleSize', 'var')
    legendTitleSize = 10;
end
if ~exist('legFontWeight', 'var')
    legFontWeight = 'bold';
end

if ~exist('axesTextSize', 'var')
    axesTextSize = 10;
end

if ~exist('groupAutotroph', 'var')
    groupAutotroph = false;
end
if ~exist('groupHeterotroph', 'var')
    groupHeterotroph = false;
end
if ~exist('groupAtlantic', 'var')
    groupAtlantic = false;
end
if ~exist('groupArctic', 'var')
    groupArctic = false;
end

% if ~exist('waterOrigin', 'var')
%     waterOrigin = [];
% end


if ~exist('meanType', 'var')
    meanType = 'arithmetic';
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
            x = dat.Value(ind_);
            y = dat.Depth(ind_);
            col = unique(cell2mat(tt_.col(ind_)), 'rows');
            if j == 1
                plt = scatter(x, -y, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
                    'MarkerFaceAlpha', pointAlpha, 'MarkerEdgeAlpha', pointAlpha);
                hold on
            else
                scatter(x, -y, ...
                    'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
                    'MarkerFaceAlpha', pointAlpha, 'MarkerEdgeAlpha', pointAlpha)
            end
            if j == length(waterMass), hold off; end
        end
        switch group
            case 'N'
                xLab = 'DIN (mmol N m^{-3})';
            case 'POC'
                xLab = 'POC (mmol C m^{-3})';
            case 'PON'
                xLab = 'PON (mmol N m^{-3})';
            case 'chl_a'
                xLab = 'chl-a (mg m^{-3})';
        end
        
        set(gca, 'FontSize', axesTextSize)

        xlabel(xLab)
        ylabel('depth (m)')
        gc = gca;
        if ~isempty(Ylim)
            gc.YLim = flip(-Ylim);
        end
        for j = 1:length(gc.YTickLabel)
            gc.YTickLabel{j} = strrep(gc.YTickLabel{j},'-','');
        end
        % Adjust y-axis limits so points do not lay on axes
        adj = 0.06;
        minDepth = -min(abs(dat.Depth(ind)));
        maxDepth = -max(abs(dat.Depth(ind)));
        yl = gc.YLim;
        yld = abs(diff(yl));
        if yl(1) == maxDepth
            yl(1) = yl(1) - adj * yld;
        end
        if yl(2) == minDepth
            yl(2) = yl(2) + adj * yld;
        end
        gc.YLim = yl;
                        
        switch includeLegend, case true
            leg = legend(waterMass, 'Location', legendPosition);
            leg.Title.String = legendTitle;
            leg.Title.FontSize  = legendTitleSize;
            leg.Title.FontWeight = legFontWeight;
            leg.FontSize = legendTextSize;            
        end
        
        
        
    case 'sizeSpectra'
        % Make a single-panel plot of size-spectra curves for auto- &
        % heterotrophs sampled from Atlantic and Arctic waters. Different
        % colours for Atlantic/Arctic; different line types for trophic
        % level.
        
        dat = Data.sizeFull;
        x = unique(dat.ESD);
        nsize = length(x);
        
        waterMasses = {'Atlantic', 'Arctic'};
        trophicLevels = {'autotroph', 'heterotroph'};
        nWaterMasses = length(waterMasses);
        nTrophicLevels = length(trophicLevels);
        
        switch group
            case 'CellConc'
                y = dat.cellDensity;
                Ylab = {'cell concentration density', '(cells m^{-3} log_{10}(ESD / 1\mum)^{-1})'};
            case 'BioVol'
                y = dat.BioVolDensity;
                Ylab = {'bio-volume density', '(\mum^3m^{-3} log_{10}(ESD / 1\mum)^{-1})'};
        end
        
        waterMass = dat.waterMass;
        event = unique(dat.Event);
        tt = table(event, waterMass);
        tt = innerjoin(tt, colTable);
        tt = sortrows(tt, 1); % sort by event to match dat
        event = dat.Event;
        tt_ = table(event);
        tt_ = innerjoin(tt_, tt);
        
%         values = nan(nsize, nWaterMasses, nTrophicLevels);
%         for i = 1:length(waterMasses)
%             waterMass = waterMasses{i};
%             ind = strcmp(tt_.waterMass, waterMass);
%             for j = 1:length(trophicLevels)
%                 trophicLevel = trophicLevels{j};
%                 ind2 = ind & strcmp(dat.trophicLevel, trophicLevel);
%                 % Average over events and depths to produce single spectra
%                 % for each water mass and trophic level
%                 events = unique(dat.Event(ind2));
%                 nevent = length(events);
%                 value = zeros(nsize, 1);
%                 for k = 1:nevent
%                     event = events(k);
%                     ind3 = ind2 & dat.Event == event;
%                     depths = unique(dat.Depth(ind3));
%                     ndepth = length(depths);
%                     value_ = zeros(nsize, 1);
%                     for l = 1:ndepth
%                         depth = depths(l);
%                         ind4 = ind3 & dat.Depth == depth;
%                         value_ = value_ + y(ind4);
%                     end
%                     value_ = value_ ./ ndepth;
%                     value = value + value_;
%                 end
%                 value = value ./ nevent; % single size spectra for each water mass & trophic level (averaged over samples and depths)
%                 values(:,i,j) = value;
%             end
%         end
                
        values = nan(nsize, nWaterMasses, nTrophicLevels);
        switch meanType
            case 'arithmetic'
                valueBase = zeros(nsize, 1);
            case 'geometric'
                valueBase = ones(nsize, 1);
        end

        for i = 1:length(waterMasses)
            waterMass = waterMasses{i};
            ind = strcmp(tt_.waterMass, waterMass);
            for j = 1:length(trophicLevels)
                trophicLevel = trophicLevels{j};
                ind2 = ind & strcmp(dat.trophicLevel, trophicLevel);
                % Average over events and depths to produce single spectra
                % for each water mass and trophic level
                events = unique(dat.Event(ind2));
                nevent = length(events);
                value = valueBase;
                for k = 1:nevent
                    event = events(k);
                    ind3 = ind2 & dat.Event == event;
                    depths = unique(dat.Depth(ind3));
                    ndepth = length(depths);
                    value_ = valueBase;
                    for l = 1:ndepth
                        depth = depths(l);
                        ind4 = ind3 & dat.Depth == depth;
                        switch meanType
                            case 'arithmetic'
                                value_ = value_ + y(ind4);
                            case 'geometric'
                                value_ = value_ .* y(ind4);
                        end
                    end
                    switch meanType
                        case 'arithmetic'
                            value_ = value_ ./ ndepth;
                            value = value + value_;
                        case 'geometric'
                            value_ = value_ .^ (1 ./ ndepth);
                            value = value .* value_;
                    end
                end
                switch meanType
                    case 'arithmetic'
                        value = value ./ nevent; % single size spectra for each water mass & trophic level (averaged over samples and depths)
                    case 'geometric'
                        value = value .^ (1 ./ nevent);
                end
                values(:,i,j) = value;
            end
        end

        Arctic = strcmp(waterMasses, 'Arctic');
        Atlantic = strcmp(waterMasses, 'Atlantic');
        autotroph = strcmp(trophicLevels, 'autotroph');
        heterotroph = strcmp(trophicLevels, 'heterotroph');
        
        plotGroups = [groupAutotroph, groupHeterotroph, ... 
            groupAtlantic, groupArctic];
        
        if ~any(plotGroups)
            plt = loglog(x,values(:,Atlantic,autotroph), ...
                'Color', colAtlantic);
            hold on
            loglog(x,values(:,Atlantic,heterotroph), ...
                '--','Color', colAtlantic);
            loglog(x,values(:,Arctic,autotroph), ...
                'Color', colArctic);
            loglog(x,values(:,Arctic,heterotroph), ...
                '--', 'Color', colArctic);
            hold off
            switch includeLegend, case true
                leg = legend('Atlantic autotrophs', 'Atlantic heterotrophs', ...
                    'Arctic autotrophs', 'Arctic heterotrophs', ...
                    'Location', legendPosition);
                leg.Title.String = legendTitle;
                leg.Title.FontSize  = legendTitleSize;
                leg.FontSize = legendTextSize;
            end
        elseif sum(plotGroups) > 1
            error('Of the 4 plotGroups options ("groupAutotroph", "groupHeterotroph", "groupAtlantic", "groupArctic") only 1 at a time may be true.')
        elseif groupAutotroph
            plt = loglog(x,values(:,Atlantic,autotroph), ...
                'Color', colAtlantic);
            hold on
            loglog(x,values(:,Arctic,autotroph), ...
                'Color', colArctic);
            hold off
            switch includeLegend, case true
                leg = legend('Atlantic autotrophs', 'Arctic autotrophs', ... 
                    'Location', legendPosition);
                leg.Title.String = legendTitle;
                leg.Title.FontSize  = legendTitleSize;
                leg.Title.FontWeight  = legFontWeight;
                leg.FontSize = legendTextSize;
            end
        elseif groupHeterotroph
            plt = loglog(x,values(:,Atlantic,heterotroph), ...
                '--','Color', colAtlantic);
            hold on
            loglog(x,values(:,Arctic,heterotroph), ...
                '--', 'Color', colArctic);
            hold off
            switch includeLegend, case true
                leg = legend('Atlantic heterotrophs', 'Arctic heterotrophs', ...
                    'Location', legendPosition);
                leg.Title.String = legendTitle;
                leg.Title.FontSize  = legendTitleSize;
                leg.FontSize = legendTextSize;
            end
        elseif groupAtlantic
            plt = loglog(x,values(:,Atlantic,autotroph), ...
                'Color', colAtlantic);
            hold on
            loglog(x,values(:,Atlantic,heterotroph), ...
                '--','Color', colAtlantic);
            hold off
            switch includeLegend, case true
                legend('Atlantic autotrophs', 'Atlantic heterotrophs', ...
                    'Location', legendPosition)
            end
        elseif groupArctic
            plt = loglog(x,values(:,Arctic,autotroph), ...
                'Color', colArctic);
            hold on
            loglog(x,values(:,Arctic,heterotroph), ...
                '--', 'Color', colArctic);
            hold off
            switch includeLegend, case true
                leg = legend('Arctic autotrophs', 'Arctic heterotrophs', ...
                    'Location', legendPosition);
                leg.Title.String = legendTitle;
                leg.Title.FontSize  = legendTitleSize;
                leg.Title.FontWeight  = legFontWeight;
                leg.FontSize = legendTextSize;
            end
        end
        
        set(gca, 'FontSize', axesTextSize)

        gc = gca;
        if isempty(Ylim)
            switch group
                case 'CellConc'
                    gc.YLim(1) = 1 / diff(log10(x(1:2))); % truncate lower y-axis at 1 cell per m^3
%                     gc.YLim(1) = 1e-4;
                    gc.YLim(2) = 10 ^ ceil(log10(max(values(:))));
                case 'BioVol'
                    gc.YLim(1) = min(dat.cellVolume) / diff(log10(x(1:2))); % truncate lower y-axis at 1 cell per m^3
%                     gc.YLim(1) = 1e-12;
                    gc.YLim(2) = 10 ^ ceil(log10(max(values(:))));
            end
        elseif length(Ylim) == 2
            gc.YLim = Ylim;
        else
            error('Optional input "Ylim" must be an increasing 2-element vector.')
        end
        
        xlabel('ESD (\mum)')
        ylabel(Ylab)
           
    case 'sizeBinned'
        % Make a single-panel plot of the binned size data (used to fit the 
        % model) for auto- & heterotrophs sampled from Atlantic and Arctic 
        % waters. Different colours for Atlantic/Arctic; different line 
        % types for trophic level.
        
        dat = Data.sizeFull.dataBinned.groupedByOrigin;
        
%         waterMasses = unique(dat.waterMass);
%         trophicLevels = unique(dat.trophicLevel);        
        Atlantic = strcmp(dat.waterMass, 'Atlantic');
        Arctic = strcmp(dat.waterMass, 'Arctic');
        autotroph = strcmp(dat.trophicLevel, 'autotroph');
        heterotroph = strcmp(dat.trophicLevel, 'heterotroph');
        
        x_auto = unique(dat.size(autotroph));
        x_hetero = unique(dat.size(heterotroph));
        
        ind = strcmp(dat.Variable, group);

        switch group
            case 'CellConc'
                Ylab = {'cell concentration', '(cells m^{-3})'};
                scale = 1;
            case 'BioVol'
                Ylab = {'bio-volume', '(\mum^3m^{-3})'};
                scale = 1;
%                 Ylab = {'bio-volume', '(mm^3m^{-3})'};
%                 scale = 1e3 ^ 3; % convert m^3 / m^3 to mm^3 / m^3
        end

        y = dat.Value(Atlantic & autotroph & ind) * scale;        
        scatter(x_auto, y, 'MarkerEdgeColor', colAtlantic, ...
            'MarkerFaceColor', colAtlantic, 'MarkerFaceAlpha', pointAlpha)
        hold on
        y = dat.Value(Atlantic & heterotroph & ind) * scale;        
        scatter(x_hetero, y, 'MarkerEdgeColor', colAtlantic)
        y = dat.Value(Arctic & autotroph & ind) * scale;        
        scatter(x_auto, y, 'MarkerEdgeColor', colArctic, ...
            'MarkerFaceColor', colArctic, 'MarkerFaceAlpha', pointAlpha)
        y = dat.Value(Arctic & heterotroph & ind) * scale;
        scatter(x_hetero, y, 'MarkerEdgeColor', colArctic)
        hold off
        
        gc = gca;
        set(gc, 'XScale', 'log', 'YScale', 'log')
        
        if isempty(Ylim)
            switch group
                case 'CellConc'
                    % Truncate lower y-axis at 1 cell / m^3, which may
                    % exclude some negligible concentrations at small/large
                    % sizes
                    gc.YLim(1) = 1;
%                     gc.YLim(1) = 10 ^ floor(log10(min(dat.Value(ind)) * scale));
                    gc.YLim(2) = 10 ^ ceil(log10(max(dat.Value(ind)) * scale));
                case 'BioVol'
                    gc.YLim(1) = min(dat.size);
%                     gc.YLim(1) = 10 ^ floor(log10(min(dat.Value(ind)) * scale));                    
                    gc.YLim(2) = 10 ^ ceil(log10(max(dat.Value(ind)) * scale));
            end
        elseif length(Ylim) == 2
            gc.YLim = Ylim;
        else
            error('Optional input "Ylim" must be an increasing 2-element vector.')
        end
        
        xlabel('ESD (\mum)')
        ylabel(Ylab)
        
        gc.XTick = unique([x_auto; x_hetero]);
        gc.XTickLabel = round(unique([x_auto; x_hetero]), 2, 'significant');
        
        switch includeLegend, case true
            leg = legend('Atlantic autotrophs', 'Atlantic heterotrophs', ...
                'Arctic autotrophs', 'Arctic heterotrophs', ...
                'Location', legendPosition);
            leg.Title.String = legendTitle;
            leg.Title.FontSize  = legendTitleSize;
            leg.Title.FontWeight  = legFontWeight;
            leg.FontSize = legendTextSize;
        end
end


function plt = plot_rawData(Type, group, Data, varargin)

extractVarargin(varargin)

if ~exist('Ylim', 'var')
    Ylim = [];
end

if ~exist('waterOrigin', 'var')
    waterOrigin = [];
end

switch Type
    case 'scalar'
        dat = Data.scalar;
        ind = strcmp(dat.Variable, group);
        x = dat.Value(ind);
        y = dat.Depth(ind);
        plt = scatter(x, -y, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
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
        xlabel(xLab)
        ylabel('depth (m)')
        gc = gca;
        
        if ~isempty(Ylim)
            gc.YLim = flip(-Ylim);
        end
        
        for j = 1:length(gc.YTickLabel)
            gc.YTickLabel{j} = strrep(gc.YTickLabel{j},'-','');
        end
        
        
    case 'sizeSpectra'
        
        if isempty(waterOrigin)
            
            dat = Data.size;
            x = unique(dat.ESD);
            
            autotroph = strcmp(dat.trophicLevel, 'autotroph');
            heterotroph = strcmp(dat.trophicLevel, 'heterotroph');
            switch group
                case 'CellConc'
                    y = dat.cellDensity;
                    Ylab = {'cell concentration density', '(cells m^{-3} log_{10}(ESD / 1\mum)^{-1})'};
                case 'BioVol'
                    y = dat.BioVolDensity;
                    Ylab = {'bio-volume density', '(m^3m^{-3} log_{10}(ESD / 1\mum)^{-1})'};
            end
            
            y1 = y(autotroph);
            y2 = y(heterotroph);
            loglog(x, y1, 'Color', [0 1 0])
            hold on
            loglog(x, y2, 'Color', [1 0 0])
            hold off
            gc = gca;
            switch group
                case 'CellConc'
                    gc.YLim(1) = 1e-4;
                    gc.YLim(2) = 10 ^ ceil(log10(max([y1;y2])));
                case 'BioVol'
                    gc.YLim(1) = 1e-12;
                    gc.YLim(2) = 10 ^ ceil(log10(max([y1;y2])));
            end
            xlabel('ESD (\mum)')
            ylabel(Ylab)
            legend({'autotrophs','heterotrophs'}, 'Location', 'SouthWest')
            
        elseif ~any(strcmp(waterOrigin, {'Arctic', 'Atlantic'}))
            error('The optional waterOrigin argument must be either "Arctic" or "Atlantic".')
        else
            
            dat = Data.sizeFull.dataBinned.groupedByOrigin;
            sizes = unique(dat.size);
            
            autotroph = strcmp(dat.trophicLevel, 'autotroph');
            heterotroph = strcmp(dat.trophicLevel, 'heterotroph');
            origin =  strcmp(dat.waterMass, waterOrigin);
            Var =  strcmp(dat.Variable, group);
            ind0 = Var & origin;
            y = dat.Value;
            y1 = y(autotroph & ind0);
            y2 = y(heterotroph & ind0);
            
            if ~exist('colAutotroph', 'var'), colAutotroph = [0 1 0]; end
            if ~exist('colHeterotroph', 'var'), colHeterotroph = [1 0 0]; end
            
%             scatter(x, y1, 'MarkerFaceColor', colAutotroph, 'MarkerEdgeColor', colAutotroph)
            scatter(sizes, y1, 'MarkerFaceColor', colAutotroph, 'MarkerEdgeColor', colAutotroph)
            hold on
%             scatter(x, y2, 'MarkerFaceColor', colHeterotroph, 'MarkerEdgeColor', colHeterotroph)
            scatter(sizes, y2, 'MarkerFaceColor', colHeterotroph, 'MarkerEdgeColor', colHeterotroph)
            hold off
            
            set(gca, 'xscale', 'log')
            set(gca, 'yscale', 'log')
            
            gc = gca;
            
            if ~exist('matchScales', 'var'), matchScales = false; end
            
            switch matchScales, case true
%                 ylim = [min(y(Var)), max(y(Var))];
                gc.YLim = [min(y(Var)), max(y(Var))];
%                 set(gca, 'YLim', ylim)
            end
            
            if ~exist('drawLegend', 'var'), drawLegend = true; end            
            if drawLegend
                legend({'autotrophs','heterotrophs'}, 'Location', 'NorthEast')
            end
            
            switch group
                case 'CellConc'
                    Ylab = {'cell concentration', '(cells m^{-3})'};
                case 'BioVol'
                    Ylab = {'bio-volume density', '(m^3m^{-3})'};
            end
            
            if ~exist('drawXlabel', 'var'), drawXlabel = true; end
            if drawXlabel
                xlabel('ESD (\mum)')
            end

            if ~exist('drawYlabel', 'var'), drawYlabel = true; end
            if drawYlabel
                ylabel(Ylab)
            end
            
            gc = gca;
            gc.XTick = sizes;
            gc.XTickLabel = round(sizes, 2, 'significant');
            
            
            if ~exist('plotTitle', 'var')
                plotTitle = waterOrigin;
            end
            
            title(plotTitle)

        end
        
end






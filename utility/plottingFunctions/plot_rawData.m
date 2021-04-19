function plt = plot_rawData(Type, group, Data, varargin)

Ylim = [];
if ~isempty(varargin)
    if any(strcmp(varargin, 'Ylim'))
        Ylim = varargin{find(strcmp(varargin, 'Ylim')) + 1};
    end
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
        dat = Data.size;
        datb = Data.size.dataBinned;
        x = unique(dat.ESD);
        
%         xb1 = unique(datb.size(strcmp(datb.trophicLevel, 'autotroph')));
%         xb2 = unique(datb.size(strcmp(datb.trophicLevel, 'heterotroph')));
        
        autotroph = strcmp(dat.trophicLevel, 'autotroph');
        heterotroph = strcmp(dat.trophicLevel, 'heterotroph');
%         autotrophb = strcmp(datb.trophicLevel, 'autotroph');
%         heterotrophb = strcmp(datb.trophicLevel, 'heterotroph');
        switch group
            case 'CellConc'
                y = dat.cellDensity;
                Ylab = {'cell concentration density', '(cells m^{-3} log_{10}(ESD / 1\mum)^{-1})'};
            case 'BioVol'
                y = dat.BioVolDensity;
                Ylab = {'bio-volume density', '(m^3m^{-3} log_{10}(ESD / 1\mum)^{-1})'};
        end
%         ib = strcmp(datb.Variable, group);
        
        y1 = y(autotroph);
        y2 = y(heterotroph);
%         yb1 = datb.Value(ib & autotrophb);
%         yb2 = datb.Value(ib & heterotrophb);
        loglog(x, y1, 'Color', [0 1 0])
        hold on
%         scatter(xb1, yb1)
        loglog(x, y2, 'Color', [1 0 0])
%         scatter(xb2, yb2)
%         loglog(x, yb1, '--', 'Color', [0 1 0])        
%         loglog(x, yb2, '--', 'Color', [1 0 0])
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
end

function plt = plot_comparePMFs(yobs, ymod, varLabel, varargin)
% Compare distributions of binned size data and the modelled equivalents

extractVarargin(varargin)

if ~exist('itraj', 'var')
    itraj = 1; % ymod is calculated over a variety of trajectories, choose one to plot
end

if ~exist('trophicLevel', 'var')
   trophicLevel = []; 
end
if ~exist('waterMass', 'var')
   waterMass = []; 
end

if ~exist('plotTitle', 'var')
    plotTitle = true;
end
if ~exist('titleText', 'var')
    titleText = 'Data distribution vs modelled equivalent';
end
if ~isempty(waterMass) || ~isempty(trophicLevel)
    titleText = {titleText, [waterMass ' ' trophicLevel 's']};
end

if ~exist('plotLegend', 'var')
    plotLegend = true;
end
if ~exist('legendPosition', 'var')
    legendPosition = 'northwest';
end

if ~exist('colMod', 'var')
    colMod = [0, 1, 0];
end
if ~exist('colDat', 'var')
    colDat = [0, 0, 0];
end
if ~exist('markerSize', 'var')
    markerSize = 6;
end

if ~exist('barplot', 'var')
    barplot = false;
end
if ~exist('capwidth', 'var')
    capwidth = 0.075;
end



if ~exist('xscale', 'var')
    xscale = [];
end

if ~exist('yLabel', 'var')
    yLabel = [];
end
if isempty(yLabel)
    switch varLabel
        case 'BioVol'
            yLabel = {'probability mass', 'bio-volume at size'};
    end
end





%%

n = size(ymod, 1);
m = size(ymod, 2);

sizeClasses = 1:n;

yobsTot = sum(yobs);
ymodTot = sum(ymod);

cdfobs = cumsum(yobs) ./ yobsTot;
pdfobs = diff([0; cdfobs]);
cdfmod = cumsum(ymod) ./ ymodTot;
pdfmod = diff([zeros(1, m); cdfmod]);

if ~isempty(itraj) && length(itraj) == 1
    pdfmod = pdfmod(:,itraj);
    
    switch barplot
        
        case false
            plt = stem(sizeClasses, pdfobs, 'Color', colDat, 'XData', sizeClasses - 0.1);
            hold on
            stem(1:n, pdfmod, 'Color', colMod, 'XData', sizeClasses + 0.1);
            hold off
            plt.BaseLine.Parent.XTick = 1:n;
            if ~isempty(xscale) && length(xscale) == n
                plt.BaseLine.Parent.XTickLabel = xscale;
            end
            
        case true
            plt = bar(sizeClasses, [pdfobs pdfmod], 'grouped');
            plt(1).FaceColor = colDat;
            plt(2).FaceColor = colMod;
            plt(1).BaseLine.Parent.XTick = 1:n;           
            if ~isempty(xscale) && length(xscale) == n
                plt(1).BaseLine.Parent.XTickLabel = xscale;
            end
    end
    
    switch plotLegend, case true
        legend('Data', 'Model', 'Location', legendPosition)
    end
    
    xlabel('ESD (\mum)')
    ylabel(yLabel)
    
    switch plotTitle, case true
        title(titleText)
    end
    
elseif isempty(itraj)
    
    switch barplot
        
        case false
            plt = scatter(sizeClasses, pdfobs, 'MarkerEdgeColor', colDat, 'XData', sizeClasses - 0.1);
            hold on
            xp = repmat(sizeClasses', [1, m]);
            yp = pdfmod;
            scatter(xp(:), yp(:), 'MarkerEdgeColor', colMod, 'XData', xp(:) + 0.1);
            hold off
            plt.Parent.XTick = 1:n;
            if ~isempty(xscale) && length(xscale) == n
                plt.Parent.XTickLabel = xscale;
            end
            
        case true
            % compare data with the median modelled values, and show the
            % range as an additional interval
            pdfmodAv = median(pdfmod, 2);
            plt = bar(sizeClasses, [pdfobs, pdfmodAv], 'grouped');
            plt(1).FaceColor = colDat;
            plt(2).FaceColor = colMod;
            plt(1).BaseLine.Parent.XTick = 1:n;           
            if ~isempty(xscale) && length(xscale) == n
                plt(1).BaseLine.Parent.XTickLabel = xscale;
            end

            pdfmodMin = min(pdfmod, [], 2);
            pdfmodMax = max(pdfmod, [], 2);
            hold on
            for i = 1:n
                ix = plt(2).XEndPoints(i);
                il = pdfmodMin(i);
                im = pdfmodMax(i);
                plot([ix, ix], [il, im], 'Color', [0 0 0])
                plot(ix + [-1, 1] .* capwidth, [il, il], 'Color', [0 0 0])
                plot(ix + [-1, 1] .* capwidth, [im, im], 'Color', [0 0 0])
            end
            hold off
    end
    
    switch plotLegend, case true
        legend('Data', 'Model', 'Location', legendPosition)
    end
    xlabel('ESD (\mum)')
    ylabel(yLabel)
    switch plotTitle, case true
        title(titleText)
    end
end




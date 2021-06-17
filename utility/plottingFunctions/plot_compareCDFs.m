function plt = plot_compareCDFs(yobs, ymod, varLabel, varargin)
% Compare distributions of scalar (nutrient) data and the modelled
% equivalents

extractVarargin(varargin)

if ~exist('itraj', 'var')
    itraj = []; % plot all trajectory selections by default... produces messy plot
end
if ~exist('showTrueDataCDF', 'var')
    showTrueDataCDF = false;
end
if ~exist('plotTitle', 'var')
    plotTitle = true;
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
    colDat = [0.5, 0.5, 0.5];
end
if ~exist('colDatLine', 'var')
    colDatLine = [0, 0, 0];
end
if ~exist('markerSize', 'var')
    markerSize = 3;
end
if ~exist('titleText', 'var')
    titleText = 'Data distribution vs modelled equivalent';
end

n = size(ymod, 1);
m = size(ymod, 2);

yobs_sort = sort(yobs); % sort observations
ymod_sort = sort(ymod);

cdf = (1:n)' ./ n;

% remove any duplicate values (more likely in model output than
% in data), retaining the largest probability values
keep = ~[diff(yobs_sort) == 0; false];
cdfobs = cdf(keep);
yobs_sort_ = yobs_sort(keep);
% store cdfmod as cell-array because vector lengths may differ
% after removing duplicates
keep = ~[diff(ymod_sort) == 0; false(1, m)];
cdfmod = cell(1, m);
ymod_sort_ = cell(1, m);
for ij = 1:m
    cdfmod{:,ij} = cdf(keep(:,ij));
    ymod_sort_{:,ij} = ymod_sort(keep(:,ij),ij);
end

% find distances between each modelled CDF value and the
% equivalent data-CDF value
% interpolate model CDFs over data query points
cdfmodi = cell(1,m);
for ij = 1:m
    cdfmodi{ij} = interp1(ymod_sort_{ij}, cdfmod{ij}, yobs_sort_, 'linear', 'extrap');
    cdfmodi{ij}(cdfmodi{ij} < 0) = 0;
    cdfmodi{ij}(cdfmodi{ij} > 1) = 1;
end

% find distances between data and model CDFs
cdfDist = cell(1,m);
% IQD = nan(1,m);
for ij = 1:m
    cdfDist{ij} = abs(cdfobs - cdfmodi{ij});
%     IQD(ij) = mean(cdfDist{ij});
end

plt = plot(yobs_sort_, cdfobs, 'Color', colDatLine, 'LineWidth', 2);
hold on
if ~isempty(itraj)
    plot(yobs_sort_, cdfmodi{itraj}, 'o', 'MarkerEdgeColor', colMod, 'MarkerFaceColor', colMod, 'MarkerSize', markerSize)
    plot(yobs_sort_, cdfmodi{itraj}, '-o', 'Color', colMod, 'MarkerFaceColor', colMod, 'MarkerSize', markerSize)
else
    for jj = 1:m
        plot(yobs_sort_, cdfmodi{jj}, 'o', 'MarkerEdgeColor', colMod, 'MarkerFaceColor', colMod, 'MarkerSize', markerSize)
        plot(yobs_sort_, cdfmodi{jj}, '-o', 'Color', colMod, 'MarkerFaceColor', colMod, 'MarkerSize', markerSize)
    end
end
hold off

switch plotLegend, case true
    legend('Data', 'Model', 'Location', legendPosition)
end

switch varLabel, case 'chl_a'
    varLabel = 'chl a';
end

xlabel(['standardised ' varLabel])
ylabel('CDF')

switch plotTitle, case true
    title(titleText)
end

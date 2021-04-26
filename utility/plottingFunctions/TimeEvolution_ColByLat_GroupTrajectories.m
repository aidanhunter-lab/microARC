function plt = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, Var, Type)

% Plot output from groups of trajectories.
% Depth dimension is omitted by plotting min/max (across depth) values.
% Find across-trajectory means and standard deviations, plot expected
% output variable against time, coloured by latitude...

nTraj = Forc.nTraj;
groups = unique(Forc.waterMass); % groups of similar trajectories
ngroups = length(groups);

for i = 1:ngroups
    % indices for each group stored as group-name
    assignin('caller', groups{i}, strcmp(Forc.waterMass, groups{i}))
end

t = yearday(Forc.t(:,1));
n = length(t);

for i = 1:ngroups
    lat = Forc.y(:,evalin('caller', groups{i})); % latitude
    assignin('caller', ['lat_' groups{i}], lat)
end


switch Var
    % extract output variable(s)
    case 'DIN'
        x = squeeze(out.N);
    case 'OM_N'
        x = squeeze(out.OM(:,:,FixedParams.OM_N_index,:,:));
    case 'OM_C'
        x = squeeze(out.OM(:,:,FixedParams.OM_C_index,:,:));
end

for i = 1:ngroups    
    % aggregate output into trajectory groups
    ind = evalin('caller', groups{i});    
    nd = ndims(x); % trajectories are indexed last
    if nd == 3, assignin('caller', ['x_' groups{i}] , x(:,:,ind)); end
    if nd == 4, assignin('caller', ['x_' groups{i}] , x(:,:,:,ind)); end
end

for i = 1:ngroups
    % min and max values across depth
    x = evalin('caller', ['x_' groups{i}]);
    switch Var
        case 'DIN'
            assignin('caller', ['xmin_' groups{i}],  squeeze(min(x)))
            assignin('caller', ['xmax_' groups{i}],  squeeze(max(x)))
        case {'OM_C', 'OM_N'}
            assignin('caller', ['xmin_' groups{i}],  squeeze(min(x, [], 2)))
            assignin('caller', ['xmax_' groups{i}],  squeeze(max(x, [], 2)))
    end
end

logScale = true;
for i = 1:ngroups
    % mean and standard deviations across all trajectories in each group --
    % use log scale to avoid negative values
    xmin = evalin('caller', ['xmin_' groups{i}]);
    xmax = evalin('caller', ['xmax_' groups{i}]);
    
    if logScale
        xmin = log(xmin);
        xmax = log(xmax);
    end
    nd = ndims(xmax);
    xminm = mean(xmin, nd); % mean across trajectories
    xmaxm = mean(xmax, nd);
    xmins = std(xmin, 1, nd); % std dev across trajectories
    xmaxs = std(xmax, 1, nd);

    assignin('caller', ['xmin_mean_' groups{i}], xminm)
    assignin('caller', ['xmax_mean_' groups{i}], xmaxm)
    assignin('caller', ['xmin_std_' groups{i}], xmins)
    assignin('caller', ['xmax_std_' groups{i}], xmaxs)
    
    lat = evalin('caller', ['lat_' groups{i}]);
    latm = mean(lat,2);
    lats = std(lat,1,2);    
    assignin('caller', ['lat_mean_' groups{i}], latm)
    assignin('caller', ['lat_std_' groups{i}], lats)
end

for i = 1:ngroups
    % reshape for plotting functions
    xmin_mean = evalin('caller', ['xmin_mean_' groups{i}]);
    xmax_mean = evalin('caller', ['xmax_mean_' groups{i}]);        
    xmin_std = evalin('caller', ['xmin_std_' groups{i}]);
    xmax_std = evalin('caller', ['xmax_std_' groups{i}]);    
    lat_mean = evalin('caller', ['lat_mean_' groups{i}]);
    lat_std = evalin('caller', ['lat_std_' groups{i}]);
    
    nd = ndims(xmax_mean);
    
    xmin_mean = reshape(xmin_mean, [nd n]);
    xmax_mean = reshape(xmax_mean, [nd n]);
    xmin_std = reshape(xmin_std, [nd n]);
    xmax_std = reshape(xmax_std, [nd n]);
    lat_mean = reshape(lat_mean, [1 n]);
    lat_std = reshape(lat_std, [1 n]);
    
    assignin('caller', ['xmin_mean_' groups{i}], xmin_mean)
    assignin('caller', ['xmax_mean_' groups{i}], xmax_mean)
    assignin('caller', ['xmin_std_' groups{i}], xmin_std)
    assignin('caller', ['xmax_std_' groups{i}], xmax_std)    
    assignin('caller', ['lat_mean_' groups{i}], lat_mean)
    assignin('caller', ['lat_std_' groups{i}], lat_std)
end

for i = 1:ngroups
    % latitude limits used to set colour bar scales
    lat = evalin('caller', ['lat_mean_' groups{i}]);
    l = [min(lat) max(lat)];
    assignin('caller', ['latLim_' groups{i}], l)
    if i ~= 1
        latLim(1) = min([latLim(1) l]);
        latLim(2) = max([latLim(2) l]);
    else
        latLim = l;
    end
end

t = reshape(t, [1 n]);
z = zeros(size(t));


% Multi-panel plot -- columns are trajectory groups

plt = figure;

nc = ngroups;
% nr = 1; % extend the plot to include more rows later... want to all grouped trajectories on a map-plot
nr = size(evalin('caller', ['xmax_mean_' groups{1}]), 1);

widthCompress = 10; % reduce subplot widths by widthCompress% to make space for colour-bar

% alpha = 0.35; % transparency
lwd = 2; % line width

cm = flip(plasma); % default colour map
% adust limits of colour map for each group
cm_ = [cm linspace(latLim(1), latLim(2), size(cm, 1))'];
for i = 1:ngroups
    l = evalin('caller', ['latLim_' groups{i}]);
    ind = (cm_(:,4) >= l(1)) & (cm_(:,4) <= l(2));
    assignin('caller', ['cm_' groups{i}], cm(ind,:))
end

yl = nan(nr,2); % y-axis limits -- same across rows

% Loop over variables
for j = 1:nr
    
    % Loop over trajectory groups -- different subplot per variable and group
    
    for i = 1:nc
        ij = (j - 1) * nc + i;
        splt = subplot(nr, nc, ij);
        
        xmax_mean = evalin('caller', ['xmax_mean_' groups{i}]);
        xmin_mean = evalin('caller', ['xmin_mean_' groups{i}]);
        xmax_std = evalin('caller', ['xmax_std_' groups{i}]);
        xmin_std = evalin('caller', ['xmin_std_' groups{i}]);
        lat_mean = evalin('caller', ['lat_mean_' groups{i}]);        
        xmax_mean = xmax_mean(j,:);
        xmin_mean = xmin_mean(j,:);
        xmax_std = xmax_std(j,:);
        xmin_std = xmin_std(j,:);
        
        xminm = [xmin_mean; xmin_mean];
        xmaxm = [xmax_mean; xmax_mean];
        xmin_lo = repmat(xmin_mean - xmin_std, [2 1]);
        xmin_hi = repmat(xmin_mean + xmin_std, [2 1]);
        xmax_lo = repmat(xmax_mean - xmax_std, [2 1]);
        xmax_hi = repmat(xmax_mean + xmax_std, [2 1]);
        
        switch logScale
            case true
                xminm = exp(xminm);
                xmaxm = exp(xmaxm);
                xmin_lo = exp(xmin_lo);
                xmin_hi = exp(xmin_hi);
                xmax_lo = exp(xmax_lo);
                xmax_hi = exp(xmax_hi);
        end
        
        switch Type
            case 'minmax'
                surface([t; t], xmaxm, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', lwd)
                hold on
                surface([t; t], xminm, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', lwd)
                surface([t; t], xmin_lo, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                surface([t; t], xmin_hi, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                surface([t; t], xmax_lo, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                surface([t; t], xmax_hi, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                hold off
            case 'max'
                surface([t; t], xmaxm, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', lwd)
                hold on
                surface([t; t], xmax_lo, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                surface([t; t], xmax_hi, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                hold off
            case 'min'
                surface([t; t], xminm, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', lwd)
                hold on
                surface([t; t], xmin_lo, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                surface([t; t], xmin_hi, [z; z], [lat_mean; lat_mean], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1, 'LineStyle', ':');
                hold off
        end
        
        % axis labels
        if j == 1, title([groups{i} ' origin']); end
        if j == nr, xlabel('day of year'); end
        if i == 1
            switch Var
                case 'DIN'
                    yunit = 'mmol N m^{-3}';
                    ylab = 'DIN';
                case 'OM_N'
                    yunit = 'mmol N m^{-3}';
                    if FixedParams.DOM_index(j), ylab = 'DON'; end
                    if FixedParams.POM_index(j), ylab = 'PON'; end
                case 'OM_C'
                    yunit = 'mmol C m^{-3}';
                    if FixedParams.DOM_index(j), ylab = 'DOC'; end
                    if FixedParams.POM_index(j), ylab = 'POC'; end
            end
            switch Type
                case 'minmax'
                    tag = 'max and min across depth';
                case 'max'
                    tag = 'maximum across depth';
                case 'min'
                    tag = 'minimum across depth';
            end
            ylabel({[ylab ' (' yunit ')'], tag})
        end
        
        splt.Colormap = evalin('caller', ['cm_' groups{i}]); % subplot colour-map
        
        gc = gca;
        
        assignin('caller', ['gc_' groups{i} '_row' num2str(j)], gc)
%         assignin('caller', ['gc_' groups{i}], gc)
        
        % store y-limits to rescale y-axis across rows
        if i == 1
            yl(j,:) = [0, max(xmax_hi(:))];
        else
            yl(j,2) = max([yl(j,2); xmax_hi(:)]);
        end
        
        % reduce subplot widths to make room for colour-bar
        pos = gc.Position;
        wr = widthCompress / 100 * pos(3); % width reduction
        pos(3) = pos(3) - wr;
        if i ~= 1
            pos(1) = pos(1) - (i-1) * wr; % shift bottom left corner
        end
        
        gc.Position = pos;
        
    end
    
end

% rescale y-axis across rows
for j = 1:nr
    for i = 1:ngroups
        gc = evalin('caller', ['gc_' groups{i} '_row' num2str(j)]);
            gc.YLim = yl(j,:);
    end
end

% colour bar
if nr == 1
    gc = evalin('caller', ['gc_' groups{end} '_row' num2str(nr)]);
    pos = gc.Position;
else
    gc0 = evalin('caller', ['gc_' groups{end} '_row1']);
    gc1 = evalin('caller', ['gc_' groups{end} '_row' num2str(nr)]);
    pos = gc1.Position;
    pos0 = gc0.Position;
    pos(2) = 0.5 * (pos(2) + pos0(2));
%     pos(4) = sqrt(nr) * pos(4);
end



% gc = evalin('caller', ['gc_' groups{end}]);
% pos = gc.Position;
cbSpacing = [0.2, 0.925];
cbBL = cbSpacing(1) + (1 - cbSpacing(1)) * (pos(1) + pos(3));
cbWidth = cbSpacing(2) - cbBL;
cbPos = [cbBL, pos(2), cbWidth, pos(4)];
cb = colorbar('Position', cbPos);
deg = char(176);
cb.Label.String = ['latitude (' deg 'N)'];

cb.Limits = latLim;
colormap(cb, cm)



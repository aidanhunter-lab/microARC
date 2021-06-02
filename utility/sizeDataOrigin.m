function Data = sizeDataOrigin(Data)
% Group size data measurements according to whether the trajectories
% originate from the Arctic or Atlantic

dat = Data.sizeFull;

events = unique(dat.Event);

Atl = strcmp(dat.waterMass, 'Atlantic');
Arc = strcmp(dat.waterMass, 'Arctic');

evAtl = events(Atl);
evArc = events(Arc);

dat_ = dat.dataBinned;
dat = rmfield(dat_, {'waterMass','AtlanticOrigin','ArcticOrigin'});
dat = struct2table(dat);

% number of size classes [autotrophs, heterotrophs]
% trophicLevels = {'autotroph','heterotroph'};
nsizes = [length(unique(dat.sizeClass(strcmp(dat.trophicLevel, 'autotroph')))), ...
    length(unique(dat.sizeClass(strcmp(dat.trophicLevel, 'heterotroph'))))];


% Atlantic
ind = ismember(dat.Event, evAtl);

if any(ind)
    
    d = dat(ind,:); % separate data from water orginating from Atlantic
    
    % For each trophic level and measurement variable, take averages over
    % sampling events and depths.
    vars = unique(d.Variable, 'stable');
    nvars = length(vars);
    nev = length(evAtl);
    
    % autotrophs
    dd = d(strcmp(d.trophicLevel, 'autotroph'),:);
    
    values = nan(nsizes(1), nvars, nev);
    depth_av = nan(nvars, nev);
    yearday_av = nan(nvars, nev);
    
    ii = 0;
    for j = 1:nvars
        ij = strcmp(dd.Variable, vars{j});
        dj = dd(ij,:);
        for k = 1:nev
            jk = dj.Event == evAtl(k);
            dk = dj(jk,:);
            yearday_av(j,k) = dk.Yearday(1);
            depths = unique(dk.Depth);
            for l = 1:length(depths)
                ii = ii +1;
                jkl = dk.Depth == depths(l);
                dl = dk(jkl,:);
                values(:,j,k,ii) = dl.Value;
                depth_av(j,k,ii) = depths(l);
            end
        end
    end
    
    values(values == 0) = nan;
    depth_av(depth_av == 0) = nan;
    
    % average over depths
    v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
        v_av = mean(v_av, ndims(v_av));
        depth_av = mean(depth_av, ndims(depth_av));
        yearday_av = round(mean(yearday_av, ndims(yearday_av)));
    end
    
    nrows = nsizes(1) * nvars;
    
    % Build data table
    dd_av = table();
    dd_av.Year = repmat(dd.Year(1), [nrows, 1]);
    x = repmat(reshape(yearday_av, [1, nvars]), [nsizes(1), 1]);
    dd_av.Yearday = x(:);
    x = repmat(reshape(depth_av, [1, nvars]), [nsizes(1), 1]);
    dd_av.Depth = x(:);
    dd_av.trophicLevel = repmat(dd.trophicLevel(1), [nrows, 1]);
    x = repmat(reshape(unique(dd.size), [1, nsizes(1)]), [nvars, 1]);
    dd_av.size = x(:);
    x = repmat(1:nsizes(1), [nvars, 1]);
    dd_av.sizeClass = x(:);
    dd_av.Variable = repmat(vars, [nsizes(1), 1]);
    x = v_av';
    dd_av.Value = x(:);
    dd_av.order = (1:height(dd_av))';
    incostfunction = unique(dd(:,{'Variable' 'inCostFunction'}));
    dd_av = innerjoin(dd_av, incostfunction);
    [~, o] = sort(dd_av.order);
    dd_av = dd_av(o,:);
    dd_av.order = [];
    
    Dat.Atlantic.autotroph = dd_av;
    
    
    % heterotrophs
    dd = d(strcmp(d.trophicLevel, 'heterotroph'),:);
    
    values = nan(nsizes(2), nvars, nev);
    depth_av = nan(nvars, nev);
    yearday_av = nan(nvars, nev);
    
    ii = 0;
    for j = 1:nvars
        ij = strcmp(dd.Variable, vars{j});
        dj = dd(ij,:);
        for k = 1:nev
            jk = dj.Event == evAtl(k);
            dk = dj(jk,:);
            yearday_av(j,k) = dk.Yearday(1);
            depths = unique(dk.Depth);
            for l = 1:length(depths)
                ii = ii +1;
                jkl = dk.Depth == depths(l);
                dl = dk(jkl,:);
                values(:,j,k,ii) = dl.Value;
                depth_av(j,k,ii) = depths(l);
            end
        end
    end
    
    values(values == 0) = nan;
    depth_av(depth_av == 0) = nan;
    
    % average over depths
    v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
        v_av = mean(v_av, ndims(v_av));
        depth_av = mean(depth_av, ndims(depth_av));
        yearday_av = round(mean(yearday_av, ndims(yearday_av)));
    end
    
    nrows = nsizes(2) * nvars;
    
    % Build data table
    dd_av = table();
    dd_av.Year = repmat(dd.Year(1), [nrows, 1]);
    x = repmat(reshape(yearday_av, [1, nvars]), [nsizes(2), 1]);
    dd_av.Yearday = x(:);
    x = repmat(reshape(depth_av, [1, nvars]), [nsizes(2), 1]);
    dd_av.Depth = x(:);
    dd_av.trophicLevel = repmat(dd.trophicLevel(1), [nrows, 1]);
    x = repmat(reshape(unique(dd.size), [1, nsizes(2)]), [nvars, 1]);
    dd_av.size = x(:);
    x = repmat(1:nsizes(2), [nvars, 1]);
    dd_av.sizeClass = x(:);
    dd_av.Variable = repmat(vars, [nsizes(2), 1]);
    x = v_av';
    dd_av.Value = x(:);
    dd_av.order = (1:height(dd_av))';
    incostfunction = unique(dd(:,{'Variable' 'inCostFunction'}));
    dd_av = innerjoin(dd_av, incostfunction);
    [~, o] = sort(dd_av.order);
    dd_av = dd_av(o,:);
    dd_av.order = [];
    
    Dat.Atlantic.heterotroph = dd_av;
    
end


% Arctic
ind = ismember(dat.Event, evArc);

if any(ind)
    
    d = dat(ind,:); % separate data from water orginating from Atlantic
    
    % For each trophic level and measurement variable, take averages over
    % sampling events and depths.
    vars = unique(d.Variable, 'stable');
    nvars = length(vars);
    nev = length(evArc);
    
    % autotrophs
    dd = d(strcmp(d.trophicLevel, 'autotroph'),:);
    
    values = nan(nsizes(1), nvars, nev);
    depth_av = nan(nvars, nev);
    yearday_av = nan(nvars, nev);
    
    ii = 0;
    for j = 1:nvars
        ij = strcmp(dd.Variable, vars{j});
        dj = dd(ij,:);
        for k = 1:nev
            jk = dj.Event == evArc(k);
            dk = dj(jk,:);
            yearday_av(j,k) = dk.Yearday(1);
            depths = unique(dk.Depth);
            for l = 1:length(depths)
                ii = ii +1;
                jkl = dk.Depth == depths(l);
                dl = dk(jkl,:);
                values(:,j,k,ii) = dl.Value;
                depth_av(j,k,ii) = depths(l);
            end
        end
    end
    
    values(values == 0) = nan;
    depth_av(depth_av == 0) = nan;
    
    % average over depths
    v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
        v_av = mean(v_av, ndims(v_av));
        depth_av = mean(depth_av, ndims(depth_av));
        yearday_av = round(mean(yearday_av, ndims(yearday_av)));
    end
    
    nrows = nsizes(1) * nvars;
    
    % Build data table
    dd_av = table();
    dd_av.Year = repmat(dd.Year(1), [nrows, 1]);
    x = repmat(reshape(yearday_av, [1, nvars]), [nsizes(1), 1]);
    dd_av.Yearday = x(:);
    x = repmat(reshape(depth_av, [1, nvars]), [nsizes(1), 1]);
    dd_av.Depth = x(:);
    dd_av.trophicLevel = repmat(dd.trophicLevel(1), [nrows, 1]);
    x = repmat(reshape(unique(dd.size), [1, nsizes(1)]), [nvars, 1]);
    dd_av.size = x(:);
    x = repmat(1:nsizes(1), [nvars, 1]);
    dd_av.sizeClass = x(:);
    dd_av.Variable = repmat(vars, [nsizes(1), 1]);
    x = v_av';
    dd_av.Value = x(:);
    dd_av.order = (1:height(dd_av))';
    incostfunction = unique(dd(:,{'Variable' 'inCostFunction'}));
    dd_av = innerjoin(dd_av, incostfunction);
    [~, o] = sort(dd_av.order);
    dd_av = dd_av(o,:);
    dd_av.order = [];
    
    Dat.Arctic.autotroph = dd_av;
    
    
    % heterotrophs
    dd = d(strcmp(d.trophicLevel, 'heterotroph'),:);
    
    values = nan(nsizes(2), nvars, nev);
    depth_av = nan(nvars, nev);
    yearday_av = nan(nvars, nev);
    
    ii = 0;
    for j = 1:nvars
        ij = strcmp(dd.Variable, vars{j});
        dj = dd(ij,:);
        for k = 1:nev
            jk = dj.Event == evArc(k);
            dk = dj(jk,:);
            yearday_av(j,k) = dk.Yearday(1);
            depths = unique(dk.Depth);
            for l = 1:length(depths)
                ii = ii +1;
                jkl = dk.Depth == depths(l);
                dl = dk(jkl,:);
                values(:,j,k,ii) = dl.Value;
                depth_av(j,k,ii) = depths(l);
            end
        end
    end
    
    values(values == 0) = nan;
    depth_av(depth_av == 0) = nan;
    
    % average over depths
    v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
        v_av = mean(v_av, ndims(v_av));
        depth_av = mean(depth_av, ndims(depth_av));
        yearday_av = round(mean(yearday_av, ndims(yearday_av)));
    end
    
    nrows = nsizes(2) * nvars;
    
    % Build data table
    dd_av = table();
    dd_av.Year = repmat(dd.Year(1), [nrows, 1]);
    x = repmat(reshape(yearday_av, [1, nvars]), [nsizes(2), 1]);
    dd_av.Yearday = x(:);
    x = repmat(reshape(depth_av, [1, nvars]), [nsizes(2), 1]);
    dd_av.Depth = x(:);
    dd_av.trophicLevel = repmat(dd.trophicLevel(1), [nrows, 1]);
    x = repmat(reshape(unique(dd.size), [1, nsizes(2)]), [nvars, 1]);
    dd_av.size = x(:);
    x = repmat(1:nsizes(2), [nvars, 1]);
    dd_av.sizeClass = x(:);
    dd_av.Variable = repmat(vars, [nsizes(2), 1]);
    x = v_av';
    dd_av.Value = x(:);
    dd_av.order = (1:height(dd_av))';
    incostfunction = unique(dd(:,{'Variable' 'inCostFunction'}));
    dd_av = innerjoin(dd_av, incostfunction);
    [~, o] = sort(dd_av.order);
    dd_av = dd_av(o,:);
    dd_av.order = [];
    
    Dat.Arctic.heterotroph = dd_av;
    
end


% recombine by including a new struct within Data.sizeFull.datBinned
arctic = isfield(Dat, 'Arctic');
atlantic = isfield(Dat, 'Atlantic');
if atlantic
    Dat.Atlantic.autotroph.waterMass = repmat({'Atlantic'}, [height(Dat.Atlantic.autotroph), 1]);
    Dat.Atlantic.heterotroph.waterMass = repmat({'Atlantic'}, [height(Dat.Atlantic.heterotroph), 1]);
end
if arctic
    Dat.Arctic.autotroph.waterMass = repmat({'Arctic'}, [height(Dat.Arctic.autotroph), 1]);
    Dat.Arctic.heterotroph.waterMass = repmat({'Arctic'}, [height(Dat.Arctic.heterotroph), 1]);
end
if atlantic && arctic
    datNew = [Dat.Atlantic.autotroph; Dat.Atlantic.heterotroph; ...
        Dat.Arctic.autotroph; Dat.Arctic.heterotroph];
elseif atlantic
    datNew = [Dat.Atlantic.autotroph; Dat.Atlantic.heterotroph];
elseif arctic
    datNew = [Dat.Arctic.autotroph; Dat.Arctic.heterotroph];
end

datNew = table2struct(datNew, 'ToScalar', true);
dat_.groupedByOrigin = datNew;
Data.sizeFull.dataBinned = dat_;




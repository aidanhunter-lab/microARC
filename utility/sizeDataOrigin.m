function Data = sizeDataOrigin(Data, varargin)
% Group size data measurements according to whether the trajectories
% originate from the Arctic or Atlantic

extractVarargin(varargin)

if ~exist('avFun', 'var')
    avFun = @geomean; % use geometric mean by default... this is less biased by outlying values
end

if ~exist('returnVariability', 'var')
    returnVariability = false; % extra code section at end of script -- not needed, or finished...
end


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
%     v_av = mean(values, ndims(values), 'omitnan');
    v_av = avFun(values, ndims(values), 'omitnan'); % average over depths -- use avFun for values
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan'); % use arithmetic mean for depths
    % average over events
    if nev > 1
%         v_av = mean(v_av, ndims(v_av));
        v_av = avFun(v_av, ndims(v_av));
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
%     v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    v_av = avFun(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
%         v_av = mean(v_av, ndims(v_av));
        v_av = avFun(v_av, ndims(v_av));
        depth_av = mean(depth_av, ndims(depth_av));
        yearday_av = round(mean(yearday_av, ndims(yearday_av)));
    end
    

    
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
%     v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    v_av = avFun(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
%         v_av = mean(v_av, ndims(v_av));
        v_av = avFun(v_av, ndims(v_av));
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
%     v_av = mean(values, ndims(values), 'omitnan'); % average over depths
    v_av = avFun(values, ndims(values), 'omitnan'); % average over depths
    depth_av = mean(depth_av, ndims(depth_av), 'omitnan');
    % average over events
    if nev > 1
%         v_av = mean(v_av, ndims(v_av));
        v_av = avFun(v_av, ndims(v_av));
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



%% Derive measurement variabilities within the 2 groupings -- Arctic & Atlantic
% This is going to be a bit awkward in terms of getting the variabilities
% into the binned data... the whole data set-up may benefit from
% re-working... i.e., leave the size data unaltered (not binned/integrated)
% and integrate the model outputs to match the data... probably a good idea
% but a total pain to reowrk the code! It would, however, make other things
% much easier, e.g., fitting to multiple discrete depths per event.

switch returnVariability, case true
    
    
    dat = Data.sizeFull;
    dat = rmfield(dat, {'dataBinned', 'obsInCostFunction', 'nSamples', ...
        'EventTraj', 'evTraj', 'waterMass','AtlanticOrigin','ArcticOrigin'});
    dat = struct2table(dat);
    
    datAtl = dat(ismember(dat.Event, evAtl),:);
    datArc = dat(ismember(dat.Event, evArc),:);
    
    n = length(unique(dat.ESD));
    
    % We have 4 groupings: autotrophs & heterotrophs in Arctic & Atlantic
    % waters. Store all unique size spectra for each group
    % cellDensityArc_autotroph = nan(n, 1);
    % cellDensityArc_heterotroph = nan(n, 1);
    % biovolDensityArc_autotroph = nan(n, 1);
    % biovolDensityArc_heterotroph = nan(n, 1);
    
    trophicLevels = unique(dat.trophicLevel);
    
    % We're interested in variability between sampling events within data
    % aggregated by water origin. Variability across sample depths is
    % disruptive to this process... average it out.
    for i = 1:length(evArc)
        ev = evArc(i);
        dati = datArc(datArc.Event == ev,:);
        for j = 1:length(trophicLevels)
            trophicLevel = trophicLevels(j);
            datj = dati(strcmp(dati.trophicLevel, trophicLevel),:);
            depths = unique(datj.Depth);
            ndepth = length(depths);
            if ndepth == 1, continue; end
            depthMean = mean(depths);
            cellDensity = nan(n,ndepth); BioVolDensity = nan(n,ndepth);
            CellConc = nan(n,ndepth); BioVol = nan(n,ndepth);
            for k = 1:ndepth
                depth = depths(k);
                datk = datj(datj.Depth == depth,:);
                cellDensity(:,k) = datk.cellDensity;
                BioVolDensity(:,k) = datk.BioVolDensity;
                CellConc(:,k) = datk.CellConc;
                BioVol(:,k) = datk.BioVol;
            end
            % Average over depth
            cellDensity = mean(cellDensity, 2); BioVolDensity = mean(BioVolDensity, 2);
            CellConc = mean(CellConc, 2); BioVol = mean(BioVol, 2);
            datk.cellDensity = cellDensity;
            datk.BioVolDensity = BioVolDensity;
            datk.CellConc = CellConc;
            datk.BioVol = BioVol;
            datk{:,'Depth'} = depthMean;
            datArc(datArc.Event == ev & strcmp(datArc.trophicLevel, trophicLevel),:) = [];
            datArc = [datArc; datk];
        end
    end
    % Repeat for Atlantic
    for i = 1:length(evAtl)
        ev = evAtl(i);
        dati = datAtl(datAtl.Event == ev,:);
        for j = 1:length(trophicLevels)
            trophicLevel = trophicLevels(j);
            datj = dati(strcmp(dati.trophicLevel, trophicLevel),:);
            depths = unique(datj.Depth);
            ndepth = length(depths);
            if ndepth == 1, continue; end
            depthMean = mean(depths);
            cellDensity = nan(n,ndepth); BioVolDensity = nan(n,ndepth);
            CellConc = nan(n,ndepth); BioVol = nan(n,ndepth);
            for k = 1:ndepth
                depth = depths(k);
                datk = datj(datj.Depth == depth,:);
                cellDensity(:,k) = datk.cellDensity;
                BioVolDensity(:,k) = datk.BioVolDensity;
                CellConc(:,k) = datk.CellConc;
                BioVol(:,k) = datk.BioVol;
            end
            % Average over depth
            cellDensity = mean(cellDensity, 2); BioVolDensity = mean(BioVolDensity, 2);
            CellConc = mean(CellConc, 2); BioVol = mean(BioVol, 2);
            datk.cellDensity = cellDensity;
            datk.BioVolDensity = BioVolDensity;
            datk.CellConc = CellConc;
            datk.BioVol = BioVol;
            datk{:,'Depth'} = depthMean;
            datAtl(datAtl.Event == ev & strcmp(datAtl.trophicLevel, trophicLevel),:) = [];
            datAtl = [datAtl; datk];
        end
    end
    % Recover row order of tables
    [~,o] = sort(datArc.rowIndex);
    datArc = datArc(o,:);
    [~,o] = sort(datAtl.rowIndex);
    datAtl = datAtl(o,:);
    
    % Now find variabilities across events...
    nev = length(evArc);
    for i = 1:length(trophicLevels)
        trophicLevel = trophicLevels(i);
        ind_i = strcmp(datArc.trophicLevel, trophicLevel);
        dati = datArc(ind_i,:);
        cellDensity = nan(n, nev); BioVolDensity = nan(n, nev);
        CellConc = nan(n, nev); BioVol = nan(n, nev);
        for j = 1:nev
            ev = evArc(j);
            datj = dati(dati.Event == ev,:);
            cellDensity(:,j) = datj.cellDensity; BioVolDensity(:,j) = datj.BioVolDensity;
            CellConc(:,j) = datj.CellConc; BioVol(:,j) = datj.BioVol;
        end
        cellDensitySD = std(cellDensity, [], 2);
        BioVolDensitySD = std(BioVolDensity, [], 2);
        CellConcSD = std(CellConc, [], 2);
        BioVolSD = std(BioVol, [], 2);
        cellDensitySE = cellDensitySD ./ sqrt(nev);
        BioVolDensitySE = BioVolDensitySD ./ sqrt(nev);
        CellConcSE = CellConcSD ./ sqrt(nev);
        BioVolSE = BioVolSD ./ sqrt(nev);
        datArc.cellDensitySD(ind_i) = repmat(cellDensitySD, [nev, 1]);
        datArc.cellDensitySE(ind_i) = repmat(cellDensitySE, [nev, 1]);
        datArc.BioVolDensitySD(ind_i) = repmat(BioVolDensitySD, [nev, 1]);
        datArc.BioVolDensitySE(ind_i) = repmat(BioVolDensitySE, [nev, 1]);
        datArc.CellConcSD(ind_i) = repmat(CellConcSD, [nev, 1]);
        datArc.CellConcSE(ind_i) = repmat(CellConcSE, [nev, 1]);
        datArc.BioVolSD(ind_i) = repmat(BioVolSD, [nev, 1]);
        datArc.BioVolSE(ind_i) = repmat(BioVolSE, [nev, 1]);
    end
    % Repeat for Atlantic
    nev = length(evAtl);
    for i = 1:length(trophicLevels)
        trophicLevel = trophicLevels(i);
        ind_i = strcmp(datAtl.trophicLevel, trophicLevel);
        dati = datAtl(ind_i,:);
        cellDensity = nan(n, nev); BioVolDensity = nan(n, nev);
        CellConc = nan(n, nev); BioVol = nan(n, nev);
        for j = 1:nev
            ev = evAtl(j);
            datj = dati(dati.Event == ev,:);
            cellDensity(:,j) = datj.cellDensity; BioVolDensity(:,j) = datj.BioVolDensity;
            CellConc(:,j) = datj.CellConc; BioVol(:,j) = datj.BioVol;
        end
        cellDensitySD = std(cellDensity, [], 2);
        BioVolDensitySD = std(BioVolDensity, [], 2);
        CellConcSD = std(CellConc, [], 2);
        BioVolSD = std(BioVol, [], 2);
        cellDensitySE = cellDensitySD ./ sqrt(nev);
        BioVolDensitySE = BioVolDensitySD ./ sqrt(nev);
        CellConcSE = CellConcSD ./ sqrt(nev);
        BioVolSE = BioVolSD ./ sqrt(nev);
        datAtl.cellDensitySD(ind_i) = repmat(cellDensitySD, [nev, 1]);
        datAtl.cellDensitySE(ind_i) = repmat(cellDensitySE, [nev, 1]);
        datAtl.BioVolDensitySD(ind_i) = repmat(BioVolDensitySD, [nev, 1]);
        datAtl.BioVolDensitySE(ind_i) = repmat(BioVolDensitySE, [nev, 1]);
        datAtl.CellConcSD(ind_i) = repmat(CellConcSD, [nev, 1]);
        datAtl.CellConcSE(ind_i) = repmat(CellConcSE, [nev, 1]);
        datAtl.BioVolSD(ind_i) = repmat(BioVolSD, [nev, 1]);
        datAtl.BioVolSE(ind_i) = repmat(BioVolSE, [nev, 1]);
    end
    
    
    
    % Include variabilities within the grouped size data
    dat = Data.sizeFull.dataBinned.groupedByOrigin;
    
    for i = 1:length(trophicLevels)
        trophicLevel = trophicLevels(i);
        ind_i = strcmp(dat.trophicLevel, trophicLevel);
        dat.ValueSD(strcmp(dat.Variable, 'CellConc') & strcmp(dat.waterMass, 'Arctic') & ind_i,:) ...
            = unique(datArc.CellConcSD(strcmp(datArc.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSD(strcmp(dat.Variable, 'BioVol') & strcmp(dat.waterMass, 'Arctic') & ind_i,:) ...
            = unique(datArc.BioVolSD(strcmp(datArc.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSD(strcmp(dat.Variable, 'CellConc') & strcmp(dat.waterMass, 'Atlantic') & ind_i,:) ...
            = unique(datAtl.CellConcSD(strcmp(datAtl.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSD(strcmp(dat.Variable, 'BioVol') & strcmp(dat.waterMass, 'Atlantic') & ind_i,:) ...
            = unique(datAtl.BioVolSD(strcmp(datAtl.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSE(strcmp(dat.Variable, 'CellConc') & strcmp(dat.waterMass, 'Arctic') & ind_i,:) ...
            = unique(datArc.CellConcSE(strcmp(datArc.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSE(strcmp(dat.Variable, 'BioVol') & strcmp(dat.waterMass, 'Arctic') & ind_i,:) ...
            = unique(datArc.BioVolSE(strcmp(datArc.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSE(strcmp(dat.Variable, 'CellConc') & strcmp(dat.waterMass, 'Atlantic') & ind_i,:) ...
            = unique(datAtl.CellConcSE(strcmp(datAtl.trophicLevel, trophicLevel),:), 'stable');
        dat.ValueSE(strcmp(dat.Variable, 'BioVol') & strcmp(dat.waterMass, 'Atlantic') & ind_i,:) ...
            = unique(datAtl.BioVolSE(strcmp(datAtl.trophicLevel, trophicLevel),:), 'stable');
    end
    
    Data.sizeFull.dataBinned.groupedByOrigin = dat;
    
end





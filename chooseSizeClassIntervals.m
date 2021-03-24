function [dat, datFull] = chooseSizeClassIntervals(sizeData, varargin)

dat = struct2table(sizeData);
dat.rowIndex = (1:height(dat))';

datFull = [];
if ~isempty(varargin) && any(strcmp(varargin, 'datFull'))
    datFull = varargin{find(strcmp(varargin, 'datFull'))+1};
    datFull = struct2table(datFull);
    datFull.rowIndex = (1:height(datFull))';    
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Omit data outside selected size range

ESDmin = min(dat.ESD);
ESDmax = max(dat.ESD);
if ~isempty(varargin)
    vESDmin = strcmp(varargin, 'ESDmin');
    vESDmax = strcmp(varargin, 'ESDmax');
    if any(vESDmin | vESDmax)
        if any(vESDmin)
            ESDmin = varargin{find(vESDmin)+1};
        end
        if any(vESDmax)
            ESDmax = varargin{find(vESDmax)+1};
        end
    end
end
if ESDmin < min(dat.ESD), ESDmin = min(dat.ESD); end
if ESDmax > max(dat.ESD), ESDmax = max(dat.ESD); end
lESDmin = log10(ESDmin);
lESDmax = log10(ESDmax);
inRange = ESDmin <= dat.ESD & dat.ESD <= ESDmax;
dat = dat(inRange,:);

if ~isempty(datFull)
    inRange = ESDmin <= datFull.ESD & datFull.ESD <= ESDmax;
    datFull = datFull(inRange,:);
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Create size class intervals

% Has the number of intervals been passed as an argument?
nsizes = [];
nsizesP = [];
nsizesZ = [];
if ~isempty(varargin)
    vnsizes = strcmp(varargin, 'nsizes');
    vnsizesP = strcmp(varargin, 'nsizesP');
    vnsizesZ = strcmp(varargin, 'nsizesZ');
    if any(vnsizes), nsizes = varargin{find(vnsizes)+1}; end
    if any(vnsizesP), nsizesP = varargin{find(vnsizesP)+1}; end
    if any(vnsizesZ), nsizesZ = varargin{find(vnsizesZ)+1}; end
end

% There are 2 cases: (1) number of intervals is unknown
%                    (2) number of size class intervals was passed as
%                        argument

if isempty(nsizes) && (isempty(nsizesP) && isempty(nsizesZ)), Case = '1'; else, Case = '2'; end

switch Case
    case '1'
        % Number of size class intervals to return has not been specified, so
        % return numerous data sets, each corresponding to a different choice
        % of size partitioning.
        % Use the aggregated phytoplankton data
        dat = dat(strcmp(dat.trophicLevel, 'autotroph'),:);
        
        nsizesMin = 2; % hard code default values of min/max number of size class intervals
        nsizesMax = 25;
        % Has the range of intervals been passed as an argument?
        if ~isempty(varargin)
            vnsizesMin = strcmp(varargin, 'nsizesMin');
            vnsizesMax = strcmp(varargin, 'nsizesMax');
            if any(vnsizesMin)
                nsizesMin = varargin{find(vnsizesMin)+1};
            end
            if any(vnsizesMax)
                nsizesMax = varargin{find(vnsizesMax)+1};
            end
        end
        
        scenarios = unique(dat.scenario);
        nscenarios = length(scenarios);
        cols = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]]; % plotting colours for different scenarios
        allESD = unique(dat.ESD);
        nESD = length(allESD);
        dat2 = table((1:nESD)', allESD); % store all unique sizes in separate table
        dat2.Properties.VariableNames = {'index','ESD'};
        dat2.log10ESD = log10(dat2.ESD);
        
        nInt = nsizesMin:nsizesMax; % filter data using different numbers, nInt, of size class intervals
        
        % Assign all measurements to size class intervals
        for i = 1:length(nInt)
            dat3 = dat2;
            nsizes = nInt(i); % number of size class intervals
            b = linspace(lESDmin, lESDmax, nsizes+1); % interval edge positions
            lm = 0.5 * (b(1:end-1) + b(2:end)); % interval midpoints
            m = 10 .^ lm; % midpoint on natural scale
            rm = round(m * 4) / 4; % round to 0.25
            for j = nsizes:-1:1
                ind = dat3.log10ESD <= b(j+1);
                dat3.sizeClass(ind) = j; % assign size classes
                dat3.size(ind) = rm(j); % and mean sizes
            end
            dat3(:,{'index','log10ESD'}) = [];
            dat3 = innerjoin(dat, dat3); % merge tables
            [~,o] = sort(dat3.rowIndex);
            dat3 = dat3(o,:); % recover original row order
            dat3.rowIndex = [];
            dat_list.(['nIntervals' num2str(nInt(i))]) = dat3; % store all data
        end
        
        % Integrate size spectra across each size class interval
        for i = 1:length(nInt)
            dat3 = dat_list.(['nIntervals' num2str(nInt(i))]);
            dat3.CellConc = nan(height(dat3), 1); % integrate cell densities to find total cell number within each size class
            dat3.BioVol = nan(height(dat3), 1); % integrate biovolume densities to find biovolume within each size class
            dat3.NConc = nan(height(dat3), 1); % integrate nitrogen densities to find nitrogen concentration within each size class
            for j = 1:nscenarios
                ind0 = strcmp(dat3.scenario, scenarios(j));
                for k = 1:nsizes
                    ind = ind0 & dat3.sizeClass == k;
                    n = sum(ind);
                    x = log10(dat3.ESD(ind));
                    xd = diff(x);
                    % cell density
                    y = dat3.cellDensity(ind);
                    ys = y(1:end-1) + y(2:end);
                    CellConc = sum(0.5 * xd .* ys);
                    dat3.CellConc(ind) = repmat(CellConc, [n, 1]);
                    % biovolume
                    y = dat3.BioVolDensity(ind);
                    ys = y(1:end-1) + y(2:end);
                    BioVol = sum(0.5 * xd .* ys);
                    dat3.BioVol(ind) = repmat(BioVol, [n, 1]);
                    % nitrogen concentration
                    y = dat3.Ndensity(ind);
                    ys = y(1:end-1) + y(2:end);
                    NConc = sum(0.5 * xd .* ys);
                    dat3.NConc(ind) = repmat(NConc, [n, 1]);
                end
            end
            dat_list.(['nIntervals' num2str(nInt(i))]) = dat3;
        end
        
        % Include reduced size spectra data set containing only the binned values
        for i = 1:length(nInt)
            dat3 = dat_list.(['nIntervals' num2str(nInt(i))]);
            dat_binned = unique(dat3(:,{'scenario', 'Year', 'YeardayFirst', ...
                'YeardayLast', 'season', 'regime', 'DepthMin', 'DepthMax', 'size', ...
                'sizeClass', 'CellConc', 'BioVol', 'NConc'}));
            % convert from short to long format
            dat_binned = stack(dat_binned, {'CellConc','BioVol','NConc'}, ...
                'IndexVariableName','Variable','NewDataVariableName','Value');
            dat_binned.Variable = cellstr(dat_binned.Variable);
            % store as structs
            dat3 = table2struct(dat3, 'ToScalar', true);
            dat_binned = table2struct(dat_binned, 'ToScalar', true);
            dat3.dataBinned = dat_binned;            
            dat_list.(['nIntervals' num2str(nInt(i))]) = dat3;
        end
        
        dat = dat_list; % output size data for all possible size class intervals
        
        % plot the intervals
        makePlot = true; % by default makePlot is true when the number of intervals nsizes was unspecified
        if ~isempty(varargin)
            vplot = strcmp(varargin, 'plotSizeClassIntervals');
            if any(vplot)
                makePlot = varargin{find(vplot)+1};
            end
        end
        if makePlot
            plotVars = {'cellDensity', 'BioVolDensity', 'Ndensity'};
            N = length(nInt); % number of subplots (one for each number of intervals)
            maxSubplots = 20; % max subplots per figure
            nfig = ceil(N / maxSubplots); % number of figures
            if nfig > 1
                nc = ceil(sqrt(maxSubplots)); % number of columns
                nr = maxSubplots / nc; % number of rows
            else
                nc = ceil(sqrt(N));
                nr = N / nc;
            end
            nsubplots = nan(nfig,1);
            for i = 1:nfig
                if i < nfig
                    nsubplots(i) = maxSubplots;
                else
                    nsubplots(i) = mod(N, maxSubplots);
                end
            end
            
            for ij = 1:length(plotVars)
                plotVar = plotVars{ij};
                for i = 1:nfig
                    figure
                    for j = 1:nsubplots(i)
                        subplot(nr,nc,j)
                        k = (i-1) * maxSubplots + j + nsizesMin - 1; % number of intervals
                        d = dat.(['nIntervals' num2str(k)]);
                        for sc = 1:nscenarios
                            ind = strcmp(d.scenario, scenarios{sc});
                            loglog(d.ESD(ind), d.(plotVars{ij})(ind), 'Color', cols(sc,:))
                            if sc == 1, hold on; end
                        end
                        gc = gca;
                        
                        if strcmp(plotVar, 'cellDensity')
                            gc.YLim(1) = 1e1;
                        end
                        if strcmp(plotVar, 'BioVolDensity')
                            gc.YLim(1) = 1e-9;
                        end
                        if strcmp(plotVar, 'Ndensity')
                            gc.YLim(1) = 1e-6;
                        end
                        gc.YLim(2) = 2 * max(d.(plotVars{ij}));
                        jr = ceil(j/nc); % current row
                        if jr == nr
                            xlabel('ESD (\mum)')
                        end
                        title([num2str(k) ' intervals'])
                        b = linspace(lESDmin, lESDmax, k + 1); % interval boundaries
                        for q = 1:length(b)
                            line(10 .^ [b(q) b(q)], gc.YLim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--')
                        end
                    end
                    if strcmp(plotVars{ij}, 'cellDensity')
                        sgtitle({'Cell concentration density spectra','partitioned into various resolutions'})
                    end
                    if strcmp(plotVars{ij}, 'BioVolDensity')
                        sgtitle({'Biovolume density spectra','partitioned into various resolutions'})
                    end
                    if strcmp(plotVars{ij}, 'Ndensity')
                        sgtitle({'Nitrogen concentration density spectra','partitioned into various resolutions'})
                    end
                end
            end
        end
        
        
    case '2'
        % Number of size class intervals to return has been specified as
        % function argument
        
        if isempty(nsizesP), nsizesP = nsizes; end
        if isempty(nsizesZ), nsizesZ = nsizes; end
        
        % Arrange the aggregated size spectra data
        scenarios = unique(dat.scenario);
        nscenarios = length(scenarios);
        trophicLevels = unique(dat.trophicLevel);
        ntrophicLevels = length(trophicLevels);
        allESD = unique(dat.ESD);
        nESD = length(allESD);
        dat2 = table((1:nESD)', allESD); % store all unique sizes in separate table
        dat2.Properties.VariableNames = {'index','ESD'};
        dat2.log10ESD = log10(dat2.ESD);
        
        % Assign all measurements to size class intervals
        for j = 1:ntrophicLevels
            if strcmp(trophicLevels{j}, 'autotroph'), nsizes = nsizesP; end
            if strcmp(trophicLevels{j}, 'heterotroph'), nsizes = nsizesZ; end
            
            datj = dat(strcmp(dat.trophicLevel, trophicLevels{j}),:);
            if j > 1
                datj.size = [];
                datj.sizeClass = [];
            end
            
            dat3 = dat2;
            if strcmp(trophicLevels{j}, 'autotroph')
                b = linspace(lESDmin, lESDmax, nsizes+1); % interval edge positions
            end
            if strcmp(trophicLevels{j}, 'heterotroph')
                % adjust edge positions of smallest grazer size
                b = [lESDmin, linspace(lESDmin + 2 .* (lESDmax - lESDmin) ./ nsizesP, lESDmax, nsizes)];
            end
            
            lm = 0.5 * (b(1:end-1) + b(2:end)); % interval midpoints
            m = 10 .^ lm; % midpoint on natural scale
            rm = round(m * 4) / 4; % round to 0.25
            for jk = nsizes:-1:1
                ind = dat3.log10ESD <= b(jk+1);
                dat3.sizeClass(ind) = jk; % assign size classes
                dat3.size(ind) = rm(jk); % and mean sizes
            end
            dat3(:,{'index','log10ESD'}) = [];
            dat3 = innerjoin(datj, dat3); % merge tables
            [~,o] = sort(dat3.rowIndex);
            dat3 = dat3(o,:); % recover original row order
            dat3.rowIndex = [];
            
            dat.sizeClass(strcmp(dat.trophicLevel, trophicLevels{j})) = dat3.sizeClass;
            dat.size(strcmp(dat.trophicLevel, trophicLevels{j})) = dat3.size;
        end

        
        % Integrate size spectra across each size class interval
        
        dat.CellConc = nan(height(dat), 1); % integrate cell densities to find total cell number within each size class
        dat.BioVol = nan(height(dat), 1); % integrate biovolume densities to find total bio-volume within each size class
        dat.NConc = nan(height(dat), 1); % integrate nitrogen densities to find nitrogen concentration within each size class
        for j = 1:nscenarios
            indj = strcmp(dat.scenario, scenarios(j));
            for jj = 1:ntrophicLevels
                indjj = indj & strcmp(dat.trophicLevel, trophicLevels(jj));
                for k = 1:nsizes
                    ind = indjj & dat.sizeClass == k;
                    n = sum(ind);
                    x = log10(dat.ESD(ind));
                    xd = diff(x);
                    % cell density (cells m^-3)
                    y = dat.cellDensity(ind);
                    ys = y(1:end-1) + y(2:end);
                    CellConc = sum(0.5 * xd .* ys);
                    dat.CellConc(ind) = repmat(CellConc, [n, 1]);
                    % biovolume (m^3 m^-3)
                    y = dat.BioVolDensity(ind);
                    ys = y(1:end-1) + y(2:end);
                    BioVol = sum(0.5 * xd .* ys);
                    dat.BioVol(ind) = repmat(BioVol, [n, 1]);
                    % nitrogen concentration (mmol N m^-3)
                    y = dat.Ndensity(ind);
                    ys = y(1:end-1) + y(2:end);
                    NConc = sum(0.5 * xd .* ys);
                    dat.NConc(ind) = repmat(NConc, [n, 1]);
                end
            end
        end
        
        % Include reduced size spectra data set containing only the binned values
        dat_binned = unique(dat(:,{'scenario', 'Year', 'YeardayFirst', ...
            'YeardayLast', 'season', 'regime', 'DepthMin', 'DepthMax', 'trophicLevel', ...
            'size', 'sizeClass', 'CellConc', 'BioVol', 'NConc'}));
        
        % convert from short to long format
        dat_binned = stack(dat_binned, {'CellConc','BioVol','NConc'}, ...
            'IndexVariableName','Variable','NewDataVariableName','Value');
        dat_binned.Variable = cellstr(dat_binned.Variable);
        
        % store as structs
        dat = table2struct(dat, 'ToScalar', true);
        dat.nSamples = length(dat.Year);
        
        dat_binned = table2struct(dat_binned, 'ToScalar', true);
        dat.dataBinned = dat_binned;

        % plot the intervals
        makePlot = true; % makePlot is true by default
        if ~isempty(varargin)
            vplot = strcmp(varargin, 'plotSizeClassIntervals');
            if any(vplot)
                makePlot = varargin{find(vplot)+1};
            end
        end
        if makePlot
            % Multipanel plot. Top row: cell concentration
            %                  Bottom row: nitrogen concentration
            %                  Columns: sampling scenarios
            figure
            plotVars0 = {'cellDensity', 'BioVolDensity', 'Ndensity'};
            plotVars = unique(dat.dataBinned.Variable, 'stable');
            nr = length(plotVars); % number of rows
            nc = nscenarios; % number of columns
            b = linspace(lESDmin, lESDmax, nsizes+1); % interval boundaries
            db = diff(b);
            Cols = [0 1 0; 1 0 0];
            
            for i = 1:nr
                pv0 = plotVars0{i};
                for j = 1:nc
                    k = (i-1) * nc + j;
                    subplot(nr, nc, k);
                    
                    for jj = 1:ntrophicLevels
                        ind0 = strcmp(dat.scenario, scenarios{j}) & ...
                            strcmp(dat.trophicLevel, trophicLevels{jj}); % indexes continuous spectra
                        loglog(dat.ESD(ind0), dat.(plotVars0{i})(ind0), '--', 'Color', Cols(jj,:)) % plot continuous spectra (densities)
                        if jj == 1
                            hold on
                        end
                    end                    
                    for jj = 1:ntrophicLevels
                        ind0 = strcmp(dat.scenario, scenarios{j}) & ...
                            strcmp(dat.trophicLevel, trophicLevels{jj}); % indexes continuous spectra
                        y = dat.(plotVars{i})(ind0) ./ db(dat.sizeClass(ind0)); % divide by interval width to get densities from absolute values
                        loglog(dat.ESD(ind0), y, 'Color', Cols(jj,:)) % plot binned spectra (densities)
                    end
                    
                    gc = gca;
                    
                    switch pv0
                        case 'cellDensity'
                            gc.YLim(1) = 1e1;
                            ylab = 'cell concentration density';
                            yunit = 'cells$\;$m$^{-3}\;\log_{10}($ESD$/1\mu m)^{-1}$';
                        case 'BioVolDensity'
                            gc.YLim(1) = 1e-9;
                            ylab = 'bio-volume density';
                            yunit = 'm$^3\;$m$^{-3}\;\log_{10}($ESD$/1\mu m)^{-1}$';
                        case 'Ndensity'
                            gc.YLim(1) = 1e-6;
                            ylab = 'nitrogen concentration density';
                            yunit = 'mmol$\,$N m$^{-3}\;\log_{10}($ESD$/1\mu m)^{-1}$';
                    end
                    
                    gc.YLim(2) = 2 * max(dat.(pv0));
                    for q = 1:length(b)
                        line(10 .^ [b(q) b(q)], gc.YLim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':')
                    end
                    
                    if i == nr, xlabel('ESD (\mum)'); end
                    if j == 1
                        ylabel({ylab, ['(' yunit ')']}, 'Interpreter', 'latex');
                    end
                    if i == 1
                        yr = unique(dat.Year(ind0));
                        title(['scenario ' scenarios{j} ': ' num2str(yr)])
                    end
                    
                    hold off
                    
                end
            end
            sgtitle({'Measured size spectra',['partitioned into ' num2str(nsizes) ' evenly spaced intervals']})
        end

        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % Arrange the full size spectra data
        
        cruises = unique(datFull.Cruise);
        ncruises = length(cruises);
        trophicLevels = unique(datFull.trophicLevel);
        ntrophicLevels = length(trophicLevels);
        allESD = unique(datFull.ESD);
        nESD = length(allESD);
        dat2 = table((1:nESD)', allESD); % store all unique sizes in separate table
        dat2.Properties.VariableNames = {'index','ESD'};
        dat2.log10ESD = log10(dat2.ESD);
        
        % Assign all measurements to size class intervals
        for j = 1:ntrophicLevels
            
            datj = datFull(strcmp(datFull.trophicLevel, trophicLevels{j}),:);
            if j > 1
                datj.size = [];
                datj.sizeClass = [];
            end
            dat3 = dat2;
            
            if strcmp(trophicLevels{j}, 'autotroph'), nsizes = nsizesP; end
            if strcmp(trophicLevels{j}, 'heterotroph'), nsizes = nsizesZ; end
            if strcmp(trophicLevels{j}, 'autotroph')
                b = linspace(lESDmin, lESDmax, nsizes+1); % interval edge positions
            end
            if strcmp(trophicLevels{j}, 'heterotroph')
                % adjust edge positions of smallest grazer size
                b = [lESDmin, linspace(lESDmin + 2 .* (lESDmax - lESDmin) ./ nsizesP, lESDmax, nsizes)];
            end
            
            lm = 0.5 * (b(1:end-1) + b(2:end)); % interval midpoints
            m = 10 .^ lm; % midpoint on natural scale
            rm = round(m * 4) / 4; % round to 0.25
            for jk = nsizes:-1:1
                ind = dat3.log10ESD <= b(jk+1);
                dat3.sizeClass(ind) = jk; % assign size classes
                dat3.size(ind) = rm(jk); % and mean sizes
            end
            
            dat3(:,{'index','log10ESD'}) = [];
            dat3 = innerjoin(datj, dat3); % merge tables
            [~,o] = sort(dat3.rowIndex);
            dat3 = dat3(o,:); % recover original row order
            dat3.rowIndex = [];
            
            datFull.sizeClass(strcmp(datFull.trophicLevel, trophicLevels{j})) = dat3.sizeClass;
            datFull.size(strcmp(datFull.trophicLevel, trophicLevels{j})) = dat3.size;
        end
        
        % Integrate size spectra across each size class interval
        
        datFull.CellConc = nan(height(datFull), 1); % integrate cell densities to find total cell number within each size class
        datFull.BioVol = nan(height(datFull), 1); % integrate biovolume densities to find total bio-volume within each size class
        
        Events = unique(datFull.Event);
        for j = 1:length(Events)
            indj = datFull.Event == Events(j);
            for jj = 1:ntrophicLevels
                indjj = indj & strcmp(datFull.trophicLevel, trophicLevels{jj});
                Depths = unique(datFull.Depth(indjj));
                nDepths = length(Depths);
                for jd = 1:nDepths
                    indjd = indjj & datFull.Depth == Depths(jd);
                    for k = 1:nsizes
                        ind = indjd & datFull.sizeClass == k;
                        n = sum(ind);
                        x = log10(datFull.ESD(ind));
                        xd = diff(x);
                        % cell density (cells m^-3)
                        y = datFull.cellDensity(ind);
                        ys = y(1:end-1) + y(2:end);
                        CellConc = sum(0.5 * xd .* ys);
                        datFull.CellConc(ind) = repmat(CellConc, [n, 1]);
                        % biovolume (m^3 m^-3)
                        y = datFull.BioVolDensity(ind);
                        ys = y(1:end-1) + y(2:end);
                        BioVol = sum(0.5 * xd .* ys);
                        datFull.BioVol(ind) = repmat(BioVol, [n, 1]);
                    end
                end
            end
        end

        
        % Include reduced size spectra data set containing only the binned values
        dat_binned = unique(datFull(:,{'Label', 'Cruise', 'Year', 'Yearday', 'Event', ...
            'regime', 'Depth', 'trophicLevel', ...
            'size', 'sizeClass', 'CellConc', 'BioVol'}), 'stable');
        
        % convert from short to long format
        dat_binned = stack(dat_binned, {'CellConc','BioVol'}, ...
            'IndexVariableName','Variable','NewDataVariableName','Value');
        dat_binned.Variable = cellstr(dat_binned.Variable);
        
        % store as structs
        datFull = table2struct(datFull, 'ToScalar', true);
        %         datFull.nSamples = length(datFull.Year);
        
        dat_binned = table2struct(dat_binned, 'ToScalar', true);
        datFull.dataBinned = dat_binned;

        
        
end

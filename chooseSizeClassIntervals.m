function dat = chooseSizeClassIntervals(sizeData, varargin)

dat = struct2table(sizeData);
dat.rowIndex = (1:height(dat))';

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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Create size class intervals

% Has the number of intervals been passed as an argument?
nsizes = [];
if ~isempty(varargin)    
    vnsizes = strcmp(varargin, 'nsizes');
    if any(vnsizes)
        nsizes = varargin{find(vnsizes)+1};
    end
end

% There are 2 cases: (1) number of intervals is unknown
%                    (2) number of size class intervals was passed as
%                        argument

if isempty(nsizes), Case = '1'; else, Case = '2'; end

switch Case
    case '1'
        % Number of size class intervals to return has not been specified, so
        % return numerous data sets, each corresponding to a different choice
        % of size partitioning.
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
        
%         wTot = lESDmax - lESDmin; % log10-scale width of full size domain
%         dx = mean(diff(dat2.log10ESD)); % measurements were evenly spaced a distance dx apart on log10 scale
        
        nInt = nsizesMin:nsizesMax; % filter data using different numbers, nInt, of size class intervals
        
        % Assign all measurements to size class intervals
        for i = 1:length(nInt)
            dat3 = dat2;
            nsizes = nInt(i); % number of size class intervals
%             w = wTot / nsizes; % log10-scale width of each interval
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
                'sizeClass', 'CellConc', 'NConc'}));
            
            % convert from short to long format
            dat_binned = stack(dat_binned, {'CellConc','NConc'}, ...
                'IndexVariableName','Variable','NewDataVariableName','Value');
            dat_binned.Variable = cellstr(dat_binned.Variable);
            
            % store as structs
            dat3 = table2struct(dat3, 'ToScalar', true);
            dat_binned = table2struct(dat_binned, 'ToScalar', true);
            dat3.dataBinned = dat_binned;
            
            %         % Include fields specifying which of these data are used in the cost
            %         % function
            %         sizeData.obsInCostFunction = {'NConc'};
            %         sizeData.dataBinned.inCostFunction = ismember(sizeData.dataBinned.Variable, ...
            %             sizeData.obsInCostFunction);
            
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
            plotVars = {'cellDensity', 'Ndensity'};
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
                    if strcmp(plotVars{ij}, 'Ndensity')
                        sgtitle({'Nitrogen concentration density spectra','partitioned into various resolutions'})
                    end
                end
            end
        end
        
        
    case '2'
        % Number of size class intervals to return has been specified as
        % function argument
        scenarios = unique(dat.scenario);
        nscenarios = length(scenarios);
        cols = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]]; % plotting colours for different scenarios
        allESD = unique(dat.ESD);
        nESD = length(allESD);
        dat2 = table((1:nESD)', allESD); % store all unique sizes in separate table
        dat2.Properties.VariableNames = {'index','ESD'};
        dat2.log10ESD = log10(dat2.ESD);
        
        wTot = lESDmax - lESDmin; % log10-scale width of full size domain
        dx = mean(diff(dat2.log10ESD)); % measurements were evenly spaced a distance dx apart on log10 scale
        
%         nInt = nsizes; % filter data using different numbers, nInt, of size class intervals
        
        % Assign all measurements to size class intervals
        dat3 = dat2;
        w = wTot / nsizes; % log10-scale width of each interval
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
        dat = dat3;
            
        % Integrate size spectra across each size class interval
        dat.CellConc = nan(height(dat), 1); % integrate cell densities to find total cell number within each size class
        dat.NConc = nan(height(dat), 1); % integrate nitrogen densities to find nitrogen concentration within each size class
        for j = 1:nscenarios
            ind0 = strcmp(dat.scenario, scenarios(j));
            for k = 1:nsizes
                ind = ind0 & dat.sizeClass == k;
                n = sum(ind);
                x = log10(dat.ESD(ind));
                xd = diff(x);
                % cell density
                y = dat.cellDensity(ind);
                ys = y(1:end-1) + y(2:end);
                CellConc = sum(0.5 * xd .* ys);
                dat.CellConc(ind) = repmat(CellConc, [n, 1]);
                % nitrogen concentration
                y = dat.Ndensity(ind);
                ys = y(1:end-1) + y(2:end);
                NConc = sum(0.5 * xd .* ys);
                dat.NConc(ind) = repmat(NConc, [n, 1]);
            end
        end
        
        % Include reduced size spectra data set containing only the binned values
        dat_binned = unique(dat(:,{'scenario', 'Year', 'YeardayFirst', ...
            'YeardayLast', 'season', 'regime', 'DepthMin', 'DepthMax', 'size', ...
            'sizeClass', 'CellConc', 'NConc'}));
        
        % convert from short to long format
        dat_binned = stack(dat_binned, {'CellConc','NConc'}, ...
            'IndexVariableName','Variable','NewDataVariableName','Value');
        dat_binned.Variable = cellstr(dat_binned.Variable);
        
        % store as structs
        dat = table2struct(dat, 'ToScalar', true);
        dat.nSamples = length(dat.Year);

        dat_binned = table2struct(dat_binned, 'ToScalar', true);
        dat.dataBinned = dat_binned;
        
        %         % Include fields specifying which of these data are used in the cost
        %         % function
        %         sizeData.obsInCostFunction = {'NConc'};
        %         sizeData.dataBinned.inCostFunction = ismember(sizeData.dataBinned.Variable, ...
        %             sizeData.obsInCostFunction);
        
        
        
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
            plotVars0 = {'cellDensity', 'Ndensity'};
            plotVars = unique(dat.dataBinned.Variable);
            nr = length(plotVars); % number of rows
            nc = nscenarios; % number of columns
            b = linspace(lESDmin, lESDmax, nsizes+1); % interval boundaries
            db = diff(b);
            
            for i = 1:nr
                for j = 1:nc
                    k = (i-1) * nc + j;
                    subplot(nr, nc, k);
                    ind0 = strcmp(dat.scenario, scenarios{j});
                    ind = strcmp(dat.dataBinned.scenario, scenarios{j}) & ...
                        strcmp(dat.dataBinned.Variable, plotVars{i});
                    loglog(dat.ESD(ind0), dat.(plotVars0{i})(ind0), 'Color', cols(j,:))
                    hold on
                    
                    y = dat.(plotVars{i})(ind0) ./ db(dat.sizeClass(ind0));
                    loglog(dat.ESD(ind0), y, 'Color', cols(j,:))
                    
                    gc = gca;
                    
                    if strcmp(plotVars0{i}, 'cellDensity'), gc.YLim(1) = 1e1; end
                    if strcmp(plotVars0{i}, 'Ndensity'), gc.YLim(1) = 1e-6; end
                    gc.YLim(2) = 2 * max(dat.(plotVars0{i}));
                    for q = 1:length(b)
                        line(10 .^ [b(q) b(q)], gc.YLim, 'Color', [0.5, 0.5, 0.5], 'LineStyle', ':')
                    end
                    
                    if i == nr, xlabel('ESD (\mum)'); end
                    if j == 1
                        if strcmp(plotVars{i}, 'NConc')
                            ylabel('Nitrogen density (mmol N m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$)', 'Interpreter', 'latex')
                        end
                        if strcmp(plotVars{i}, 'CellConc')
                            ylabel('cell conc. density (cells m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$)', 'Interpreter', 'latex')
                        end
                    end
                    if i == 1
                        yr = unique(dat.Year(ind0));
%                         se = unique(dat.season(ind0));
                        title(['scenario ' scenarios{j} ': ' num2str(yr)])
%                         title(['scenario ' scenarios{j} ': ' se{1} ' ' num2str(yr)])
                    end
                    
                    hold off
                    
                end
            end
            sgtitle({'Measured size spectra',['partitioned into ' num2str(nsizes) ' evenly spaced intervals']})
        end
        
end





















% % For each possible number of size class intervals, loop through all
% % possible configurations storing a measure of misfit to data, then store
% % the best fitting approximations. AIC could be used to choose between
% % different numbers of size class intervals...
% % Hmm, this may be overcomplicating things. If widths of the first and last
% % intervals are different from the rest, and midpoints are to be equally
% % spaced, then the true midpoints of the first and last intervals will
% % differ from what the model uses! JUST SET SIMPLE EQUALLY SPACED INTERVALS
% % ALL OF EQUALL WIDTH, THEN CHOOSE NUMBER OF INTERVALS.
% 
% 
% % I SHOULD SIMPLYFY THIS FURTHER!!! DO NOT AUTOMATICALLY CHOOSE THE NUMBER
% % OF SIZE INTERVALS. INSTEAD JUST PASS THE INFORMATION IN AS ARGUMENTS...
% % I'LL NEED: MIN AND MAX SIZES TO INCLUDE FROM DATA, AND NUMBER OF MODELLED
% % SIZE CLASSES. IF NUMBER OF MODELLED SIZE CLASSES IS NOT SPECIFIED THEN
% % RETURN A LIST OF DATA SETS, WHERE EACH LIST ELEMENT IS THE DATA PARTITION
% % BY A DIFFERENT NUMBER OF INTERVALS...
% 
% 
% dx = mean(diff(log10(dat_size2.ESD))); % distance between points in (log scale) size spectra
% lESDmin = log10(ESDmin);
% lESDmax = log10(ESDmax);
% 
% logSize = log10(dat_size2.ESD);
% logVal = log10(dat_size2.(bvar));
% 
% AIC = nan(nbinMax-nbinMin+1,1); % use AIC metric to choose number of size class intervals
% 
% for i = nbinMin:nbinMax % loop over number of size class intervals
%     figure
%     plot(logSize, logVal, 'LineWidth', 2)
%     gc = gca;
%     gc.XLim = [lESDmin, lESDmax];
%     hold on
%     
%     w = (lESDmax - lESDmin) / i; % interval width
%     b0 = linspace(lESDmin, lESDmax, i+1); % boundaries
%     m0 = 0.5 * (b0(1:end-1) + b0(2:end)); % midpoints
%     
%     for j = 1:length(b0)
%         line([b0(j), b0(j)], gc.YLim);
%     end
%     
%     d = dat_size2;
%     d.sizeClass = nan(height(d),1);
%     for j = i:-1:1
%         d.sizeClass(logSize <= b0(j+1)) = j;
%     end
%     
%     meanVals = nan(1,i);
%     for j = 1:i
%        meanVals(j) = mean(logVal(d.sizeClass == j));
%     end
%     
%     % add piece-wise approx to plot
%     for j = 1:i
%         line([b0(j) b0(j+1)], [meanVals(j) meanVals(j)], 'Color', [1 0 0])
%     end
%     
%     % calculate goodness-of-fit metric    
%     logLik = nan(1,i);
%     for j = 1:i
%         err = logVal(d.sizeClass == j) - meanVals(j);
%         err_s = std(err);
%         logLik(j) = (1 ./ (floor(w/dx))-1) .* sum(-0.5*log(2*pi) - log(err_s) -0.5 .* (err ./ err_s) .^ 2);
%     end
%     AIC(i-nbinMin+1) = 2*i - 2*sum(logLik);
%     
%     
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

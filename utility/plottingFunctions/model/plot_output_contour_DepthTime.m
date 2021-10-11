function plot_output_contour_DepthTime(var, Out, AuxVars, ...
    fixedParams, forcing, varargin)

if forcing.nTraj ~= size(Out.N, ndims(Out.N))
    error('Number of trajectories contained in forcing data struct must match the number of trajectories contained in model output struct.')
end

extractVarargin(varargin)

if ~exist('xlab', 'var')
    xlab = 'year-day';
end
if ~exist('ylab', 'var')
    ylab = 'depth (m)';
end
if ~exist('smooth', 'var')
    smooth = 'linear';
end
if ~exist('ensurePositive', 'var')
    ensurePositive = true; % interpolation can produce negative values which are reset
end

% Index trajectories to plot
if exist('waterOrigin', 'var')
    waterOrigin = eval('waterOrigin');
    if ~contains(waterOrigin, forcing.waterMass)
        error('Optional argument "waterOrigin" must match one of the "waterMass" values in forcing data.')
    end
    traj = find(strcmp(forcing.waterMass, waterOrigin));
elseif exist('traj', 'var')
    traj = eval('traj');
    if ~islogical(traj)
        if any(~ismember(traj, 1:forcing.nTraj))
            error('Values in optional indexing argument "traj" must not exceed number of trajectories stored in forcing data.')
        end
    else
        if length(traj) ~= forcing.nTraj
            error('Length of logical indexing argument "traj" must equal number of trajectories in forcing data.')
        end
        traj = find(traj);
    end
else
    traj = 1:forcing.nTraj;
end


% Extract outputs
N = squeeze(Out.N(:,:,:,traj));
P = squeeze(Out.P(:,:,:,:,traj));
Z = squeeze(Out.Z(:,:,:,:,traj));
OM = squeeze(Out.OM(:,:,:,:,traj));

biovolume = AuxVars.biovolume;
P_biovolume = biovolume(1:fixedParams.nPP_size,:);
Z_biovolume = biovolume(fixedParams.nPP_size+1:end,:);
P_biovolume = reshape(P_biovolume, [fixedParams.nPP_size, size(biovolume, 2:ndims(biovolume))]);
Z_biovolume = reshape(Z_biovolume, [fixedParams.nZP_size, size(biovolume, 2:ndims(biovolume))]);


% If multiple trajectories selected then take averages
ntraj = length(traj);
if ntraj > 1
    N = mean(N, ndims(N), 'omitnan');
    P = mean(P, ndims(P), 'omitnan');
    Z = mean(Z, ndims(Z), 'omitnan');
    OM = mean(OM, ndims(OM), 'omitnan');
    P_biovolume = mean(P_biovolume, ndims(P_biovolume), 'omitnan');
    Z_biovolume = mean(Z_biovolume, ndims(Z_biovolume), 'omitnan');
end

% fixedParams.nPP_size = double(fixedParams.nPP_size);
% fixedParams.nZP_size = double(fixedParams.nZP_size);


% Time-depth grid for interpolation
[depth, time] = ndgrid(abs(fixedParams.z), 1:fixedParams.nt);
[depthGrid, timeGrid] = ndgrid(1:1:abs(fixedParams.zw(end)), 1:fixedParams.nt);

if exist('xLim', 'var')
    depthGrid = depthGrid(:,xLim(1):xLim(2));
    timeGrid = timeGrid(:,xLim(1):xLim(2));
end
if exist('yLim', 'var')
    depthGrid = depthGrid(yLim(1):yLim(2),:);
    timeGrid = timeGrid(yLim(1):yLim(2),:);
end


switch var

    case 'DIN'
        x = N;
        F = griddedInterpolant(depth, time, x, smooth);
        Fsmooth = flip(F(depthGrid, timeGrid));
        if ensurePositive, Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0)); end
        contourf(Fsmooth)
        cb = colorbar;
        if ~exist('ColourBarLabel', 'var')
           ColourBarLabel =  'mmol N / m^3';
        end
        cb.Label.String = ColourBarLabel;
        xlabel(xlab)
        ylabel(ylab)
        xticks(100:100:fixedParams.nt)
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(abs(linspace(fixedParams.zw(end),0,7)))
        if ~exist('Title', 'var')
            Title = 'DIN';
            if exist('waterOrigin', 'var')
                Title = [Title ': ' waterOrigin ' waters'];
            end
        end
        title(Title)

        %------------------------------------------------------------------
        
    case {'DOC','DON','POC','PON'}
        
        switch var
            case 'DOC', x = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:));
            case 'DON', x = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:));
            case 'POC', x = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:));
            case 'PON', x = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:));
        end
        
        F = griddedInterpolant(depth, time, x, smooth);
        Fsmooth = flip(F(depthGrid, timeGrid));
        if ensurePositive, Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0)); end
        contourf(Fsmooth)
        cb = colorbar;
        if ~exist('ColourBarLabel', 'var')
            ColourBarLabel =  'mmol N / m^3';
        end
        cb.Label.String = ColourBarLabel;
        xlabel(xlab)
        ylabel(ylab)
        xticks(100:100:fixedParams.nt)
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(abs(linspace(fixedParams.zw(end),0,7)))
        if ~exist('Title', 'var')
            Title = var;
            if exist('waterOrigin', 'var')
                Title = [Title ': ' waterOrigin ' waters'];
            end
        end
        title(Title)
        
    case {'P_N', 'P_Chl', 'P_C', 'P_biovolume', 'Z_N', 'Z_C', 'Z_biovolume', 'P_N_C', 'P_Chl_N', 'Z_N_C'}
        if ~exist('sizeClass', 'var')
            error('Optional argument "sizeClass" must be specified')
        else
            if ~ismember(sizeClass, 1:fixedParams.nPP_size)
                error('Value of optional argument "sizeClass" must be within interval [1, FixedParams.nPP_size]')
            end
        end
        switch var
            case 'P_N'
                x = squeeze(P(sizeClass,:,fixedParams.PP_N_index,:));
                ColourBarLabel_ =  'mmol N m^{-3}';
            case 'P_Chl'
                x = squeeze(P(sizeClass,:,fixedParams.PP_Chl_index,:));
                ColourBarLabel_ =  'mg Chl a m^{-3}';
            case 'P_C'
                x = squeeze(P(sizeClass,:,fixedParams.PP_C_index,:));
                ColourBarLabel_ =  'mmol C m^{-3}';
            case 'Z_N'
                x = squeeze(Z(sizeClass,:,fixedParams.ZP_N_index,:));
                ColourBarLabel_ =  'mmol N m^{-3}';
            case 'Z_C'
                x = squeeze(Z(sizeClass,:,fixedParams.ZP_C_index,:));
                ColourBarLabel_ =  'mmol C m^{-3}';
            case 'P_N_C'
                x = squeeze(P(sizeClass,:,fixedParams.PP_N_index,:) ./ ...
                    P(sizeClass,:,fixedParams.PP_C_index,:));
                ColourBarLabel_ =  'mmol N (mmol C)^{-1}';
            case 'P_Chl_N'
                x = squeeze(P(sizeClass,:,fixedParams.PP_Chl_index,:) ./ ...
                    P(sizeClass,:,fixedParams.PP_N_index,:));
                ColourBarLabel_ =  'mg Chl a (mmol N)^{-1}';
            case 'Z_N_C'
                x = squeeze(Z(sizeClass,:,fixedParams.ZP_N_index,:) ./ ...
                    Z(sizeClass,:,fixedParams.ZP_C_index,:));
                ColourBarLabel_ =  'mmol N (mmol C)^{-1}';
            case 'P_biovolume'
                x = 1e-9 .* squeeze(P_biovolume(sizeClass,:,:));
                ColourBarLabel_ =  'mm^3 m^{-3}';
            case 'Z_biovolume'
                x = 1e-9 .* squeeze(Z_biovolume(sizeClass,:,:));
                ColourBarLabel_ =  'mm^3 m^{-3}';
        end
        F = griddedInterpolant(depth, time, x, smooth);
        Fsmooth = flip(F(depthGrid, timeGrid));
        if ensurePositive, Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0)); end
        contourf(Fsmooth)
        cb = colorbar;
        if ~exist('ColourBarLabel', 'var')
            ColourBarLabel =  ColourBarLabel_;
        end
        cb.Label.String = ColourBarLabel;
        if exist('xLim', 'var')
            xt = get(gca, 'XTick');
            xt = xt + xLim(1);
            set(gca, 'XTickLabel', xt)
        end
        
        if exist('yLim', 'var')
            yt = get(gca, 'YTick');
            yt = yLim(2) - yt;
            set(gca, 'YTickLabel', yt)
        end
        
        xlabel(xlab)
        ylabel(ylab)
        
        trophicLevel = var(1);
        switch trophicLevel
            case 'P', sizeInterval = fixedParams.PPdia_intervals(sizeClass:sizeClass+1);
            case 'Z', sizeInterval = fixedParams.ZPdia_intervals(sizeClass:sizeClass+1);
        end
        sizeInterval = round(sizeInterval, 2, 'significant');
        if ~exist('Title', 'var')
            Title = ['[' num2str(sizeInterval(1)) ', ' num2str(sizeInterval(2)) '] \mum'];
        end
        title(Title)
        
        if exist('sizeLab', 'var') && eval('sizeLab')
            if ~exist('sizeLabPosition', 'var')
                sizeLabPosition = [0.1, 0.8, 0.3, 0.15];
            end
            xl = get(gca, 'XLim');
            yl = get(gca, 'YLim');
            xd = diff(xl); yd = diff(yl);
            recPos = [xl(1) + sizeLabPosition(1) * xd, ...
                yl(1) + sizeLabPosition(2) * yd, ...
                sizeLabPosition(3) * xd, ...
                sizeLabPosition(4) * yd];
            rectangle('Position', recPos, 'FaceColor', [1 1 1]);
            recTxt = ['[' num2str(sizeInterval(1)) ', ' num2str(sizeInterval(2)) '] \mum '];
            txtPos = [recPos(1) + recPos(3), recPos(2) + 0.5 * recPos(4)];
            txt = text(txtPos(1), txtPos(2), recTxt);
            set(txt, {'HorizontalAlignment', 'VerticalAlignment'}, {'right', 'middle'})
        end
        
        
        
        
        
end
        
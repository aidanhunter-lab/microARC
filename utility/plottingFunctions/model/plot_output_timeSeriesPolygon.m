function plot_output_timeSeriesPolygon(var, Out, AuxVars, fixedParams, ...
    forcing, dat, varargin)

extractVarargin(varargin)

if ~exist('depth', 'var')
    depth = 'mean';
end
if ~exist('polygonColour', 'var')
    polygonColour = [0 0 0.5];
end
if ~exist('lineColour', 'var')
    lineColour = [0 0 0];
end

if ~exist('highlightColour', 'var')
    highlightColour = [0 1 0]; % green as default
end
if ~exist('plotNew', 'var')
    plotNew = true;
end
if ~exist('fixedYaxis', 'var')
    fixedYaxis = false;
end

% if ~exist('axesTextSize', 'var')
%     axesTextSize = 12;
% end
if ~exist('titleTextSize', 'var')
    titleTextSize = 10;
end
if ~exist('legendTextSize', 'var')
    legendTextSize = 10;
end
if ~exist('legendPointSize', 'var')
    legendPointSize = 36;
end

if exist('waterMass', 'var')
    waterMass = eval('waterMass');
    if ~any(strcmp(waterMass, forcing.waterMass))
        error('Optional argument "waterMass" does not match values in forcing data.')
    end
    traj = find(strcmp(waterMass, forcing.waterMass));
else
%     waterMass = [];
    if ~exist('traj', 'var')
        traj = 1:forcing.nTraj;
    end
end
singleTrajectory = length(traj) == 1;

if singleTrajectory
    warning('Only a single trajectory has been selected. Multiple trajectories are required to define polygon.')
end


%%

nt = fixedParams.nt;
nz = fixedParams.nz;
nPP_size = double(fixedParams.nPP_size);
nZP_size = double(fixedParams.nZP_size);

x = yearday(forcing.t(:,1)); % yearday on x-axis

% if ~isempty(event)
%     etime = unique(dat.scalar.Yearday(dat.scalar.Event == event,:)); % yearday of sampling event
% end

% Extract outputs
N = squeeze(Out.N(:,:,:,traj));
P = Out.P(:,:,:,:,traj);
Z = Out.Z(:,:,:,:,traj);
OM = Out.OM(:,:,:,:,traj);

biovolume = 1e-9 .* AuxVars.biovolume(:,:,:,traj); % convert mu m^3 to mm^3
P_biovolume = biovolume(1:fixedParams.nPP_size,:);
Z_biovolume = biovolume(fixedParams.nPP_size+1:end,:);
P_biovolume = reshape(P_biovolume, [fixedParams.nPP_size, size(biovolume, 2:ndims(biovolume))]);
Z_biovolume = reshape(Z_biovolume, [fixedParams.nZP_size, size(biovolume, 2:ndims(biovolume))]);



stackedPlot = contains(var, 'stacked');

switch stackedPlot

    case false

        switch var
            case 'DIN'
                y = N;
                ylab = 'DIN (mmol N m^{-3})';
            case 'DON'
                y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:,:));
                ylab = 'DON (mmol N m^{-3})';
            case 'DOC'
                y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:,:));
                ylab = 'DOC (mmol C m^{-3})';
            case 'PON'
                y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:,:));
                ylab = 'PON (mmol N m^{-3})';
            case 'POC'
                y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:,:));
                ylab = 'POC (mmol C m^{-3})';
            case 'P_N'
                y = P(:,:,fixedParams.PP_N_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'autotroph biomass (mmol N m^{-3})';
            case 'P_C'
                y = P(:,:,fixedParams.PP_C_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'autotroph biomass (mmol C m^{-3})';
            case 'P_Chl'
                y = P(:,:,fixedParams.PP_Chl_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'autotroph biomass (mg chl a m^{-3})';
            case 'Z_N'
                y = Z(:,:,fixedParams.ZP_N_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'heterotroph biomass (mmol N m^{-3})';
            case 'Z_C'
                y = Z(:,:,fixedParams.ZP_C_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'heterotroph biomass (mmol C m^{-3})';
            case 'P_N_C'
                y = P(:,:,fixedParams.PP_N_index,:,:) ./ P(:,:,fixedParams.PP_C_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'autotroph quota (mmol N (mmol C)^{-1})';
            case 'Z_N_C'
                y = Z(:,:,fixedParams.ZP_N_index,:,:) ./ Z(:,:,fixedParams.ZP_C_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'heterotroph quota (mmol N (mmol C)^{-1})';
            case 'P_Chl_N'
                y = P(:,:,fixedParams.PP_Chl_index,:,:) ./ P(:,:,fixedParams.PP_N_index,:,:);
                if ~exist('sizeClass', 'var')
                    y = squeeze(sum(y));
                    warning('Optional argument "sizeClass" was not specified so the sum over size classes was plotted.')
                else
                    y = squeeze(y(sizeClass,:,:,:,:));
                end
                ylab = 'autotroph quota (mg Chl (mmol N)^{-1})';
        end
        
        y = depthAggregate(y, depth, fixedParams);
        
        switch singleTrajectory
            case true % line plot
                plot(x, y, 'Color', [0 0 0])
            case false % polygon
                lo = min(y, [], ndims(y));
                hi = max(y, [], ndims(y));
                ym = median(y, ndims(y), 'omitnan');
                lo = lo(:); hi = hi(:); ym = ym(:);
                tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
                pgon = polyshape(tt{:,1}, tt{:,2});
                plot(pgon, 'FaceColor', polygonColour);
                set(gca, {'XLim', 'YLim'}, {[min(x) max(x)] , [min(lo) max(hi)]})
                hold on
                plot(x,ym,'Color',lineColour)
                hold off
        end
        
        if ~exist('xlab', 'var')
            xlab = 'year-day';
        end
        if strcmp(depth, 'surface')
            ylab = {ylab, 'surface layer'};
        elseif strcmp(depth, 'mean')
            ylab = {ylab, 'depth averaged'};
        elseif strcmp(depth, 'max')
            ylab = {ylab, 'max over depths'};
        end
        xlabel(xlab)
        ylabel(ylab)
        
        if exist('sizeClass', 'var')
            switch var(1)
                case 'P', sizeInterval = fixedParams.PPdia_intervals(sizeClass:sizeClass+1);
                case 'Z', sizeInterval = fixedParams.ZPdia_intervals(sizeClass:sizeClass+1);
                otherwise, sizeInterval = [];
            end
            if ~isempty(sizeInterval)
                sizeInterval = round(sizeInterval, 2, 'significant');
            end
        end
        
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
        
    case true
        
        switch var
            
            case 'PZ_C_stacked'
                Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
                Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
                Y = squeeze(sum(Y,2,'omitnan')); % sum over depth
                YZ = squeeze(Z(:,:,fixedParams.ZP_C_index,:,:));
                if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
                YZ = repmat(reshape(fixedParams.zwidth, [1 nz]), [nZP_size 1]) .* YZ;
                YZ = squeeze(sum(YZ,2,'omitnan'));
                if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
                Y = cat(1, Y, YZ);
                
                if size(Y, ndims(Y)) == length(traj)
                    % average over trajectories
                    Y = mean(Y, ndims(Y));
                end
                
                
                % polygon vertices
                xpgon = [x(:); flip(x(:))];
                ypgon = zeros(nPP_size + nZP_size + 1, nt);
                ypgon(2:end,:) = Y;
                ypgonc = cumsum(ypgon);
                % 2 colour scales: green for autotrophs & red for heterotrophs
                cols = zeros(nPP_size+nZP_size, 3);
                if ~exist('colsP', 'var')
                    colsP = [0 100 0; 127 255 212] ./ 255;
                end
                if ~exist('colsZ', 'var')
                    colsZ = [0 100 0; 127 255 212] ./ 255;
                end
                cols(1:nPP_size, 1) = linspace(colsP(1,1), colsP(2,1), nPP_size);
                cols(1:nPP_size, 2) = linspace(colsP(1,2), colsP(2,2), nPP_size);
                cols(1:nPP_size, 3) = linspace(colsP(1,3), colsP(2,3), nPP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 1) = linspace(colsZ(1,1), colsZ(2,1), nZP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 2) = linspace(colsZ(1,2), colsZ(2,2), nZP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 3) = linspace(colsZ(1,3), colsZ(2,3), nZP_size);
                fill(xpgon, [ypgonc(1,:), flip(ypgonc(2,:))], cols(1,:))
                xlim([min(x) max(x)])
                ylim([0 max(ypgonc(:))])
                xl = xlim; yl = ylim;
                yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.25*diff(yl), max([nPP_size, nZP_size]) + 2);
                hold on
                
                text(xl(1)+0.12*diff(xl), yleg(1), 'cell diameter (\mum)', 'FontSize', legendTextSize)
                text(xl(1)+0.05*diff(xl), yleg(2), 'autotrophs', 'FontSize', legendTextSize)
                text(xl(1)+0.25*diff(xl), yleg(2), 'heterotrophs', 'FontSize', legendTextSize)
                
                ps = scatter(xl(1)+0.05*diff(xl), yleg(3), 'filled');
                ps.MarkerFaceColor = cols(1,:);
                ps.SizeData = legendPointSize;
                text(xl(1)+0.07*diff(xl), yleg(3), ...
                    num2str(round(fixedParams.PPdia(1),2,'significant')), ...
                    'FontSize', legendTextSize)
                for i = 2:nPP_size
                    fill(xpgon, [ypgonc(i,:), flip(ypgonc(i+1,:))], cols(i,:))
                    ps = scatter(xl(1)+0.05*diff(xl), yleg(i+2), 'filled');
                    ps.MarkerFaceColor = cols(i,:);
                    ps.SizeData = legendPointSize;
                    text(xl(1)+0.07*diff(xl), yleg(i+2), ...
                        num2str(round(fixedParams.PPdia(i),2,'significant')), ...
                        'FontSize', legendTextSize)
                end
                for i = 1:nZP_size
                    fill(xpgon, [ypgonc(i+nPP_size,:), flip(ypgonc(i+nPP_size+1,:))], cols(i+nPP_size,:))
                    ps = scatter(xl(1)+0.25*diff(xl), yleg(i+2), 'filled');
                    %             ps = scatter(xl(1)+0.15*diff(xl), yleg(i+2), 'filled');
                    ps.MarkerFaceColor = cols(i+nPP_size,:);
                    ps.SizeData = legendPointSize;
                    text(xl(1)+0.27*diff(xl), yleg(i+2), ...
                        num2str(round(fixedParams.ZPdia(i),2,'significant')), ...
                        'FontSize', legendTextSize)
                end
                hold off
                
                if exist('axisTextSize', 'var')
                    set(gca, 'FontSize', axesTextSize)
                end
                
                
                if ~exist('xlab', 'var')
                    xlab = 'year-day';
                end
                if ~exist('ylab', 'var')
                    ylab = 'biomass (mmol C)';
                end
                xlabel(xlab)
                ylabel(ylab)
                
                if ~exist('Title', 'var')
                    Title = ['plankton carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'];
                    if exist('waterMass', 'var')
                        if strcmp(waterMass, 'Arctic/Atlantic')
                            waterMass = 'Arctic & Atlantic';
                        end
                        Title = {Title, [waterMass, ' water']};
                    end
                end
                
                title(Title, 'FontSize', titleTextSize)

            case 'PZ_N_stacked'
                Y = squeeze(P(:,:,fixedParams.PP_N_index,:,:));
                Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
                Y = squeeze(sum(Y,2,'omitnan')); % sum over depth
                YZ = squeeze(Z(:,:,fixedParams.ZP_N_index,:,:));
                if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
                YZ = repmat(reshape(fixedParams.zwidth, [1 nz]), [nZP_size 1]) .* YZ;
                YZ = squeeze(sum(YZ,2,'omitnan'));
                if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
                Y = cat(1, Y, YZ);
                
                if size(Y, ndims(Y)) == length(traj)
                    % average over trajectories
                    Y = mean(Y, ndims(Y));
                end
                
                
                % polygon vertices
                xpgon = [x(:); flip(x(:))];
                ypgon = zeros(nPP_size + nZP_size + 1, nt);
                ypgon(2:end,:) = Y;
                ypgonc = cumsum(ypgon);
                % 2 colour scales: green for autotrophs & red for heterotrophs
                cols = zeros(nPP_size+nZP_size, 3);
                if ~exist('colsP', 'var')
                    colsP = [0 100 0; 127 255 212] ./ 255;
                end
                if ~exist('colsZ', 'var')
                    colsZ = [0 100 0; 127 255 212] ./ 255;
                end
                cols(1:nPP_size, 1) = linspace(colsP(1,1), colsP(2,1), nPP_size);
                cols(1:nPP_size, 2) = linspace(colsP(1,2), colsP(2,2), nPP_size);
                cols(1:nPP_size, 3) = linspace(colsP(1,3), colsP(2,3), nPP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 1) = linspace(colsZ(1,1), colsZ(2,1), nZP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 2) = linspace(colsZ(1,2), colsZ(2,2), nZP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 3) = linspace(colsZ(1,3), colsZ(2,3), nZP_size);
                fill(xpgon, [ypgonc(1,:), flip(ypgonc(2,:))], cols(1,:))
                xlim([min(x) max(x)])
                ylim([0 max(ypgonc(:))])
                xl = xlim; yl = ylim;
                yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.25*diff(yl), max([nPP_size, nZP_size]) + 2);
                hold on
                
                text(xl(1)+0.12*diff(xl), yleg(1), 'cell diameter (\mum)', 'FontSize', legendTextSize)
                text(xl(1)+0.05*diff(xl), yleg(2), 'autotrophs', 'FontSize', legendTextSize)
                text(xl(1)+0.25*diff(xl), yleg(2), 'heterotrophs', 'FontSize', legendTextSize)
                
                ps = scatter(xl(1)+0.05*diff(xl), yleg(3), 'filled');
                ps.MarkerFaceColor = cols(1,:);
                ps.SizeData = legendPointSize;
                text(xl(1)+0.07*diff(xl), yleg(3), ...
                    num2str(round(fixedParams.PPdia(1),2,'significant')), ...
                    'FontSize', legendTextSize)
                for i = 2:nPP_size
                    fill(xpgon, [ypgonc(i,:), flip(ypgonc(i+1,:))], cols(i,:))
                    ps = scatter(xl(1)+0.05*diff(xl), yleg(i+2), 'filled');
                    ps.MarkerFaceColor = cols(i,:);
                    ps.SizeData = legendPointSize;
                    text(xl(1)+0.07*diff(xl), yleg(i+2), ...
                        num2str(round(fixedParams.PPdia(i),2,'significant')), ...
                        'FontSize', legendTextSize)
                end
                for i = 1:nZP_size
                    fill(xpgon, [ypgonc(i+nPP_size,:), flip(ypgonc(i+nPP_size+1,:))], cols(i+nPP_size,:))
                    ps = scatter(xl(1)+0.25*diff(xl), yleg(i+2), 'filled');
                    %             ps = scatter(xl(1)+0.15*diff(xl), yleg(i+2), 'filled');
                    ps.MarkerFaceColor = cols(i+nPP_size,:);
                    ps.SizeData = legendPointSize;
                    text(xl(1)+0.27*diff(xl), yleg(i+2), ...
                        num2str(round(fixedParams.ZPdia(i),2,'significant')), ...
                        'FontSize', legendTextSize)
                end
                hold off
                
                if exist('axextextSize', 'var')
                    set(gca, 'FontSize', axesTextSize)
                end
                
                if ~exist('xlab', 'var')
                    xlab = 'year-day';
                end
                if ~exist('ylab', 'var')
                    ylab = 'biomass (mmol N)';
                end
                xlabel(xlab)
                ylabel(ylab)
                
                if ~exist('Title', 'var')
                    Title = ['plankton nitrogen in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'];
                    if exist('waterMass', 'var')
                        if strcmp(waterMass, 'Arctic/Atlantic')
                            waterMass = 'Arctic & Atlantic';
                        end
                        Title = {Title, [waterMass, ' water']};
                    end
                end
                
                title(Title, 'FontSize', titleTextSize)
                
            case 'PZ_biomass_stacked'
                
                Y = P_biovolume;
                Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
                Y = squeeze(sum(Y,2,'omitnan')); % sum over depth
                YZ = Z_biovolume;
                if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
                YZ = repmat(reshape(fixedParams.zwidth, [1 nz]), [nZP_size 1]) .* YZ;
                YZ = squeeze(sum(YZ,2,'omitnan'));
                if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
                Y = cat(1, Y, YZ);
                if size(Y, ndims(Y)) == length(traj)
                    % average over trajectories
                    Y = mean(Y, ndims(Y));
                end
                
                
                % polygon vertices
                xpgon = [x(:); flip(x(:))];
                ypgon = zeros(nPP_size + nZP_size + 1, nt);
                ypgon(2:end,:) = Y;
                ypgonc = cumsum(ypgon);
                % 2 colour scales: green for autotrophs & red for heterotrophs
                cols = zeros(nPP_size+nZP_size, 3);
                if ~exist('colsP', 'var')
                    colsP = [0 100 0; 127 255 212] ./ 255;
                end
                if ~exist('colsZ', 'var')
                    colsZ = [0 100 0; 127 255 212] ./ 255;
                end
                cols(1:nPP_size, 1) = linspace(colsP(1,1), colsP(2,1), nPP_size);
                cols(1:nPP_size, 2) = linspace(colsP(1,2), colsP(2,2), nPP_size);
                cols(1:nPP_size, 3) = linspace(colsP(1,3), colsP(2,3), nPP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 1) = linspace(colsZ(1,1), colsZ(2,1), nZP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 2) = linspace(colsZ(1,2), colsZ(2,2), nZP_size);
                cols(nPP_size+1:nPP_size+nZP_size, 3) = linspace(colsZ(1,3), colsZ(2,3), nZP_size);
                fill(xpgon, [ypgonc(1,:), flip(ypgonc(2,:))], cols(1,:))
                xlim([min(x) max(x)])
                ylim([0 max(ypgonc(:))])
                xl = xlim; yl = ylim;
                yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.25*diff(yl), max([nPP_size, nZP_size]) + 2);
                hold on
                
                text(xl(1)+0.12*diff(xl), yleg(1), 'cell diameter (\mum)', 'FontSize', legendTextSize)
                text(xl(1)+0.05*diff(xl), yleg(2), 'autotrophs', 'FontSize', legendTextSize)
                text(xl(1)+0.25*diff(xl), yleg(2), 'heterotrophs', 'FontSize', legendTextSize)
                
                ps = scatter(xl(1)+0.05*diff(xl), yleg(3), 'filled');
                ps.MarkerFaceColor = cols(1,:);
                ps.SizeData = legendPointSize;
                text(xl(1)+0.07*diff(xl), yleg(3), ...
                    num2str(round(fixedParams.PPdia(1),2,'significant')), ...
                    'FontSize', legendTextSize)
                for i = 2:nPP_size
                    fill(xpgon, [ypgonc(i,:), flip(ypgonc(i+1,:))], cols(i,:))
                    ps = scatter(xl(1)+0.05*diff(xl), yleg(i+2), 'filled');
                    ps.MarkerFaceColor = cols(i,:);
                    ps.SizeData = legendPointSize;
                    text(xl(1)+0.07*diff(xl), yleg(i+2), ...
                        num2str(round(fixedParams.PPdia(i),2,'significant')), ...
                        'FontSize', legendTextSize)
                end
                for i = 1:nZP_size
                    fill(xpgon, [ypgonc(i+nPP_size,:), flip(ypgonc(i+nPP_size+1,:))], cols(i+nPP_size,:))
                    ps = scatter(xl(1)+0.25*diff(xl), yleg(i+2), 'filled');
                    %             ps = scatter(xl(1)+0.15*diff(xl), yleg(i+2), 'filled');
                    ps.MarkerFaceColor = cols(i+nPP_size,:);
                    ps.SizeData = legendPointSize;
                    text(xl(1)+0.27*diff(xl), yleg(i+2), ...
                        num2str(round(fixedParams.ZPdia(i),2,'significant')), ...
                        'FontSize', legendTextSize)
                end
                hold off
                
                if exist('axextextSize', 'var')
                    set(gca, 'FontSize', axesTextSize)
                end
                
                if ~exist('xlab', 'var')
                    xlab = 'year-day';
                end
                if ~exist('ylab', 'var')
                    ylab = 'biovolume (mm^3 m^{-2})';
                end
                xlabel(xlab)
                ylabel(ylab)
                
                if ~exist('Title', 'var')
                    Title = ['plankton biovolume in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'];
                    if exist('waterMass', 'var')
                        if strcmp(waterMass, 'Arctic/Atlantic')
                            waterMass = 'Arctic & Atlantic';
                        end
                        Title = {Title, [waterMass, ' water']};
                    end
                end
                
                title(Title, 'FontSize', titleTextSize)

                
        end

end






end
        
%     case 'DON'
%         if plotNew
%             plt = figure;
%             plt.Units = 'inches';
%             plt.Position = [0 0 12 4];
%         end
%         y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:,:));
%         y = depthAggregate(y, depth, fixedParams);
%         lo = min(y, [], ndims(y));
%         hi = max(y, [], ndims(y));
%         ym = median(y, ndims(y));
%         lo = lo(:); hi = hi(:); ym = ym(:);
%         tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
%         pgon = polyshape(tt{:,1}, tt{:,2});
%         plot(pgon)
%         xlim([min(x) max(x)])
%         ylim([min(lo) max(hi)])
%         yl = ylim;
%         hold on
%         plot(x,ym,'k')
%         plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
%         text(etime,yl(2)-0.1*diff(yl), ...
%             ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
%             'Color', highlightColour);
%         hold off
%         xlabel('year-day')
%         if isempty(depth)
%             ylabel('DON (mmol N m^{-3})')
%         elseif strcmp(depth, 'surface')
%             ylabel({'DON (mmol N m^{-3})', 'surface layer'})
%         elseif strcmp(depth, 'mean')
%             ylabel({'DON (mmol N m^{-3})', 'depth averaged'})
%         elseif strcmp(depth, 'max')
%             ylabel({'DON (mmol N m^{-3})', 'max over depths'})
%         end
%         
%     case 'DOC'
%         if plotNew
%             plt = figure;
%             plt.Units = 'inches';
%             plt.Position = [0 0 12 4];
%         end
%         y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:,:));
%         y = depthAggregate(y, depth, fixedParams);
%         lo = min(y, [], ndims(y));
%         hi = max(y, [], ndims(y));
%         ym = median(y, ndims(y));
%         lo = lo(:); hi = hi(:); ym = ym(:);
%         tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
%         pgon = polyshape(tt{:,1}, tt{:,2});
%         plot(pgon)
%         xlim([min(x) max(x)])
%         ylim([min(lo) max(hi)])
%         yl = ylim;
%         hold on
%         plot(x,ym,'k')
%         plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
%         text(etime,yl(2)-0.1*diff(yl), ...
%             ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
%             'Color', highlightColour);
%         hold off
%         xlabel('year-day')
%         if isempty(depth)
%             ylabel('DOC (mmol C m^{-3})')
%         elseif strcmp(depth, 'surface')
%             ylabel({'DOC (mmol C m^{-3})', 'surface layer'})
%         elseif strcmp(depth, 'mean')
%             ylabel({'DOC (mmol C m^{-3})', 'depth averaged'})
%         elseif strcmp(depth, 'max')
%             ylabel({'DOC (mmol C m^{-3})', 'max over depths'})
%         end
%         
%     case 'PON'
%         if plotNew
%             plt = figure;
%             plt.Units = 'inches';
%             plt.Position = [0 0 12 4];
%         end
%         y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:,:));
%         y = depthAggregate(y, depth, fixedParams);
%         lo = min(y, [], ndims(y));
%         hi = max(y, [], ndims(y));
%         ym = median(y, ndims(y));
%         lo = lo(:); hi = hi(:); ym = ym(:);
%         tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
%         pgon = polyshape(tt{:,1}, tt{:,2});
%         plot(pgon)
%         xlim([min(x) max(x)])
%         ylim([min(lo) max(hi)])
%         yl = ylim;
%         hold on
%         plot(x,ym,'k')
%         plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
%         text(etime,yl(2)-0.1*diff(yl), ...
%             ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
%             'Color', highlightColour);
%         hold off
%         xlabel('year-day')
%         if isempty(depth)
%             ylabel('PON (mmol N m^{-3})')
%         elseif strcmp(depth, 'surface')
%             ylabel({'PON (mmol N m^{-3})', 'surface layer'})
%         elseif strcmp(depth, 'mean')
%             ylabel({'PON (mmol N m^{-3})', 'depth averaged'})
%         elseif strcmp(depth, 'max')
%             ylabel({'PON (mmol N m^{-3})', 'max over depths'})
%         end
%         
%     case 'POC'
%         if plotNew
%             plt = figure;
%             plt.Units = 'inches';
%             plt.Position = [0 0 12 4];
%         end
%         y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:,:));
%         y = depthAggregate(y, depth, fixedParams);
%         lo = min(y, [], ndims(y));
%         hi = max(y, [], ndims(y));
%         ym = median(y, ndims(y));
%         lo = lo(:); hi = hi(:); ym = ym(:);
%         tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
%         pgon = polyshape(tt{:,1}, tt{:,2});
%         plot(pgon)
%         xlim([min(x) max(x)])
%         ylim([min(lo) max(hi)])
%         yl = ylim;
%         hold on
%         plot(x,ym,'k')
%         plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
%         text(etime,yl(2)-0.1*diff(yl), ...
%             ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
%             'Color', highlightColour);
%         hold off
%         xlabel('year-day')
%         if isempty(depth)
%             ylabel('POC (mmol C m^{-3})')
%         elseif strcmp(depth, 'surface')
%             ylabel({'POC (mmol C m^{-3})', 'surface layer'})
%         elseif strcmp(depth, 'mean')
%             ylabel({'POC (mmol C m^{-3})', 'depth averaged'})
%         elseif strcmp(depth, 'max')
%             ylabel({'POC (mmol C m^{-3})', 'max over depths'})
%         end
%         
%         %----------------------------------------------------------
%         
%     case 'phytoplankton_C'
%         plt = figure;
%         plt.Units = 'inches';
%         plt.Position = [0 0 24 16];
%         Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
%         Y_ndims = ndims(Y);
%         Y = permute(Y, [2 1 3:Y_ndims]); % permute depth to 1st dimension
%         Y = depthAggregate(Y, depth, fixedParams);
%         Y_ndims = Y_ndims - 1;
%         nc = floor(sqrt(nPP_size));
%         nr = ceil(nPP_size / nc);
%         index = reshape([1:nPP_size, ...
%             zeros(1, nc * nr - nPP_size)], [nc nr])';
%         YL = [0 max(Y(:))]; % same y-axis for all size classes
%         for ii = 1:nPP_size
%             subplot(nr,nc,index(ii))
%             y = Y(ii,:,:);
%             lo = min(y, [], Y_ndims);
%             hi = max(y, [], Y_ndims);
%             ym = median(y, Y_ndims);
%             lo = lo(:); hi = hi(:); ym = ym(:);            
%             tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
%             pgon = polyshape(tt{:,1}, tt{:,2});
%             plot(pgon)
%             xlim([min(x) max(x)])
%             if fixedYaxis
%                 ylim(YL)
%             else
%                 ylim([0 max(hi)])
%             end
%             yl = ylim;
%             hold on
%             plot(x,ym,'k')
%             plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
%             if ii==1
%                 text(min(x)+0.2*(max(x)-min(x)),yl(2)-0.1*diff(yl), ...
%                     ['sampling event ' num2str(event)], 'Color', highlightColour)
%             end
%             hold off
%             if ismember(index(ii),index(end,:))
%                 xlabel('year-day')
%             else
%                 xlabel([])
%             end
%             
%             if ismember(index(ii),index(:,1))
%                 
%                 if isempty(depth)
%                     ylabel('mmol C m^{-3}')
%                 elseif strcmp(depth, 'surface')
%                     ylabel({'mmol C m^{-3}', 'surface layer'})
%                 elseif strcmp(depth, 'mean')
%                     ylabel({'mmol C m^{-3}', 'depth averaged'})
%                 elseif strcmp(depth, 'max')
%                     ylabel({'mmol C m^{-3}', 'max over depths'})
%                 end
%             else
%                 ylabel([])
%             end
%             
%             title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
%             sgtitle('Phytoplankton concentration time series')
%             
%         end
%         
%         %----------------------------------------------------------
% 
%     case 'phytoplanktonStacked'
%         
%         plt = figure;
%         plt.Units = 'inches';
%         plt.Position = [0 0 12 4];
%         
%         Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
%         Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
%         Y = squeeze(sum(Y, 2));
% 
%         % polygon vertices
%         xpgon = [x(:); flip(x(:))];
%         ypgon = zeros(nPP_size + 1, nt);
%         ypgon(2:end,:) = mean(Y, ndims(Y));
%         ypgonc = cumsum(ypgon);
%         % color scale - red --> blue
%         rb_hsv = rgb2hsv([1 0 0; 0 0 1]);
%         r2b_hsv =  [linspace(rb_hsv(1,1), rb_hsv(2,1), nPP_size)' ...
%             linspace(rb_hsv(1,2), rb_hsv(2,2), nPP_size)' ...
%             linspace(rb_hsv(1,3), rb_hsv(2,3), nPP_size)'];
%         r2b = hsv2rgb(r2b_hsv);
%         cols = r2b;
%         pgon = polyshape(xpgon, [ypgonc(1,:) flip(ypgonc(2,:))]');
%         ps = plot(pgon);
%         ps.FaceColor = cols(1,:);
%         xlim([min(x) max(x)])
%         ylim([0 max(ypgonc(:))])
%         xl = xlim; yl = ylim;
%         yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.5*diff(yl), nPP_size + 1);
%         hold on
%         text(xl(1)+0.05*diff(xl), yleg(1), 'cell diameter')
%         ps = scatter(xl(1)+0.05*diff(xl), yleg(2), 'filled');
%         ps.MarkerFaceColor = cols(1,:);
%         text(xl(1)+0.07*diff(xl), yleg(2), ...
%             [num2str(round(fixedParams.PPdia(1),2,'significant')) ' \mum'])
%         for i = 2:nPP_size
%             pgon = polyshape(xpgon, [ypgonc(i,:) flip(ypgonc(i+1,:))]');
%             ps = plot(pgon);
%             ps.FaceColor = cols(i,:);
%             ps = scatter(xl(1)+0.05*diff(xl), yleg(i+1), 'filled');
%             ps.MarkerFaceColor = cols(i,:);
%             text(xl(1)+0.07*diff(xl), yleg(i+1), ...
%                 [num2str(round(fixedParams.PPdia(i),2,'significant')) ' \mum'])
%         end
%         hold off
%         xlabel('year-day')
%         ylabel('abundance (mmol C)')
%         title(['phytoplankton carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'])
% 
%         %----------------------------------------------------------
%         
%     case 'phytoZooPlanktonStacked'
%         
%         plt = figure;
%         plt.Units = 'inches';
%         plt.Position = [0 0 10 4*5/6];
%         
%         Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
%         Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
%         Y = squeeze(sum(Y,2,'omitnan')); % sum over depth
%         
%         YZ = squeeze(Z(:,:,fixedParams.ZP_C_index,:,:));
%         if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
%         YZ = repmat(reshape(fixedParams.zwidth, [1 nz]), [nZP_size 1]) .* YZ;
%         YZ = squeeze(sum(YZ,2,'omitnan'));
%         if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
%         
%         Y = cat(1, Y, YZ);
%         
%         % polygon vertices
%         xpgon = [x(:); flip(x(:))];
%         ypgon = zeros(nPP_size + nZP_size + 1, nt);
%         ypgon(2:end,:) = mean(Y, ndims(Y));
%         ypgonc = cumsum(ypgon);
%         % 2 colour scales: green for autotrophs & red for heterotrophs
%         cols = zeros(nPP_size+nZP_size, 3);
%         greens = [0 100 0; 127 255 212] ./ 255; % dark green to aquamarine
%         reds = [128 0 0; 255 255 0] ./ 255;   % maroon to yellow
%         cols(1:nPP_size, 1) = linspace(greens(1,1), greens(2,1), nPP_size);
%         cols(1:nPP_size, 2) = linspace(greens(1,2), greens(2,2), nPP_size);
%         cols(1:nPP_size, 3) = linspace(greens(1,3), greens(2,3), nPP_size);        
%         cols(nPP_size+1:nPP_size+nZP_size, 1) = linspace(reds(1,1), reds(2,1), nZP_size);
%         cols(nPP_size+1:nPP_size+nZP_size, 2) = linspace(reds(1,2), reds(2,2), nZP_size);
%         cols(nPP_size+1:nPP_size+nZP_size, 3) = linspace(reds(1,3), reds(2,3), nZP_size);        
%         fill(xpgon, [ypgonc(1,:), flip(ypgonc(2,:))], cols(1,:))
%         xlim([min(x) max(x)])
%         ylim([0 max(ypgonc(:))])
%         xl = xlim; yl = ylim;        
%         yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.25*diff(yl), max([nPP_size, nZP_size]) + 2);
% %         yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.5*diff(yl), max([nPP_size, nZP_size]) + 2);
%         hold on
%         
%         text(xl(1)+0.15*diff(xl), yleg(1), 'cell diameter (\mum)', 'FontSize', legendTextSize)
%         text(xl(1)+0.05*diff(xl), yleg(2), 'autotrophs', 'FontSize', legendTextSize)        
%         text(xl(1)+0.25*diff(xl), yleg(2), 'heterotrophs', 'FontSize', legendTextSize)                
% %         text(xl(1)+0.1*diff(xl), yleg(1), 'cell diameter')
% %         text(xl(1)+0.05*diff(xl), yleg(2), 'autotrophs')        
% %         text(xl(1)+0.15*diff(xl), yleg(2), 'heterotrophs')        
% 
%         ps = scatter(xl(1)+0.05*diff(xl), yleg(3), 'filled');
% %         ps = scatter(xl(1)+0.05*diff(xl), yleg(3), 'filled');
%         ps.MarkerFaceColor = cols(1,:);
%         ps.SizeData = legendPointSize;
%         text(xl(1)+0.07*diff(xl), yleg(3), ...
%             num2str(round(fixedParams.PPdia(1),2,'significant')), ...
%             'FontSize', legendTextSize)
% %         text(xl(1)+0.07*diff(xl), yleg(3), ...
% %             [num2str(round(fixedParams.PPdia(1),2,'significant')) ' \mum'], ...
% %             'FontSize', legendTextSize)
%         for i = 2:nPP_size
%             fill(xpgon, [ypgonc(i,:), flip(ypgonc(i+1,:))], cols(i,:))
%             ps = scatter(xl(1)+0.05*diff(xl), yleg(i+2), 'filled');
%             ps.MarkerFaceColor = cols(i,:);
%             ps.SizeData = legendPointSize;
%             text(xl(1)+0.07*diff(xl), yleg(i+2), ...
%                 num2str(round(fixedParams.PPdia(i),2,'significant')), ...
%                 'FontSize', legendTextSize)
% %             text(xl(1)+0.07*diff(xl), yleg(i+2), ...
% %                 [num2str(round(fixedParams.PPdia(i),2,'significant')) ' \mum'], ...
% %                 'FontSize', legendTextSize)
%         end        
%         for i = 1:nZP_size
%             fill(xpgon, [ypgonc(i+nPP_size,:), flip(ypgonc(i+nPP_size+1,:))], cols(i+nPP_size,:))
%             ps = scatter(xl(1)+0.25*diff(xl), yleg(i+2), 'filled');
% %             ps = scatter(xl(1)+0.15*diff(xl), yleg(i+2), 'filled');
%             ps.MarkerFaceColor = cols(i+nPP_size,:);
%             ps.SizeData = legendPointSize;
%             text(xl(1)+0.27*diff(xl), yleg(i+2), ...
%                 num2str(round(fixedParams.ZPdia(i),2,'significant')), ...
%                 'FontSize', legendTextSize)
% %             text(xl(1)+0.27*diff(xl), yleg(i+2), ...
% %                 [num2str(round(fixedParams.ZPdia(i),2,'significant')) ' \mum'], ...
% %                 'FontSize', legendTextSize)
% %             text(xl(1)+0.17*diff(xl), yleg(i+2), ...
% %                 [num2str(round(fixedParams.ZPdia(i),2,'significant')) ' \mum'])
%         end
%         hold off
%         
%         set(gca, 'FontSize', axesTextSize)
%         
%         xlabel('year-day')
%         ylabel('biomass (mmol C)')
%         
%         if ~exist('waterMass', 'var') || isempty(waterMass)
%             waterMass = dat.scalar.waterMass{event};
%         end
%         if strcmp(waterMass, 'Arctic/Atlantic')
%             waterMass = 'Arctic & Atlantic';
%         end
%         
%         if ~isempty(event)
%             title({['plankton carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'], ...
%                 [waterMass ' water mass: sample event ' num2str(event)]})
%         else
%             title({['plankton carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'], ...
%                 [waterMass ' water']}, 'FontSize', titleTextSize)
%         end
% end
% 
% end


function v = depthAggregate(x, depth, fixedParams)
xsize = size(x);
if strcmp(depth, 'surface')
    v = reshape(x(1,:), xsize(2:end));
end
if strcmp(depth, 'mean')
    v = squeeze(sum(fixedParams.zwidth .* x) ./ fixedParams.Htot);
end
if strcmp(depth, 'max')
    v = squeeze(max(x));
end
end

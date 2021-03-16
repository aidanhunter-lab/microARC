function plt = plot_timeSeries_trajectoryPolygon(var, event, traj, out, ...
    auxVars, fixedParams, forcing, dat, varargin)

plt = [];
depth = 'mean';
plotNew = true;
fixedYaxis = false;

if ~isempty(varargin)
    if any(strcmp(varargin, 'depth'))
        depth = varargin{find(strcmp(varargin, 'depth'))+1};
    end
    
    if any(strcmp(varargin, 'highlightColour'))
        highlightColour = varargin{find(strcmp(varargin, 'highlightColour'))+1};
    else
        highlightColour = [0 1 0]; % green as default
    end
    
    if any(strcmp(varargin, 'plotNew'))
        plotNew = varargin{find(strcmp(varargin, 'plotNew'))+1};
    end
    
    if any(strcmp(varargin, 'fixedYaxis'))
        fixedYaxis = varargin{find(strcmp(varargin, 'fixedYaxis'))+1};
    end
    
end

nt = fixedParams.nt;
nz = fixedParams.nz;
nPP_size = fixedParams.nPP_size;
nZP_size = fixedParams.nZP_size;

% Extract times
x = yearday(forcing.t(:,1)); % yearday on x-axis
etime = unique(dat.scalar.Yearday(dat.scalar.Event == event,:)); % yearday of sampling event

% Extract outputs
N = squeeze(out.N(:,:,:,traj));
% P = squeeze(out.P(:,:,:,:,traj));
% Z = squeeze(out.Z(:,:,:,:,traj));
P = out.P(:,:,:,:,traj);
Z = out.Z(:,:,:,:,traj);
OM = squeeze(out.OM(:,:,:,:,traj));


switch var
    
    % Forcing data
    case 'temperature'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = forcing.T(:,:,traj);
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y)); % min and max over trajectories
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y)); % median over trajectories
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon);
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour)
        hold off
        if isempty(depth)
            ylabel(['temperature (' char(176) 'C)'])
        elseif strcmp(depth, 'surface')
            ylabel({['temperature (' char(176) 'C)'], 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({['temperature (' char(176) 'C)'], 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({['temperature (' char(176) 'C)'], 'max over depths'})
        end
        
    case 'diffusivity'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = forcing.K_center(:,:,traj);
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        if isempty(depth)
            ylabel('diffusivity (m^2 day^{-1})')
        elseif strcmp(depth, 'surface')
            ylabel({'diffusivity (m^2 day^{-1})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'diffusivity (m^2 day^{-1})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'diffusivity (m^2 day^{-1})', 'max over depths'})
        end
        
    case 'PAR'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = forcing.PARsurf(:,:,traj);
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        xlabel('year-day')
        if isempty(depth)
            ylabel('PAR (\muEin day^{-1} m^{-2})')
        elseif strcmp(depth, 'surface')
            ylabel({'PAR (\muEin day^{-1} m^{-2})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'PAR (\muEin day^{-1} m^{-2})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'PAR (\muEin day^{-1} m^{-2})', 'max over depths'})
        end
        
        %----------------------------------------------------------
        
    case 'DIN'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = depthAggregate(N, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        xlabel('year-day')
        if isempty(depth)
            ylabel('DIN (mmol N m^{-3})')
        elseif strcmp(depth, 'surface')
            ylabel({'DIN (mmol N m^{-3})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'DIN (mmol N m^{-3})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'DIN (mmol N m^{-3})', 'max over depths'})
        end
        
        %----------------------------------------------------------
        
    case 'DON'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:,:));
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        xlabel('year-day')
        if isempty(depth)
            ylabel('DON (mmol N m^{-3})')
        elseif strcmp(depth, 'surface')
            ylabel({'DON (mmol N m^{-3})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'DON (mmol N m^{-3})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'DON (mmol N m^{-3})', 'max over depths'})
        end
        
    case 'DOC'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:,:));
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        xlabel('year-day')
        if isempty(depth)
            ylabel('DOC (mmol C m^{-3})')
        elseif strcmp(depth, 'surface')
            ylabel({'DOC (mmol C m^{-3})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'DOC (mmol C m^{-3})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'DOC (mmol C m^{-3})', 'max over depths'})
        end
        
    case 'PON'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:,:));
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        xlabel('year-day')
        if isempty(depth)
            ylabel('PON (mmol N m^{-3})')
        elseif strcmp(depth, 'surface')
            ylabel({'PON (mmol N m^{-3})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'PON (mmol N m^{-3})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'PON (mmol N m^{-3})', 'max over depths'})
        end
        
    case 'POC'
        if plotNew
            plt = figure;
            plt.Units = 'inches';
            plt.Position = [0 0 12 4];
        end
        y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:,:));
        y = depthAggregate(y, depth, fixedParams);
        lo = min(y, [], ndims(y));
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y));
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon)
        xlim([min(x) max(x)])
        ylim([min(lo) max(hi)])
        yl = ylim;
        hold on
        plot(x,ym,'k')
        plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
        text(etime,yl(2)-0.1*diff(yl), ...
            ['sampling event ' num2str(event)], 'HorizontalAlignment', 'right', ...
            'Color', highlightColour);
        hold off
        xlabel('year-day')
        if isempty(depth)
            ylabel('POC (mmol C m^{-3})')
        elseif strcmp(depth, 'surface')
            ylabel({'POC (mmol C m^{-3})', 'surface layer'})
        elseif strcmp(depth, 'mean')
            ylabel({'POC (mmol C m^{-3})', 'depth averaged'})
        elseif strcmp(depth, 'max')
            ylabel({'POC (mmol C m^{-3})', 'max over depths'})
        end
        
        %----------------------------------------------------------
        
    case 'phytoplankton_C'
        plt = figure;
        plt.Units = 'inches';
        plt.Position = [0 0 24 16];
        Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
        Y_ndims = ndims(Y);
        Y = permute(Y, [2 1 3:Y_ndims]); % permute depth to 1st dimension
        Y = depthAggregate(Y, depth, fixedParams);
        Y_ndims = Y_ndims - 1;
        nc = floor(sqrt(nPP_size));
        nr = ceil(nPP_size / nc);
        index = reshape([1:nPP_size, ...
            zeros(1, nc * nr - nPP_size)], [nc nr])';
        YL = [0 max(Y(:))]; % same y-axis for all size classes
        for ii = 1:nPP_size
            subplot(nr,nc,index(ii))
            y = Y(ii,:,:);
            lo = min(y, [], Y_ndims);
            hi = max(y, [], Y_ndims);
            ym = median(y, Y_ndims);
            lo = lo(:); hi = hi(:); ym = ym(:);            
            tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
            pgon = polyshape(tt{:,1}, tt{:,2});
            plot(pgon)
            xlim([min(x) max(x)])
            if fixedYaxis
                ylim(YL)
            else
                ylim([0 max(hi)])
            end
            yl = ylim;
            hold on
            plot(x,ym,'k')
            plot([etime etime],[yl(1) yl(2)],'--', 'Color', highlightColour)
            if ii==1
                text(min(x)+0.2*(max(x)-min(x)),yl(2)-0.1*diff(yl), ...
                    ['sampling event ' num2str(event)], 'Color', highlightColour)
            end
            hold off
            if ismember(index(ii),index(end,:))
                xlabel('year-day')
            else
                xlabel([])
            end
            
            if ismember(index(ii),index(:,1))
                
                if isempty(depth)
                    ylabel('mmol C m^{-3}')
                elseif strcmp(depth, 'surface')
                    ylabel({'mmol C m^{-3}', 'surface layer'})
                elseif strcmp(depth, 'mean')
                    ylabel({'mmol C m^{-3}', 'depth averaged'})
                elseif strcmp(depth, 'max')
                    ylabel({'mmol C m^{-3}', 'max over depths'})
                end
            else
                ylabel([])
            end
            
            title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
            sgtitle('Phytoplankton concentration time series')
            
        end
        
        %----------------------------------------------------------

    case 'phytoplanktonStacked'
        
        plt = figure;
        plt.Units = 'inches';
        plt.Position = [0 0 12 4];
        
        Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
        Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
        Y = squeeze(sum(Y, 2));

        % polygon vertices
        xpgon = [x(:); flip(x(:))];
        ypgon = zeros(nPP_size + 1, nt);
        ypgon(2:end,:) = mean(Y, ndims(Y));
        ypgonc = cumsum(ypgon);
        % color scale - red --> blue
        rb_hsv = rgb2hsv([1 0 0; 0 0 1]);
        r2b_hsv =  [linspace(rb_hsv(1,1), rb_hsv(2,1), nPP_size)' ...
            linspace(rb_hsv(1,2), rb_hsv(2,2), nPP_size)' ...
            linspace(rb_hsv(1,3), rb_hsv(2,3), nPP_size)'];
        r2b = hsv2rgb(r2b_hsv);
        cols = r2b;
        pgon = polyshape(xpgon, [ypgonc(1,:) flip(ypgonc(2,:))]');
        ps = plot(pgon);
        ps.FaceColor = cols(1,:);
        xlim([min(x) max(x)])
        ylim([0 max(ypgonc(:))])
        xl = xlim; yl = ylim;
        yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.5*diff(yl), nPP_size + 1);
        hold on
        text(xl(1)+0.05*diff(xl), yleg(1), 'cell diameter')
        ps = scatter(xl(1)+0.05*diff(xl), yleg(2), 'filled');
        ps.MarkerFaceColor = cols(1,:);
        text(xl(1)+0.07*diff(xl), yleg(2), ...
            [num2str(round(fixedParams.PPdia(1),2,'significant')) ' \mum'])
        for i = 2:nPP_size
            pgon = polyshape(xpgon, [ypgonc(i,:) flip(ypgonc(i+1,:))]');
            ps = plot(pgon);
            ps.FaceColor = cols(i,:);
            ps = scatter(xl(1)+0.05*diff(xl), yleg(i+1), 'filled');
            ps.MarkerFaceColor = cols(i,:);
            text(xl(1)+0.07*diff(xl), yleg(i+1), ...
                [num2str(round(fixedParams.PPdia(i),2,'significant')) ' \mum'])
        end
        hold off
        xlabel('year-day')
        ylabel('abundance (mmol C)')
        title(['phytoplankton carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'])

        %----------------------------------------------------------

    case 'phytoZooPlanktonStacked'
        
        plt = figure;
        plt.Units = 'inches';
        plt.Position = [0 0 12 4];
        
        Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
        Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
        Y = squeeze(sum(Y,2));
        
        YZ = squeeze(Z(:,:,fixedParams.ZP_C_index,:,:));
        if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
        YZ = repmat(reshape(fixedParams.zwidth, [1 nz]), [nZP_size 1]) .* YZ;
        YZ = squeeze(sum(YZ,2));
        if nZP_size == 1, YZ = reshape(YZ, [1 size(YZ)]); end
        
        Y = cat(1, Y, YZ);
        
        % polygon vertices
        xpgon = [x(:); flip(x(:))];
        ypgon = zeros(nPP_size + 2, nt);
        ypgon(2:end,:) = mean(Y,3);
        ypgonc = cumsum(ypgon);
        % color scale - red --> blue
        rb_hsv = rgb2hsv([1 0 0; 0 0 1]);
        r2b_hsv =  [linspace(rb_hsv(1,1), rb_hsv(2,1), nPP_size)' ...
            linspace(rb_hsv(1,2), rb_hsv(2,2), nPP_size)' ...
            linspace(rb_hsv(1,3), rb_hsv(2,3), nPP_size)'];
        r2b = hsv2rgb(r2b_hsv);
        cols = r2b;
        cols(nPP_size + 1,:) = [0 0 0]; % black for zooplankton
        pgon = polyshape(xpgon, [ypgonc(1,:) flip(ypgonc(2,:))]');
        ps = plot(pgon);
        ps.FaceColor = cols(1,:);
        xlim([min(x) max(x)])
        ylim([0 max(ypgonc(:))])
        xl = xlim; yl = ylim;
        
        yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.5*diff(yl), nPP_size + 2);
        hold on
        text(xl(1)+0.05*diff(xl), yleg(1), 'cell diameter')
        ps = scatter(xl(1)+0.05*diff(xl), yleg(2), 'filled');
        ps.MarkerFaceColor = cols(1,:);
        text(xl(1)+0.07*diff(xl), yleg(2), ...
            [num2str(round(fixedParams.PPdia(1),2,'significant')) ' \mum'])
        for i = 2:nPP_size+1
            pgon = polyshape(xpgon, [ypgonc(i,:) flip(ypgonc(i+1,:))]');
            ps = plot(pgon);
            ps.FaceColor = cols(i,:);
            ps = scatter(xl(1)+0.05*diff(xl), yleg(i+1), 'filled');
            ps.MarkerFaceColor = cols(i,:);
            if i < nPP_size+1
                text(xl(1)+0.07*diff(xl), yleg(i+1), ...
                    [num2str(round(fixedParams.PPdia(i),2,'significant')) ' \mum'])
            else
                text(xl(1)+0.07*diff(xl), yleg(i+1), 'zooplankton')
            end
        end
        hold off
        xlabel('year-day')
        ylabel('abundance (mmol C)')
        title(['plankton carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'])

        
end

end


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

function plt = outputPlot(type, varargin)

switch type

    case 'contour_DepthTime'
        var = varargin{1}; % variables to plot
        itraj = varargin{2}; % trajectory
        out = varargin{3}; % model output
        fixedParams = varargin{4};
        forcing = varargin{5};
        auxVars = varargin{6};
        smooth = varargin{7};
        
        % Extract outputs
        N = squeeze(out.N(:,:,:,itraj));
        P = squeeze(out.P(:,:,:,:,itraj));
        Z = squeeze(out.Z(:,:,:,itraj));
        OM = squeeze(out.OM(:,:,:,:,itraj));
        
        % Time-depth grid for interpolation
        [depth, time] = ndgrid(abs(fixedParams.z), 1:fixedParams.nt);
        [depthGrid, timeGrid] = ndgrid(1:1:abs(fixedParams.zw(end)), 1:fixedParams.nt);
        % Should plots be smoothed by linear interpolation?
                
        switch var
            case 'forcing'
                plt = figure;                
                plt.Units = 'inches';
                plt.Position = [0 0 8 9];                
                subplot(3,1,1) % temperature
                x = forcing.T(:,:,itraj);
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = '\circC';
                title('Temperature')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                
                subplot(3,1,2); % diffusivity
                x = forcing.K_center(:,:,itraj);
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                Fsmooth(Fsmooth<=0) = nan;
                contourf(log10(Fsmooth))
                cb = colorbar;
                for ii = 1:length(cb.TickLabels), cb.TickLabels{ii} = string(round(10 ^ str2num(cb.TickLabels{ii}),2,'significant')); end
                cb.Label.String = 'm^2 day^{-1}';
                title('Diffusivity')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))                
                %                 set(gca,'colorscale','log')
                
                subplot(3,1,3) % PAR
                x = auxVars.PAR(:,:,itraj);
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                Fsmooth(Fsmooth<=0) = nan;
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = '\muEin day^{-1} m^{-2}';
                title('PAR')
                xlabel('year-day')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                colormap plasma
                
                %----------------------------------------------------------

            case 'inorganicNutrient'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 8 3];
                x = N;
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = 'mmol N / m^3';
                title('DIN')
                xlabel('year-day')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                colormap plasma
                
                %------------------------------------------------------------------

            case 'DOM_POM'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 2*8 2*3];
                subplot(2,2,1)
                x = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:));
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = 'mmol C / m^3';
                title('DOC')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                
                subplot(2,2,2)
                x = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:));
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = 'mmol N / m^3';
                title('DON')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                
                subplot(2,2,3)
                x = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:));
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = 'mmol C / m^3';
                title('POC')
                xlabel('year-day')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                
                subplot(2,2,4)
                x = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:));
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = 'mmol N / m^3';
                title('PON')
                xlabel('year-day')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                colormap plasma                
                
                %------------------------------------------------------------------
                
            case 'phytoplankton_N'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];                
                index = reshape([1:fixedParams.nPP_size, ... 
                    zeros(1, nc * nr - fixedParams.nPP_size)], [nc nr])';
                for ii = 1:fixedParams.nPP_size
                    subplot(nr,nc,index(ii))
                    x = squeeze(P(ii,:,fixedParams.PP_N_index,:));
                    F = griddedInterpolant(depth, time, x, smooth);
                    Fsmooth = flip(F(depthGrid, timeGrid));
                    Fsmooth(Fsmooth<0) = 0;
                    contourf(Fsmooth)
                    cb = colorbar;
                    cb.Label.String = 'mmol N / m^3';
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
%                     title([num2str(round(fixedParams.PPsize(ii),2,'significant')) ' \mum^3'])
                    xlabel('year-day')
                    ylabel('depth (m)')
                    xticks(100:100:fixedParams.nt)
                    xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                    yticks(linspace(0,abs(fixedParams.zw(end)),7))
                    yticklabels(linspace(fixedParams.zw(end),0,7))
                end
                sgtitle('phytoplankton nitrogen concentration given cell diameter')
%                 sgtitle('phytoplankton abundance given cell volume (mmol N / m^3)')
                colormap plasma
                
                %------------------------------------------------------------------
            
            case 'phytoplankton_Chl'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                index = reshape([1:fixedParams.nPP_size, ... 
                    zeros(1, nc * nr - fixedParams.nPP_size)], [nc nr])';
                for ii = 1:fixedParams.nPP_size
                    subplot(nr,nc,index(ii))
                    x = squeeze(P(ii,:,fixedParams.PP_Chl_index ,:));
                    F = griddedInterpolant(depth, time, x, smooth);
                    Fsmooth = flip(F(depthGrid, timeGrid));
                    Fsmooth(Fsmooth<0) = 0;
                    contourf(Fsmooth)
                    cb = colorbar;
                    cb.Label.String = 'mg Chl / m^3';
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
%                     title([num2str(round(fixedParams.PPsize(ii),2,'significant')) ' \mum^3'])
                    xlabel('year-day')
                    ylabel('depth (m)')
                    xticks(100:100:fixedParams.nt)
                    xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                    yticks(linspace(0,abs(fixedParams.zw(end)),7))
                    yticklabels(linspace(fixedParams.zw(end),0,7))
                end
                sgtitle('phytoplankton chlorophyll_a concentration given cell diameter')
%                 sgtitle('phytoplankton abundance given cell volume (mmol N / m^3)')
                colormap plasma
                
                %------------------------------------------------------------------
                
            case 'phytoplankton_C'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                index = reshape([1:fixedParams.nPP_size, ... 
                    zeros(1, nc * nr - fixedParams.nPP_size)], [nc nr])';
                for ii = 1:fixedParams.nPP_size
                    subplot(nr,nc,index(ii))
                    x = squeeze(P(ii,:,fixedParams.PP_C_index ,:));
                    F = griddedInterpolant(depth, time, x, smooth);
                    Fsmooth = flip(F(depthGrid, timeGrid));
                    Fsmooth(Fsmooth<0) = 0;
                    contourf(Fsmooth)
                    cb = colorbar;
                    cb.Label.String = 'mmol C / m^3';
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
                    xlabel('year-day')
                    ylabel('depth (m)')
                    xticks(100:100:fixedParams.nt)
                    xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                    yticks(linspace(0,abs(fixedParams.zw(end)),7))
                    yticklabels(linspace(fixedParams.zw(end),0,7))
                end
                sgtitle('phytoplankton carbon concentration given cell diameter')
                colormap plasma
                                
                %------------------------------------------------------------------

            case 'phytoplankton_N_C'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                index = reshape([1:fixedParams.nPP_size, ... 
                    zeros(1, nc * nr - fixedParams.nPP_size)], [nc nr])';
                for ii = 1:fixedParams.nPP_size
                    subplot(nr,nc,index(ii))
                    x = squeeze(P(ii,:,fixedParams.PP_N_index,:) ./ P(ii,:,fixedParams.PP_C_index,:));
                    F = griddedInterpolant(depth, time, x, smooth);
                    Fsmooth = flip(F(depthGrid, timeGrid));
                    Fsmooth(Fsmooth<0) = 0;
                    contourf(Fsmooth)
                    cb = colorbar;
                    cb.Label.String = 'mmol N / mmol C';
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
                    xlabel('year-day')
                    ylabel('depth (m)')
                    xticks(100:100:fixedParams.nt)
                    xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                    yticks(linspace(0,abs(fixedParams.zw(end)),7))
                    yticklabels(linspace(fixedParams.zw(end),0,7))
                end
                sgtitle('phytoplankton N/C ratio given cell diameter')
                colormap plasma
                
                %------------------------------------------------------------------
            
            case 'phytoplankton_Chl_N'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                index = reshape([1:fixedParams.nPP_size, ... 
                    zeros(1, nc * nr - fixedParams.nPP_size)], [nc nr])';
                for ii = 1:fixedParams.nPP_size
                    subplot(nr,nc,index(ii))
                    x = squeeze(P(ii,:,fixedParams.PP_Chl_index,:) ./ ... 
                        P(ii,:,fixedParams.PP_N_index,:));
                    F = griddedInterpolant(depth, time, x, smooth);
                    Fsmooth = flip(F(depthGrid, timeGrid));
                    Fsmooth(Fsmooth<0) = 0;
                    contourf(Fsmooth)
                    cb = colorbar;
                    cb.Label.String = 'mg Chl / mmol N';
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
                    xlabel('year-day')
                    ylabel('depth (m)')
                    xticks(100:100:fixedParams.nt)
                    xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                    yticks(linspace(0,abs(fixedParams.zw(end)),7))
                    yticklabels(linspace(fixedParams.zw(end),0,7))
                end
                sgtitle('phytoplankton Chl/N ratio given cell diameter')
                colormap plasma
                
                %------------------------------------------------------------------

            case 'zooplankton'                
                plt = figure;                
                plt.Units = 'inches';
                plt.Position = [0 0 8 3];
                x = Z;
                F = griddedInterpolant(depth, time, x, smooth);
                Fsmooth = flip(F(depthGrid, timeGrid));
                contourf(Fsmooth)
                cb = colorbar;
                cb.Label.String = 'mmol C / m^3';
                title('zooplankton abundance')
                xlabel('year-day')
                ylabel('depth (m)')
                xticks(100:100:fixedParams.nt)
                xticklabels(yearday(forcing.t(100:100:fixedParams.nt,itraj)))
                yticks(linspace(0,abs(fixedParams.zw(end)),7))
                yticklabels(linspace(fixedParams.zw(end),0,7))
                colormap plasma                
                plt = gcf;
        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
    case 'trajectoryLine_LatLong'
        
        var = varargin{1}; % variables to plot
        ie = varargin{2}; % event number
        itraj = varargin{3}; % trajectories
        out = varargin{4}; % model output
        fixedParams = varargin{5};
        forcing = varargin{6};
        auxVars = varargin{7};
        dat = varargin{8};
        alpha = varargin{9};
        
        % Extract positions and times
        lon = forcing.x(:,itraj);
        lat = forcing.y(:,itraj);
        day = yearday(forcing.t(:,itraj));
        x = lon; y = lat;
        x = [nan(1,size(x,2)); x; nan(1,size(x,2))]; % separate trajectories with nans
        y = [nan(1,size(y,2)); y; nan(1,size(y,2))];
        x = reshape(x, [1 numel(x)]);
        y = reshape(y, [1 numel(y)]);        
        z = zeros(size(x));

        % Extract outputs
        N = squeeze(out.N(:,:,:,itraj));
        P = squeeze(out.P(:,:,:,:,itraj));
        Z = squeeze(out.Z(:,:,:,itraj));
        OM = squeeze(out.OM(:,:,:,:,itraj));
        
        % Event position
        elat = dat.scalar.Latitude(dat.scalar.Event == ie);
        elon = dat.scalar.Longitude(dat.scalar.Event == ie);

        switch var
            case 'direction'                
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 6 4];
                col = day;
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);                
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = 'year-day';
                title('Trajectory directions')
                xlabel([char(176) 'E'])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off                
                colormap inferno
                
                %----------------------------------------------------------
                
            case 'forcing'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 6 12];

                subplot(3,1,1) % temperature
                col = squeeze(forcing.T(1,:,itraj));
                if length(itraj) == 1, col = col'; end 
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);                
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;                
                cb = colorbar;
                cb.Label.String = [char(176) 'C'];
                title('temperature')
                ylabel([char(176) 'N'])                
                hold on
                scatter(elon,elat,'pr','filled')
                hold off

                subplot(3,1,2) % diffusivity                
                col = squeeze(forcing.K(1,:,itraj));
                if length(itraj) == 1, col = col'; end
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);                
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;                
                cb = colorbar;
                cb.Label.String = 'm^2 day^{-1}';
                title('diffusivity')
                ylabel([char(176) 'N'])                
                hold on
                scatter(elon,elat,'pr','filled')
                hold off

                subplot(3,1,3)
                col = squeeze(forcing.PARsurf(1,:,itraj));
                if length(itraj) == 1, col = col'; end
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = '\muEin day^{-1} m^{-2}';
                title('PAR')
                xlabel([char(176) 'E'])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off
                sgtitle('Forcing data - surface layer')
                colormap plasma
                
                %----------------------------------------------------------
                
            case 'inorganicNutrient'
                plt = figure;                
                plt.Units = 'inches';
                plt.Position = [0 0 6 4];
                col = squeeze(sum(fixedParams.zwidth .* N) / fixedParams.Htot);
                if length(itraj) == 1, col = col'; end
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);                
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;                
                cb = colorbar;
                cb.Label.String = 'mmol N / m^3';
                title('DIN')
                xlabel([char(176) 'E'])
                ylabel([char(176) 'N'])                
                hold on
                scatter(elon,elat,'pr','filled')
                hold off                
                colormap plasma
                
                %----------------------------------------------------------

            case 'DOM_POM'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 2*6 2*4];
                subplot(2,2,1)
                col = squeeze(sum(fixedParams.zwidth .* squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:,:))) / fixedParams.Htot);
                if length(itraj) == 1, col = col'; end
                % col = squeeze(sum(OM));
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = 'mmol C / m^3';
                title('DOC')
                xlabel([])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off
                
                subplot(2,2,2)
                col = squeeze(sum(fixedParams.zwidth .* squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:,:))) / fixedParams.Htot);
                if length(itraj) == 1, col = col'; end
                % col = squeeze(sum(OM));
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = 'mmol N / m^3';
                title('DON')
                xlabel([])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off

                subplot(2,2,3)
                col = squeeze(sum(fixedParams.zwidth .* squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:,:))) / fixedParams.Htot);
                if length(itraj) == 1, col = col'; end
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = 'mmol C / m^3';
                title('POC')
                xlabel([char(176) 'E'])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off
                
                subplot(2,2,4)
                col = squeeze(sum(fixedParams.zwidth .* squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:,:))) / fixedParams.Htot);
                if length(itraj) == 1, col = col'; end
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = 'mmol N / m^3';
                title('PON')
                xlabel([char(176) 'E'])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off
                colormap plasma

                %----------------------------------------------------------
                
            case 'phytoplankton_N'
                plt = figure;
                nTraj = length(itraj);
                ind_nut = fixedParams.PP_N_index;
                Col = squeeze(max(P(:,:,ind_nut,:,:),[],2));
%                 Col = squeeze(max(P(:,:,ind_nut,:),[],2));
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 6*nc 4*nr];
                index = reshape([1:fixedParams.nPP_size, ...
                    zeros(1, nc * nr - fixedParams.nPP_size)], [nc nr])';
                for ii = 1:fixedParams.nPP_size
                    subplot(nr,nc,index(ii))
                    col = squeeze(Col(ii,:,:));
                    if nTraj == 1, col=col'; end
                    col = [nan(1,nTraj); col; nan(1,nTraj)];
                    col = reshape(col, [1 numel(col)]);
                    s = surface([x;x],[y;y],[z;z],[col;col], ...
                        'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                    s.EdgeAlpha = alpha;
                    cb = colorbar;
                    cb.Label.String = 'mmol N / m^3';
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
%                     title([num2str(round(fixedParams.PPsize(ii),2,'significant')) ' \mum^3'])
                    if ismember(index(ii),index(end,:)), xlabel([char(176) 'E']); else, xlabel([]); end
                    if ismember(index(ii),index(:,1)), ylabel([char(176) 'N']); else, ylabel([]); end
                    hold on
                    scatter(elon,elat,'pr','filled')
                    hold off
                end
                sgtitle('phytoplankton abundance / cell diameter')
%                 sgtitle('phytoplankton abundance / cell volume')
                colormap plasma
                
                %----------------------------------------------------------
                
            case 'zooplankton'                
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 6 4];
                col = squeeze(max(Z));
                if length(itraj) == 1, col = col'; end
                col = [nan(1,size(col,2)); col; nan(1,size(col,2))];
                col = reshape(col, [1 numel(col)]);
                s = surface([x;x],[y;y],[z;z],[col;col], ...
                    'facecol', 'no', 'edgecol', 'interp', 'linew', 1);
                s.EdgeAlpha = alpha;
                cb = colorbar;
                cb.Label.String = 'mmol C / m^3';
                title('zooplankton abundance')
                xlabel([char(176) 'E'])
                ylabel([char(176) 'N'])
                hold on
                scatter(elon,elat,'pr','filled')
                hold off
                colormap plasma

        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
    case 'trajectoryPolygon_TimeSeries'
        
        var = varargin{1}; % variables to plot
        ie = varargin{2}; % event number
        itraj = varargin{3}; % trajectories
        out = varargin{4}; % model output
        fixedParams = varargin{5};
        forcing = varargin{6};
        auxVars = varargin{7};
        dat = varargin{8};
        alpha = varargin{9};
        
        nt = fixedParams.nt;
        nz = fixedParams.nz;
        nPP_size = fixedParams.nPP_size;
        
        % Extract times
        day = yearday(forcing.t(:,itraj));
        x = day(:,1);
        etime = unique(dat.scalar.Yearday(dat.scalar.Event == ie,:));
        
        % Extract outputs
        N = squeeze(out.N(:,:,:,itraj));
        P = squeeze(out.P(:,:,:,:,itraj));
        Z = squeeze(out.Z(:,:,:,itraj));
        OM = squeeze(out.OM(:,:,:,:,itraj));

        
        switch var
            case 'forcing'
                plt = figure;                
                plt.Units = 'inches';
                plt.Position = [0 0 12 12];

                subplot(3,1,1)
                y = forcing.T(1,:,itraj);
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon);
%                 % remove duplicate vertices from polygon... probably not
%                 % necessary but MatLab complains with warning messages
%                 % otherwise...
%                 dupl = diff(lo) == 0 | diff(hi) == 0;
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                text(max(x)-0.1*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                plot([max(x)-0.15*(max(x)-min(x)), max(x)-0.11*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                hold off
                ylabel(['temperature (' char(176) 'C)'])

                subplot(3,1,2)
                y = forcing.K(1,:,itraj);
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)                
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                hold off
                ylabel('diffusivity (m^2 day^{-1})')

                subplot(3,1,3)
                y = forcing.PARsurf(1,:,itraj);
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                hold off
                xlabel('year-day')
                ylabel('PAR (\muEin day^{-1} m^{-2})')
                
                %----------------------------------------------------------

            case 'DIN'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 12 4];
                
                y = sum(fixedParams.zwidth .* N) / fixedParams.Htot;
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                text(max(x)-0.1*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                plot([max(x)-0.15*(max(x)-min(x)), max(x)-0.11*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                hold off
                xlabel('year-day')
                ylabel('DIN (mmol N m^{-3})')
                
                %----------------------------------------------------------

            case 'DOM_POM'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 2*12 2*4];
                
                subplot(2,2,1)                
                y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_C_index,:,:));                
                y = sum(fixedParams.zwidth .* y) / fixedParams.Htot;
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)                
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                text(max(x)-0.1*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                plot([max(x)-0.15*(max(x)-min(x)), max(x)-0.11*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                hold off
                ylabel('DOC (mmol C m^{-3})')

                subplot(2,2,2)
                y = squeeze(OM(fixedParams.DOM_index,:,fixedParams.OM_N_index,:,:));
                y = sum(fixedParams.zwidth .* y) / fixedParams.Htot;
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                text(max(x)-0.1*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                plot([max(x)-0.15*(max(x)-min(x)), max(x)-0.11*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                hold off
                ylabel('DON (mmol N m^{-3})')
                
                subplot(2,2,3)
                y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_C_index,:,:));
                y = sum(fixedParams.zwidth .* y) / fixedParams.Htot;
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)                
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                text(max(x)-0.1*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                plot([max(x)-0.15*(max(x)-min(x)), max(x)-0.11*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                hold off
                xlabel('year-day')
                ylabel('POC (mmol C m^{-3})')

                subplot(2,2,4)
                y = squeeze(OM(fixedParams.POM_index,:,fixedParams.OM_N_index,:,:));
                y = sum(fixedParams.zwidth .* y) / fixedParams.Htot;
                lo = min(y,[],3);
                hi = max(y,[],3);
                ym = median(y,3);
                lo = lo(:); hi = hi(:); ym = ym(:);
                pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                plot(pgon)
                xlim([min(x) max(x)])
                ylim([min(lo) max(hi)])
                yl = ylim;
                hold on
                plot(x,ym,'k')
                plot([etime etime],[yl(1) yl(2)],':k')
                text(max(x)-0.1*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                plot([max(x)-0.15*(max(x)-min(x)), max(x)-0.11*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                hold off
                xlabel('year-day')
                ylabel('PON (mmol N m^{-3})')

                sgtitle('organic matter')
                

                %----------------------------------------------------------
                
            case 'phytoplankton_C'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 12 12];                
                Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));                
                Y = squeeze(max(Y,[],2));
                nc = floor(sqrt(nPP_size));
                nr = ceil(nPP_size / nc);
                index = reshape([1:nPP_size, ... 
                    zeros(1, nc * nr - nPP_size)], [nc nr])';
                for ii = 1:nPP_size
                    subplot(nr,nc,index(ii))
                    y = Y(ii,:,:);
                    lo = min(y,[],3);
                    hi = max(y,[],3);
                    ym = median(y,3);
                    lo = lo(:); hi = hi(:); ym = ym(:);
                    pgon = polyshape([x; flip(x); x(1)], [lo; flip(hi); lo(1)]);
                    plot(pgon)
                    xlim([min(x) max(x)])
                    ylim([min(lo) max(hi)])
                    yl = ylim;
                    hold on
                    plot(x,ym,'k')
                    plot([etime etime],[yl(1) yl(2)],':k')
                    if ii==1
                        text(min(x)+0.2*(max(x)-min(x)),yl(2)-0.1*diff(yl),['sample ' num2str(ie)])
                        plot([min(x)+0.11*(max(x)-min(x)), min(x)+0.19*(max(x)-min(x))], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], ':k')
                    end
                    hold off
                    if ismember(index(ii),index(end,:)), xlabel('year-day'); else, xlabel([]); end
                    if ismember(index(ii),index(:,1)), ylabel('mmol C m^{-3}'); else, ylabel([]); end
                    title([num2str(round(fixedParams.PPdia(ii),2,'significant')) ' \mum'])
                end
                
                %----------------------------------------------------------
                
            case 'phytoplanktonStacked'
                
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 12 4];
                
                Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));
                Y = squeeze(sum(repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y, 2));
                % polygon vertices
                xpgon = [x(:); flip(x(:))];
                ypgon = zeros(nPP_size + 1, nt);
                ypgon(2:end,:) = mean(Y,3);
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
                Y = squeeze(sum(repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y, 2));
                Y = cat(1, Y, sum(fixedParams.zwidth .* Z));
                
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
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
    case 'barplot_TimeSeries'
        var = varargin{1}; % variables to plot
        ie = varargin{2}; % event number
        itraj = varargin{3}; % trajectories
        out = varargin{4}; % model output
        fixedParams = varargin{5};
        forcing = varargin{6};
        auxVars = varargin{7};
        dat = varargin{8};
        alpha = varargin{9};

        nt = fixedParams.nt;
        nz = fixedParams.nz;
        nPP_size = fixedParams.nPP_size;
        
        % Extract outputs
        N = squeeze(out.N(:,:,:,itraj));
        P = squeeze(out.P(:,:,:,:,itraj));
        Z = squeeze(out.Z(:,:,:,itraj));
        OM = squeeze(out.OM(:,:,:,:,itraj));

        switch var
            case 'phytoZooPlankton'                
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 8 6];
                % abundance summed over depth                
                Y = P(:,:,fixedParams.PP_C_index,:,:);                
                Y = squeeze(sum(repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y, 2));
                Y = cat(1, Y, sum(fixedParams.zwidth .* Z));
                mY = median(Y, 3);
                % color scale - red --> blue
                rb_hsv = rgb2hsv([1 0 0; 0 0 1]);
                r2b_hsv =  [linspace(rb_hsv(1,1), rb_hsv(2,1), nPP_size)' ...
                    linspace(rb_hsv(1,2), rb_hsv(2,2), nPP_size)' ...
                    linspace(rb_hsv(1,3), rb_hsv(2,3), nPP_size)'];
                r2b = hsv2rgb(r2b_hsv);
                cols = r2b;
                cols(nPP_size+1,:) = [0 0 0]; % black for zooplankton
                % group output into intervals of nd days
                nd = 50;
                int = 0:nd:ceil(nt/nd)*nd;
                intn = length(int)-1;
                y = nan(intn, nPP_size+1);
                for i = 1:intn
                    j0 = int(i)+1; j1 = min(int(i+1), nt);
                    iY = mY(:,j0:j1);
                    y(i,:) = mean(iY,2);
                end                
                bp = bar(y,'grouped');
                xl = xlim; yl = ylim;
                yleg = linspace(yl(1)+0.95*diff(yl), yl(1)+0.5*diff(yl), nPP_size+2);
                hold on
                text(xl(1)+0.8*diff(xl), yleg(1), 'cell diameter')
                for i = 1:nPP_size+1
                    bp(i).FaceColor = cols(i,:);
                    ps = scatter(xl(1)+0.8*diff(xl), yleg(i+1), 'filled');
                    ps.MarkerFaceColor = cols(i,:);
                    if i < nPP_size+1
                        text(xl(1)+0.82*diff(xl), yleg(i+1), ...
                            [num2str(round(fixedParams.PPdia(i),2,'significant')) ' \mum'])
                    else
                        text(xl(1)+0.82*diff(xl), yleg(i+1), 'zooplankton')
                    end
                end
                hold off
                xticklabels(strcat(string(int(1:end-1)), '-', string(int(2:end))))
                xlabel('year-day')
                ylabel('abundance (mmol C)')
                title(['planktonic carbon in ' num2str(fixedParams.Htot) 'm deep (1m^2) water column'])
        end
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
    
    case 'outputVsData_depth_boxplot'
        
        var = varargin{1}; % variables to plot
        out = varargin{2}; % model output
        Data = varargin{3};
        FixedParams = varargin{4};
        plotStd = varargin{5}; % index [0, 1, 2] controlling whether raw or standardised data are plotted
        
        switch var
            case 'DIN'
                ind = strcmp('N', Data.scalar.Variable); % index data
                v = 'N'; % label used in model output
                xlab = 'DIN (mmol N m^{-3})';
                xlab_std = 'standardised DIN';
                scaleFun = 'scaleFun_N';
            case 'PON'
                ind = strcmp('PON', Data.scalar.Variable); % index data
                v = 'OM'; % label used in model output
                xlab = 'PON (mmol N m^{-3})';
                xlab_std = 'standardised PON';
                scaleFun = 'scaleFun_PON';
            case 'POC'
                ind = strcmp('POC', Data.scalar.Variable); % index data
                v = 'OM'; % label used in model output
                xlab = 'POC (mmol C m^{-3})';
                xlab_std = 'standardised POC';
                scaleFun = 'scaleFun_POC';
            case 'Chl'
                ind = strcmp('chl_a', Data.scalar.Variable); % index data
                v = 'P'; % label used in model output
                xlab = 'Chl_a (mg Chl_a m^{-3})';
                xlab_std = 'standardised Chl_a';
                scaleFun = 'scaleFun_chl_a';
        end
        
        figure
                
        % data
        event = Data.scalar.Event(ind);
        depth = Data.scalar.Depth(ind);
        value = Data.scalar.Value(ind);
        value_scaled = Data.scalar.scaled_Value(ind);
        
        uevent = unique(event);
        udepth = unique(depth);
        nevent = length(uevent);
        ndepth = length(udepth);
        
        % equivalent model output
        depthMod = [];
        valueMod = [];
        valueMod_scaled = [];
        time = Data.scalar.Yearday(ind);
        
        for i = 1:nevent
            ev = uevent(i);
            ind_ev = ind & Data.scalar.Event == ev;
            evTime = unique(time(event == ev));
            evDepths = depth(event == ev);
            traj = find(Data.scalar.EventTraj(ev,:));
            ntraj = length(traj);
            depthMod = [depthMod; repmat(evDepths, [ntraj, 1])];
            x = out.(v);
            switch var
                case 'DIN', x = squeeze(x(:,:,evTime,traj)); % modelled values [depth,trajectory]
                case 'PON', x = squeeze(x(FixedParams.POM_index,:,FixedParams.OM_N_index,evTime,traj));
                case 'POC', x = squeeze(x(FixedParams.POM_index,:,FixedParams.OM_C_index,evTime,traj));
                case 'Chl', x = squeeze(sum(x(:,:,FixedParams.PP_Chl_index,evTime,traj)));
            end
            x = interp1(abs(FixedParams.z), x, evDepths); % interpolate to match sampling depths
            valueMod = [valueMod; x(:)];
            x = Data.scalar.(scaleFun)(Data.scalar.scale_mu(ind_ev), Data.scalar.scale_sig(ind_ev), x);
%             x = Data.scalar.scaleFun_N(Data.scalar.scale_mu(ind_ev), Data.scalar.scale_sig(ind_ev), x);
            valueMod_scaled = [valueMod_scaled; x(:)];
        end
        for i = ndepth:-1:1
            x = value(depth == udepth(i));
            allValues(1:length(x),2*ndepth-2*i+1) = x;
            x = valueMod(depthMod == udepth(i));
            allValues(1:length(x),2*ndepth-2*i+2) = x;
            x = value_scaled(depth == udepth(i));
            allValues_scaled(1:length(x),2*ndepth-2*i+1) = x;
            x = valueMod_scaled(depthMod == udepth(i));
            allValues_scaled(1:length(x),2*ndepth-2*i+2) = x;
        end
        allValues(~isnan(allValues) & allValues == 0) = nan;
        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
        
        dataCol = [0 0 0];
        modCol = [0 1 0];
        
        if ismember(plotStd, [0 1])
            if plotStd == 0, bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7); end
            if plotStd == 1, bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7); end
            for i = 1:size(allValues,2)
                bp.out(i).Marker = '.';
            end
            for i = 1:size(allValues,2) / 2
                bp.out(2*i-1).MarkerFaceColor = dataCol;
                bp.out(2*i-1).MarkerEdgeColor = dataCol;
                bp.out(2*i).MarkerFaceColor = modCol;
                bp.out(2*i).MarkerEdgeColor = modCol;
                bp.box(2*i-1).Color = dataCol;
                bp.box(2*i).Color = modCol;
                bp.med(2*i-1).Color = [0 0 0];
                bp.med(2*i).Color = [0 0 0];
            end
            gc = gca;
            tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
            gc.YTick = mean(tt);
            gc.YTickLabel = num2str(flip(udepth));
            gc.TickLength = 0.5 * gc.TickLength;
            xl = gc.XLim; yl = gc.YLim;
            tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
            hold on
            for i = 1:length(tt)
                line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
            end
            hold off
            if plotStd == 0, xlabel(xlab); end
            if plotStd == 1, xlabel(xlab_std); end
            ylabel('Depth (m)')

        else
            
            j = 1;
            while j < 3
                subplot(1,2,j)
                if j == 1, bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7); end
                if j == 2, bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7); end
                for i = 1:size(allValues,2)
                    bp.out(i).Marker = '.';
                end
                for i = 1:size(allValues,2) / 2
                    bp.out(2*i-1).MarkerFaceColor = dataCol;
                    bp.out(2*i-1).MarkerEdgeColor = dataCol;
                    bp.out(2*i).MarkerFaceColor = modCol;
                    bp.out(2*i).MarkerEdgeColor = modCol;
                    bp.box(2*i-1).Color = dataCol;
                    bp.box(2*i).Color = modCol;
                    bp.med(2*i-1).Color = [0 0 0];
                    bp.med(2*i).Color = [0 0 0];
                end
                gc = gca;
                tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                gc.YTick = mean(tt);
                gc.YTickLabel = num2str(flip(udepth));
                gc.TickLength = 0.5 * gc.TickLength;
                xl = gc.XLim; yl = gc.YLim;
                tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                hold on
                for i = 1:length(tt)
                    line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                end
                hold off
                if j == 1, xlabel(xlab); end
                if j == 2, xlabel(xlab_std); end
                ylabel('Depth (m)')
                j = j + 1;
            end
        end
        % legend
        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', dataCol)
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', modCol)

        %------------------------------------------------------------------
        %------------------------------------------------------------------
    
    case 'outputVsData_event_boxplot'

        var = varargin{1}; % variables to plot
        out = varargin{2}; % model output
        Data = varargin{3};
        FixedParams = varargin{4};
        plotStd = varargin{5}; % index [0, 1, 2] controlling whether raw or standardised data are plotted
        
        switch var
            case 'DIN'
                ind = strcmp('N', Data.scalar.Variable); % index data
                v = 'N'; % label used in model output
                xlab = 'DIN (mmol N m^{-3})';
                xlab_std = 'standardised DIN';
                scaleFun = 'scaleFun_N';
            case 'PON'
                ind = strcmp('PON', Data.scalar.Variable); % index data
                v = 'OM'; % label used in model output
                xlab = 'PON (mmol N m^{-3})';
                xlab_std = 'standardised PON';
                scaleFun = 'scaleFun_PON';
            case 'POC'
                ind = strcmp('POC', Data.scalar.Variable); % index data
                v = 'OM'; % label used in model output
                xlab = 'POC (mmol C m^{-3})';
                xlab_std = 'standardised POC';
                scaleFun = 'scaleFun_POC';
            case 'Chl'
                ind = strcmp('chl_a', Data.scalar.Variable); % index data
                v = 'P'; % label used in model output
                xlab = 'Chl_a (mg Chl_a m^{-3})';
                xlab_std = 'standardised Chl_a';
                scaleFun = 'scaleFun_chl_a';
        end
        
        figure
                
        % data
        event = Data.scalar.Event(ind);
        depth = Data.scalar.Depth(ind);
        value = Data.scalar.Value(ind);
        value_scaled = Data.scalar.scaled_Value(ind);
        
        uevent = unique(event);
        udepth = unique(depth);
        nevent = length(uevent);
        ndepth = length(udepth);
        
        % equivalent model output
        depthMod = [];
        eventMod = [];
        valueMod = [];
        valueMod_scaled = [];
        time = Data.scalar.Yearday(ind);
        
        for i = 1:nevent
            ev = uevent(i);
            ind_ev = ind & Data.scalar.Event == ev;
            evTime = unique(time(event == ev));
            evDepths = depth(event == ev);
            traj = find(Data.scalar.EventTraj(ev,:));
            ntraj = length(traj);
            depthMod = [depthMod; repmat(evDepths, [ntraj, 1])];
            eventMod = [eventMod; repmat(ev, [length(evDepths)*ntraj 1])];
            x = out.(v);
            switch var
                case 'DIN', x = squeeze(x(:,:,evTime,traj)); % modelled values [depth,trajectory]
                case 'PON', x = squeeze(x(FixedParams.POM_index,:,FixedParams.OM_N_index,evTime,traj));
                case 'POC', x = squeeze(x(FixedParams.POM_index,:,FixedParams.OM_C_index,evTime,traj));
                case 'Chl', x = squeeze(sum(x(:,:,FixedParams.PP_Chl_index,evTime,traj)));
            end
            x = interp1(abs(FixedParams.z), x, evDepths); % interpolate to match sampling depths
            valueMod = [valueMod; x(:)];
            x = Data.scalar.(scaleFun)(Data.scalar.scale_mu(ind_ev), Data.scalar.scale_sig(ind_ev), x);
%             x = Data.scalar.scaleFun_N(Data.scalar.scale_mu(ind_ev), Data.scalar.scale_sig(ind_ev), x);
            valueMod_scaled = [valueMod_scaled; x(:)];
        end
        
        for i = 1:nevent
            ev = uevent(i);
            x = value(event == ev);
            allValues(1:length(x),2*i-1) = x;
            x = valueMod(eventMod == ev);
            allValues(1:length(x),2*i) = x;
            x = value_scaled(event == ev);
            allValues_scaled(1:length(x),2*i-1) = x;
            x = valueMod_scaled(eventMod == ev);
            allValues_scaled(1:length(x),2*i) = x;
        end
        allValues(~isnan(allValues) & allValues == 0) = nan;
        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
        
        dataCol = [0 0 0];
        modCol = [0 1 0];
        
        if ismember(plotStd, [0 1])
            if plotStd == 0, bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7); end
            if plotStd == 1, bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7); end
            for i = 1:size(allValues,2)
                bp.out(i).Marker = '.';
            end
            for i = 1:size(allValues,2) / 2
                bp.out(2*i-1).MarkerFaceColor = dataCol;
                bp.out(2*i-1).MarkerEdgeColor = dataCol;
                bp.out(2*i).MarkerFaceColor = modCol;
                bp.out(2*i).MarkerEdgeColor = modCol;
                bp.box(2*i-1).Color = dataCol;
                bp.box(2*i).Color = modCol;
                bp.med(2*i-1).Color = [0 0 0];
                bp.med(2*i).Color = [0 0 0];
            end
            gc = gca;
            tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
            gc.YTick = mean(tt);
            gc.YTickLabel = num2str(flip(udepth));
            gc.TickLength = 0.5 * gc.TickLength;
            xl = gc.XLim; yl = gc.YLim;
            tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
            hold on
            for i = 1:length(tt)
                line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
            end
            hold off
            if plotStd == 0, xlabel(xlab); end
            if plotStd == 1, xlabel(xlab_std); end
            ylabel('Depth (m)')

        else
            
            j = 1;
            while j < 3
                subplot(1,2,j)
                if j == 1, bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7); end
                if j == 2, bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7); end
                for i = 1:size(allValues,2)
                    bp.out(i).Marker = '.';
                end
                for i = 1:size(allValues,2) / 2
                    bp.out(2*i-1).MarkerFaceColor = dataCol;
                    bp.out(2*i-1).MarkerEdgeColor = dataCol;
                    bp.out(2*i).MarkerFaceColor = modCol;
                    bp.out(2*i).MarkerEdgeColor = modCol;
                    bp.box(2*i-1).Color = dataCol;
                    bp.box(2*i).Color = modCol;
                    bp.med(2*i-1).Color = [0 0 0];
                    bp.med(2*i).Color = [0 0 0];
                end
                gc = gca;
                tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                gc.YTick = mean(tt);
                gc.YTickLabel = num2str((uevent));
                gc.TickLength = 0.5 * gc.TickLength;
                xl = gc.XLim; yl = gc.YLim;
                tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                hold on
                for i = 1:length(tt)
                    line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                end
                hold off
                if j == 1, xlabel(xlab); end
                if j == 2, xlabel(xlab_std); end
                ylabel('Event')
                j = j + 1;
            end
        end
        % legend
        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', dataCol)
        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', modCol)
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
    

    case 'outputVsData_summaryPlots'
        
        v = varargin{1}; % variable to plot
        Data = varargin{2}; % model output
        modData = varargin{3};
        FixedParams = varargin{4};
        costfun = FixedParams.costFunction;
        
        % Different selection of plots for different cost functions
        switch costfun
           
            case {'polyLikelihood', 'polyLikelihood2'}
                switch v
                    % summary plots displaying model fit to scalar data
                    case {'N','PON','POC','chl_a'}
                        
                        coldat = [0 0 0];
                        colmod = [0 1 0];
                        
                        plt = figure;
                        plt.Units = 'inches';
                        plt.Position = [0 0 16 12];
                        
                        % Model fit to polynomial coefficients representing data
                        subplot(2,3,1)
                        
                        cdat = flip(Data.scalar.(['polyCoefs_' v]));
                        cmod = flip(modData.scalar.(['polyCoefs_' v]));
                        nc = length(cdat);
                        x = 1:nc;
                        mu_c = mean(cmod, 2);
                        cmin = min(cmod,[],2);
                        cmax = max(cmod,[],2);
                        
                        yl = [min([cmin(:);cdat(:)]), max([cmax(:);cdat(:)])];
                        yl(1) = floor(yl(1));
                        yl(2) = ceil(yl(2));
                        
                        scatter(x-0.125, cdat, 'MarkerEdgeColor', coldat, 'MarkerFaceColor', coldat)
                        hold on
                        scatter(x+0.125, mu_c, 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                        for i = 1:nc
                            line([x(i) x(i)]+0.125, [cmin(i) cmax(i)], 'Color', colmod)
                            if i < nc
                                line([x(i) x(i)]+0.5, yl, 'Color', [0.5 0.5 0.5], 'LineStyle', ':')
                            end
                        end
                        
                        gc = gca;
                        gc.XLim = [0.5, nc+0.5];
                        xl = gc.XLim;
                        line(xl, [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
                        
                        gc.XTick = 1:nc;
                        gc.XTickLabel{1} = 'intercept';
                        
                        for i = 1:nc
                            if i > 1, l = ['\beta_' num2str(i-1)];
                            else, l = ['\beta_' num2str(i-1) ' (intercept)'];
                            end
                            gc.XTickLabel{i} = l;
                        end
                        
                        xl = gc.XLim;
                        xl(2) = xl(2) + diff(xl) / 5;
                        gc.XLim = xl;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        switch v
                            case 'N'
                                Var = 'DIN';
                            case 'PON'
                                Var = 'PON';
                            case 'POC'
                                Var = 'POC';
                            case 'chl_a'
                                Var = 'chl-a';
                        end
                        title(['Polynomial coefficients representing ' Var  ' data'])
                        hold off
                        
                        
                        % Sorted data and polynomial curves
                        subplot(2,3,4)
                        
                        ind = strcmp(Data.scalar.Variable, v);
                        x = Data.scalar.(['polyXvals_' v]);
                        o = Data.scalar.(['sortOrder_' v]);
                        ymod = modData.scalar.scaled_Value(ind,:);
                        ydat = Data.scalar.scaled_Value(ind);
                        ymod = ymod(o,:);
                        ydat = ydat(o);
                        
                        for i = 1:size(ymod, 2)
                            %     scatter(x, ymod(:,i), 'Marker', '.', 'MarkerEdgeColor', [0.7 0.7 0.7]);
                            scatter(x, ymod(:,i), 'Marker', '.', 'MarkerEdgeColor', colmod);
                            if i == 1, hold on; end
                            modpoly = polyval(modData.scalar.(['polyCoefs_' v])(:,i), x);
                            plot(x, modpoly, 'Color', colmod)
                        end
                        scatter(x, ydat, 'MarkerEdgeColor', coldat)
                        datpoly = polyval(Data.scalar.(['polyCoefs_' v]), x);
                        plot(x, datpoly, 'Color', coldat)
                        hold off
                        
                        gc = gca;
                        xl = gc.XLim; yl = gc.YLim;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        ylabel(['sorted standardised ' Var])
                        title(['Polynomials representing ' Var ' data'])
                        
                        
                        
                        % depth vs standardised variable boxplot
                        subplot(2,3,2)
                        
                        switch v
                            case 'N'
                                xlab = 'DIN (mmol N m^{-3})';
                                xlab_std = 'standardised DIN';
                            case 'PON'
                                xlab = 'PON (mmol N m^{-3})';
                                xlab_std = 'standardised PON';
                            case 'POC'
                                xlab = 'POC (mmol C m^{-3})';
                                xlab_std = 'standardised POC';
                            case 'chl_a'
                                xlab = 'chl-a (mg chl-a m^{-3})';
                                xlab_std = 'standardised chl-a';
                        end
                        
                        ind = strcmp(Data.scalar.Variable, v);
                        
                        event = Data.scalar.Event(ind);
                        depth = Data.scalar.Depth(ind);
                        value = Data.scalar.Value(ind);
                        value_scaled = Data.scalar.scaled_Value(ind);
                        
                        uevent = unique(event);
                        udepth = unique(depth);
                        nevent = length(uevent);
                        ndepth = length(udepth);
                        
                        valueMod = modData.scalar.Value(ind,:);
                        valueMod_scaled = modData.scalar.scaled_Value(ind,:);
                        
                        clear allValues allValues_scaled
                        for i = ndepth:-1:1
                            di = depth == udepth(i);
                            x = value(di);
                            allValues(1:numel(x),2*ndepth-2*i+1) = x(:);
                            x = valueMod(di,:);
                            allValues(1:numel(x),2*ndepth-2*i+2) = x(:);
                            x = value_scaled(di);
                            allValues_scaled(1:numel(x),2*ndepth-2*i+1) = x(:);
                            x = valueMod_scaled(di,:);
                            allValues_scaled(1:numel(x),2*ndepth-2*i+2) = x(:);
                        end
                        allValues(~isnan(allValues) & allValues == 0) = nan;
                        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
                        
                        bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(flip(udepth));
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        gc.XLim = xl;
                        
                        xlabel(xlab_std)
                        ylabel('depth (m)')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        
                        % depth vs raw variable boxplot
                        subplot(2,3,5)
                        
                        bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(flip(udepth));
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        xlabel(xlab)
                        ylabel('depth (m)')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        % event vs standardised variable boxplot
                        subplot(2,3,3)
                        
                        clear allValues allValues_scaled
                        for i = 1:nevent
                            ei = event == uevent(i);
                            x = value(ei);
                            allValues(1:numel(x),2*i-1) = x(:);
                            x = valueMod(ei,:);
                            allValues(1:numel(x),2*i) = x(:);
                            x = value_scaled(ei);
                            allValues_scaled(1:numel(x),2*i-1) = x(:);
                            x = valueMod_scaled(ei);
                            allValues_scaled(1:numel(x),2*i) = x(:);
                        end
                        allValues(~isnan(allValues) & allValues == 0) = nan;
                        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
                        
                        bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(uevent);
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        xlabel(xlab_std)
                        ylabel('sampling event')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        % event vs raw variable boxplot
                        subplot(2,3,6)
                        
                        bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(uevent);
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        xlabel(xlab)
                        ylabel('sampling event')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        % summary plots displaying model fit to size spectra data
                    case {'NConc', 'CellConc'}
                        
%                         Var = v;
                        
                        coldat = [0 0 0];
                        colmod = [0 1 0];
                        
                        plt = figure;
                        plt.Units = 'inches';
                        plt.Position = [0 0 12 12];
                        
                        % Model fit to polynomial coefficients representing data
                        subplot(2,2,1)
                        
                        cdat = flip(Data.size.(['polyCoefs_' v]));
                        cmod = flip(modData.size.(['polyCoefs_' v]));
                        nc = length(cdat);
                        x = 1:nc;
                        mu_c = mean(cmod, 2);
                        cmin = min(cmod,[],2);
                        cmax = max(cmod,[],2);
                        
                        yl = [min([cmin(:);cdat(:)]), max([cmax(:);cdat(:)])];
                        yl(1) = floor(yl(1));
                        yl(2) = ceil(yl(2));
                        
                        scatter(x-0.125, cdat, 'MarkerEdgeColor', coldat, 'MarkerFaceColor', coldat)
                        hold on
                        scatter(x+0.125, mu_c, 'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod)
                        for i = 1:nc
                            line([x(i) x(i)]+0.125, [cmin(i) cmax(i)], 'Color', colmod)
                            if i < nc
                                line([x(i) x(i)]+0.5, yl, 'Color', [0.5 0.5 0.5], 'LineStyle', ':')
                            end
                        end
                        
                        gc = gca;
                        gc.XLim = [0.5, nc+0.5];
                        xl = gc.XLim;
                        line(xl, [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
                        
                        gc.XTick = 1:nc;
                        gc.XTickLabel{1} = 'intercept';
                        
                        for i = 1:nc
                            if i > 1, l = ['\beta_' num2str(i-1)];
                            else, l = ['\beta_' num2str(i-1) ' (intercept)'];
                            end
                            gc.XTickLabel{i} = l;
                        end
                        
                        xl = gc.XLim;
                        xl(2) = xl(2) + diff(xl) / 5;
                        gc.XLim = xl;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        title('Polynomial coefficients representing size spectra data')
                        hold off
                        
                        
                        % Sorted data and polynomial curves
                        subplot(2,2,3)
                        x = Data.size.(['polyXvals_' v]);
                        o = Data.size.(['sortOrder_' v]);
                        
                        ymod = modData.size.scaled_Value(strcmp(modData.size.Variable, v));
                        ydat = Data.size.dataBinned.scaled_Value( ...
                            strcmp(Data.size.dataBinned.Variable, v));
%                         ydat = Data.size.dataBinned.(['scaled_' Var]);
                        ymod = ymod(o,:);
                        ydat = ydat(o);
                        
                        for i = 1:size(ymod, 2)
                            scatter(x, ymod(:,i), 'Marker', '.', 'MarkerEdgeColor', colmod);
                            if i == 1, hold on; end
%                             modpoly = polyval(modData.size.polyCoefs_Nconc(:,i), x);
                            modpoly = polyval(modData.size.(['polyCoefs_' v])(:,i), x);
                            plot(x, modpoly, 'Color', colmod)
                        end
                        scatter(x, ydat, 'MarkerEdgeColor', coldat)
                        datpoly = polyval(Data.size.(['polyCoefs_' v]), x);
                        plot(x, datpoly, 'Color', coldat)
                        hold off
                        
                        gc = gca;
                        xl = gc.XLim; yl = gc.YLim;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        switch v
                            case 'Nconc'
                                ylabel('sorted standardised N-at-size spectra')
                            case 'cellTot'
                                ylabel('sorted standardised cell abundance-at-size')
                        end
                        
                        title('Polynomials representing size spectra data')
                        
                        
                        
                        % standardised variable vs size
                        subplot(2,2,2)
                        ind = strcmp(Data.size.dataBinned.Variable, v);
%                         Size = Data.size.dataBinned.size;
                        Size = Data.size.dataBinned.size(ind);
%                         SizeClass = Data.size.dataBinned.sizeClass;
                        SizeClass = Data.size.dataBinned.sizeClass(ind);
%                         value = Data.size.dataBinned.Ntot;
                        value = Data.size.dataBinned.Value(ind);
%                         value_scaled = Data.size.dataBinned.scaled_Ntot;
                        value_scaled = Data.size.dataBinned.scaled_Value(ind);
                        
                        usize = unique(Size);
                        usizeClass = unique(SizeClass);
                        nsize = length(usize);
                        
                        valueMod = modData.size.Value(strcmp(modData.size.Variable, v));
                        valueMod_scaled = modData.size.scaled_Value(strcmp(modData.size.Variable, v));
%                         valueMod_scaled = modData.size.(['scaled_' v]);
                        
                        semilogx(Size, value_scaled, '-o', ...
                            'MarkerEdgeColor', coldat, 'MarkerFaceColor', coldat, ...
                            'Color', coldat);
                        hold on
                        for i = 1:size(valueMod_scaled, 2)
                            semilogx(Size(:), valueMod_scaled(:,i), '-o', ...
                                'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod, ...
                                'Color', colmod)
                        end
                        gc = gca;
                        xl = gc.XLim;
                        if xl(1) == Size(1), xl(1) = 0.5 * xl(1); end
                        xl(2) = ceil(xl(2) / 100) * 100;
                        gc.XLim = xl;
                        yl = gc.YLim;
                        
                        %                 legend('data', 'model')
                        
                        xlabel('cell diameter (\mum)')
                        switch v
                            case 'Nconc'
                                ylabel('standardised N size spectra')
                            case 'cellTot'
                                ylabel('standardised cell abundance-at-size')
                        end
                        
                        % legend
                        lxl = log10(xl);
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.05*diff(yl), 'data')
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.1*diff(yl), 'model')
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        hold off
                        
                        
                        % raw variable vs size
                        subplot(2,2,4)
                        
                        semilogx(Size, value, '-o', ...
                            'MarkerEdgeColor', coldat, 'MarkerFaceColor', coldat, ...
                            'Color', coldat);
                        hold on
                        for i = 1:size(valueMod, 2)
                            semilogx(Size(:), valueMod(:,i), '-o', ...
                                'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod, ...
                                'Color', colmod)
                        end
                        gc = gca;
                        xl = gc.XLim;
                        if xl(1) == Size(1), xl(1) = 0.5 * xl(1); end
                        xl(2) = ceil(xl(2) / 100) * 100;
                        gc.XLim = xl;
                        yl = gc.YLim;
                        
                        %                 legend('data', 'model')
                        
                        xlabel('cell diameter (\mum)')
                        ylabel('N size spectra')
                        
                        % legend
                        lxl = log10(xl);
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.05*diff(yl), 'data')
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.1*diff(yl), 'model')
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        hold off
                        
                        
                end
                
            case 'LSS'
                switch v
                    % summary plots displaying model fit to scalar data
                    case {'N','PON','POC','chl_a'}
                        
                        coldat = [0 0 0];
                        colmod = [0 1 0];
                        
                        plt = figure;
                        plt.Units = 'inches';
                        plt.Position = [0 0 12 12];
                        
                        switch v
                            case 'N'
                                Var = 'DIN';
                                xlab = 'DIN (mmol N m^{-3})';
                                xlab_std = 'standardised DIN';
                            case 'PON'
                                Var = 'PON';
                                xlab = 'PON (mmol N m^{-3})';
                                xlab_std = 'standardised PON';
                            case 'POC'
                                Var = 'POC';
                                xlab = 'POC (mmol C m^{-3})';
                                xlab_std = 'standardised POC';
                            case 'chl_a'
                                Var = 'chl-a';
                                xlab = 'chl-a (mg chl-a m^{-3})';
                                xlab_std = 'standardised chl-a';
                        end

                        % Sorted standardised data
                        subplot(2,3,1)
                        
                        ind = strcmp(Data.scalar.Variable, v);
                        x = Data.scalar.(['polyXvals_' v]);
                        o = Data.scalar.(['sortOrder_' v]);
                        ymod = modData.scalar.scaled_Value(ind,:);
                        ydat = Data.scalar.scaled_Value(ind);
                        ymod = ymod(o,:);
                        ydat = ydat(o);
                        x2d = repmat(x, [1 size(ymod, 2)]);
                        
                        scatter(x2d(:), ymod(:), 'Marker', '.', 'MarkerEdgeColor', colmod)
                        hold on
                        scatter(x, ydat, 'MarkerEdgeColor', coldat)
                        hold off
                        
                        gc = gca;
                        xl = gc.XLim; yl = gc.YLim;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        ylabel(['standardised ' Var])
                        
                        
                        % Sorted raw data
                        subplot(2,3,4)
                        
                        ymod = modData.scalar.Value(ind,:);
                        ydat = Data.scalar.Value(ind);
                        [ydat,J] = sort(ydat);
                        ymod = ymod(J,:);
%                         ymod = ymod(o,:);
%                         ydat = ydat(o);
                        
                        scatter(x2d(:), ymod(:), 'Marker', '.', 'MarkerEdgeColor', colmod)
                        hold on
                        scatter(x, ydat, 'MarkerEdgeColor', coldat)
                        hold off
                        
                        gc = gca;
                        xl = gc.XLim; yl = gc.YLim;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        ylabel(xlab)
                        
                        
                        % depth vs standardised variable boxplot
                        subplot(2,3,2)
                                                
                        ind = strcmp(Data.scalar.Variable, v);
                        
                        event = Data.scalar.Event(ind);
                        depth = Data.scalar.Depth(ind);
                        value = Data.scalar.Value(ind);
                        value_scaled = Data.scalar.scaled_Value(ind);
                        
                        uevent = unique(event);
                        udepth = unique(depth);
                        nevent = length(uevent);
                        ndepth = length(udepth);
                        
                        valueMod = modData.scalar.Value(ind,:);
                        valueMod_scaled = modData.scalar.scaled_Value(ind,:);
                        
                        clear allValues allValues_scaled
                        for i = ndepth:-1:1
                            di = depth == udepth(i);
                            x = value(di);
                            allValues(1:numel(x),2*ndepth-2*i+1) = x(:);
                            x = valueMod(di,:);
                            allValues(1:numel(x),2*ndepth-2*i+2) = x(:);
                            x = value_scaled(di);
                            allValues_scaled(1:numel(x),2*ndepth-2*i+1) = x(:);
                            x = valueMod_scaled(di,:);
                            allValues_scaled(1:numel(x),2*ndepth-2*i+2) = x(:);
                        end
                        allValues(~isnan(allValues) & allValues == 0) = nan;
                        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
                        
                        bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(flip(udepth));
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        gc.XLim = xl;
                        
                        xlabel(xlab_std)
                        ylabel('depth (m)')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        
                        % depth vs raw variable boxplot
                        subplot(2,3,5)
                        
                        bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(flip(udepth));
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        xlabel(xlab)
                        ylabel('depth (m)')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        % event vs standardised variable boxplot
                        subplot(2,3,3)
                        
                        clear allValues allValues_scaled
                        for i = 1:nevent
                            ei = event == uevent(i);
                            x = value(ei);
                            allValues(1:numel(x),2*i-1) = x(:);
                            x = valueMod(ei,:);
                            allValues(1:numel(x),2*i) = x(:);
                            x = value_scaled(ei);
                            allValues_scaled(1:numel(x),2*i-1) = x(:);
                            x = valueMod_scaled(ei);
                            allValues_scaled(1:numel(x),2*i) = x(:);
                        end
                        allValues(~isnan(allValues) & allValues == 0) = nan;
                        allValues_scaled(~isnan(allValues_scaled) & allValues_scaled == 0) = nan;
                        
                        bp = boxplot2(allValues_scaled, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(uevent);
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        xlabel(xlab_std)
                        ylabel('sampling event')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        
                        % event vs raw variable boxplot
                        subplot(2,3,6)
                        
                        bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7);
                        
                        for i = 1:size(allValues,2)
                            bp.out(i).Marker = '.';
                        end
                        for i = 1:size(allValues,2) / 2
                            bp.out(2*i-1).MarkerFaceColor = coldat;
                            bp.out(2*i-1).MarkerEdgeColor = coldat;
                            bp.out(2*i).MarkerFaceColor = colmod;
                            bp.out(2*i).MarkerEdgeColor = colmod;
                            bp.box(2*i-1).Color = coldat;
                            bp.box(2*i).Color = colmod;
                            bp.med(2*i-1).Color = [0 0 0];
                            bp.med(2*i).Color = [0 0 0];
                        end
                        
                        gc = gca;
                        tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
                        gc.YTick = mean(tt);
                        gc.YTickLabel = num2str(uevent);
                        gc.TickLength = 0.5 * gc.TickLength;
                        xl = gc.XLim; yl = gc.YLim;
                        tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));
                        hold on
                        for i = 1:length(tt)
                            line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
                        end
                        hold off
                        
                        xlabel(xlab)
                        ylabel('sampling event')
                        
                        % legend
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        % summary plots displaying model fit to size spectra data
                    case 'N_at_size'
                        
                        Var = 'Ntot';
                        
                        coldat = [0 0 0];
                        colmod = [0 1 0];
                        
                        plt = figure;
                        plt.Units = 'inches';
                        plt.Position = [0 0 12 12];
                        
                        % Model fit to polynomial coefficients representing data
                        subplot(2,2,1)
                        
                        % Sorted standardised data
                        subplot(2,2,1)
                        x = Data.size.polyXvals;
                        o = Data.size.sortOrder;
                        ymod = modData.size.(['scaled_' Var]);
                        ydat = Data.size.dataBinned.(['scaled_' Var]);
                        ymod = ymod(o,:);
                        ydat = ydat(o);
                        x2d = repmat(x, [1 size(ymod, 2)]);
                        
                        scatter(x2d(:), ymod(:), 'Marker', '.', 'MarkerEdgeColor', colmod);
                        hold on
                        scatter(x, ydat, 'MarkerEdgeColor', coldat)
                        hold off
                        
                        gc = gca;
                        xl = gc.XLim; yl = gc.YLim;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        ylabel('standardised N size spectra')
                        
                        
                        % Sorted raw data
                        subplot(2,2,3)
                        ymod = modData.size.(Var);
                        ydat = Data.size.dataBinned.(Var);
                        [~,J] = sort(ydat);
                        ymod = ymod(J,:);
                        ydat = ydat(J);
                        scatter(x2d(:), ymod(:), 'Marker', '.', 'MarkerEdgeColor', colmod);
                        hold on
                        scatter(x, ydat, 'MarkerEdgeColor', coldat)
                        
                        gc = gca;
                        xl = gc.XLim; yl = gc.YLim;
                        text(xl(2)-0.15*diff(xl), yl(2)-0.05*diff(yl), 'data')
                        text(xl(2)-0.15*diff(xl), yl(2)-0.1*diff(yl), 'model')
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([xl(2)-0.22*diff(xl), xl(2)-0.17*diff(xl)], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        ylabel('size spectra (mmol N m^{-3})')
                        
                        
                        % standardised variable vs size
                        subplot(2,2,2)
                        
                        Size = Data.size.dataBinned.size;
                        SizeClass = Data.size.dataBinned.sizeClass;
                        value = Data.size.dataBinned.Ntot;
                        value_scaled = Data.size.dataBinned.scaled_Ntot;
                        
                        usize = unique(Size);
                        
                        valueMod = modData.size.(Var);
                        valueMod_scaled = modData.size.(['scaled_' Var]);
                        
                        semilogx(Size, value_scaled, '-o', ...
                            'MarkerEdgeColor', coldat, 'MarkerFaceColor', coldat, ...
                            'Color', coldat);
                        hold on
                        for i = 1:size(valueMod_scaled, 2)
                            semilogx(Size(:), valueMod_scaled(:,i), '-o', ...
                                'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod, ...
                                'Color', colmod)
                        end
                        gc = gca;
                        xl = gc.XLim;
                        if xl(1) == Size(1), xl(1) = 0.5 * xl(1); end
                        xl(2) = ceil(xl(2) / 100) * 100;
                        gc.XLim = xl;
                        yl = gc.YLim;
                        
                        xlabel('cell diameter (\mum)')
                        ylabel('standardised N size spectra')
                        
                        % legend
                        lxl = log10(xl);
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.05*diff(yl), 'data')
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.1*diff(yl), 'model')
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        hold off
                        
                        
                        % raw variable vs size
                        subplot(2,2,4)
                        
                        semilogx(Size, value, '-o', ...
                            'MarkerEdgeColor', coldat, 'MarkerFaceColor', coldat, ...
                            'Color', coldat);
                        hold on
                        for i = 1:size(valueMod, 2)
                            semilogx(Size(:), valueMod(:,i), '-o', ...
                                'MarkerEdgeColor', colmod, 'MarkerFaceColor', colmod, ...
                                'Color', colmod)
                        end
                        gc = gca;
                        xl = gc.XLim;
                        if xl(1) == Size(1), xl(1) = 0.5 * xl(1); end
                        xl(2) = ceil(xl(2) / 100) * 100;
                        gc.XLim = xl;
                        yl = gc.YLim;
                        
                        %                 legend('data', 'model')
                        
                        xlabel('cell diameter (\mum)')
                        ylabel('size spectra (mmol N m^{-3})')

                        % legend
                        lxl = log10(xl);
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.05*diff(yl), 'data')
                        text(10^(lxl(2)-0.15*diff(lxl)), yl(2)-0.1*diff(yl), 'model')
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.05*diff(yl), [1 2]), 'Color', coldat)
                        line([10^(lxl(2)-0.22*diff(lxl)), 10^(lxl(2)-0.17*diff(lxl))], repmat(yl(2)-0.1*diff(yl), [1 2]), 'Color', colmod)
                        
                        hold off
                        
                end
                
        end
        
end

end





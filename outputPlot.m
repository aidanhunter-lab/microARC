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
                if ~mod(fixedParams.nPP_size,2) % even number of size classes
                    index = reshape(1:fixedParams.nPP_size, [nc nr])';
                else
                    index = reshape([1:fixedParams.nPP_size 0], [nc nr])';
                end                
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
                suptitle('phytoplankton nitrogen concentration given cell diameter')
%                 suptitle('phytoplankton abundance given cell volume (mmol N / m^3)')
                colormap plasma
                
                %------------------------------------------------------------------
            
            case 'phytoplankton_Chl'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                if ~mod(fixedParams.nPP_size,2) % even number of size classes
                    index = reshape(1:fixedParams.nPP_size, [nc nr])';
                else
                    index = reshape([1:fixedParams.nPP_size 0], [nc nr])';
                end
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
                suptitle('phytoplankton chlorophyll_a concentration given cell diameter')
%                 suptitle('phytoplankton abundance given cell volume (mmol N / m^3)')
                colormap plasma
                
                %------------------------------------------------------------------
                
            case 'phytoplankton_C'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                if ~mod(fixedParams.nPP_size,2) % even number of size classes
                    index = reshape(1:fixedParams.nPP_size, [nc nr])';
                else
                    index = reshape([1:fixedParams.nPP_size 0], [nc nr])';
                end
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
                suptitle('phytoplankton carbon concentration given cell diameter')
                colormap plasma
                                
                %------------------------------------------------------------------

            case 'phytoplankton_N_C'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                if ~mod(fixedParams.nPP_size,2) % even number of size classes
                    index = reshape(1:fixedParams.nPP_size, [nc nr])';
                else
                    index = reshape([1:fixedParams.nPP_size 0], [nc nr])';
                end
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
                suptitle('phytoplankton N/C ratio given cell diameter')
                colormap plasma
                
                %------------------------------------------------------------------
            
            case 'phytoplankton_Chl_N'
                plt = figure;
                nc = floor(sqrt(fixedParams.nPP_size));
                nr = ceil(fixedParams.nPP_size / nc);
                plt.Units = 'inches';
                plt.Position = [0 0 8*nc 3*nr];
                if ~mod(fixedParams.nPP_size,2) % even number of size classes
                    index = reshape(1:fixedParams.nPP_size, [nc nr])';
                else
                    index = reshape([1:fixedParams.nPP_size 0], [nc nr])';
                end
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
                suptitle('phytoplankton Chl/N ratio given cell diameter')
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
                suptitle('Forcing data - surface layer')
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
                if ~mod(fixedParams.nPP_size, 2)
                    index = reshape(1:fixedParams.nPP_size, [nc nr])';
                else
                    index = reshape([1:fixedParams.nPP_size 0], [nc nr])';
                end
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
                suptitle('phytoplankton abundance / cell diameter')
%                 suptitle('phytoplankton abundance / cell volume')
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

                suptitle('organic matter')
                

                %----------------------------------------------------------
                
            case 'phytoplankton_C'
                plt = figure;
                plt.Units = 'inches';
                plt.Position = [0 0 12 12];                
                Y = squeeze(P(:,:,fixedParams.PP_C_index,:,:));                
                Y = squeeze(max(Y,[],2));
                nc = floor(sqrt(nPP_size));
                nr = ceil(nPP_size / nc);
                if ~mod(nPP_size, 2)
                    index = reshape(1:nPP_size, [nc nr])';
                else
                    index = reshape([1:nPP_size 0], [nc nr])';
                end
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
        
end





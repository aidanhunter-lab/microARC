function plt = plot_contour_DepthTime(var, traj, out, auxVars, ...
    fixedParams, forcing, smooth, varargin)

if ~isempty(varargin)
    if any(strcmp(varargin, 'Event'))
        sampleEvent = varargin{find(strcmp(varargin, 'Event'))+1};
    end
    if any(strcmp(varargin, 'waterOrigin'))
        waterOrigin = varargin{find(strcmp(varargin, 'waterOrigin'))+1};
    end
end

fixedParams.nPP_size = double(fixedParams.nPP_size);
fixedParams.nZP_size = double(fixedParams.nZP_size);


% Extract outputs
N = squeeze(out.N(:,:,:,traj));
P = squeeze(out.P(:,:,:,:,traj));
Z = squeeze(out.Z(:,:,:,:,traj));
OM = squeeze(out.OM(:,:,:,:,traj));

% If multiple trajectories selected then take averages
ntraj = length(traj);
if ntraj > 1
    N = mean(N, ndims(N));
    P = mean(P, ndims(P));
    Z = mean(Z, ndims(Z));
    OM = mean(OM, ndims(OM));
end

% Time-depth grid for interpolation
[depth, time] = ndgrid(abs(fixedParams.z), 1:fixedParams.nt);
[depthGrid, timeGrid] = ndgrid(1:1:abs(fixedParams.zw(end)), 1:fixedParams.nt);

switch var
    
    case 'forcing'
        
        plt = figure;
        plt.Units = 'inches';
        plt.Position = [0 0 8 9];
        
        % Temperaturele
        subplot(3,1,1)
        x = forcing.T(:,:,traj);
        if ntraj > 1
            x = mean(x, ndims(x));
        end
        F = griddedInterpolant(depth, time, x, smooth);
        Fsmooth = flip(F(depthGrid, timeGrid));
        contourf(Fsmooth)
        cb = colorbar;
        cb.Label.String = '\circC';
        title('Temperature')
        ylabel('depth (m)')
        xticks(100:100:fixedParams.nt)
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(linspace(fixedParams.zw(end),0,7))
        
        % Diffusivity
        subplot(3,1,2);
        x = forcing.K_center(:,:,traj);
        if ntraj > 1
            x = mean(x, ndims(x));
        end
        F = griddedInterpolant(depth, time, x, smooth);
        Fsmooth = flip(F(depthGrid, timeGrid));
        Fsmooth(Fsmooth<=0) = nan;
        contourf(log10(Fsmooth))
        cb = colorbar;
        for ii = 1:length(cb.TickLabels), cb.TickLabels{ii} = string(round(10 ^ str2double(cb.TickLabels{ii}),2,'significant')); end
        cb.Label.String = 'm^2 day^{-1}';
        title('Diffusivity')
        ylabel('depth (m)')
        xticks(100:100:fixedParams.nt)
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(linspace(fixedParams.zw(end),0,7))
        %                 set(gca,'colorscale','log')
        
        % PAR
        subplot(3,1,3)
        x = squeeze(auxVars.I(:,:,:,traj));
        if ntraj > 1
            x = mean(x, ndims(x));
        end
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
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(linspace(fixedParams.zw(end),0,7))
        colormap plasma
        sgtitle(['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'])
        
        %----------------------------------------------------------
    
    case 'DIN'
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
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(linspace(fixedParams.zw(end),0,7))
        colormap plasma
        title(['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'])
        subtitle('DIN')

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
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
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
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
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
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
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
        xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
        yticks(linspace(0,abs(fixedParams.zw(end)),7))
        yticklabels(linspace(fixedParams.zw(end),0,7))
        colormap plasma
        sgtitle(['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'])

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
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'phytoplankton nitrogen concentration given cell diameter'})
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
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'phytoplankton chlorophyll_a concentration given cell diameter'})
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
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'phytoplankton carbon concentration given cell diameter'})
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
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'phytoplankton N/C ratio given cell diameter'})
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
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'phytoplankton Chl/N ratio given cell diameter'})
        colormap plasma
        
        %------------------------------------------------------------------
        

    
    
    
    
    case 'zooplankton_C'
        plt = figure;        
        nc = floor(sqrt(fixedParams.nZP_size));
        nr = ceil(fixedParams.nZP_size / nc);
        plt.Units = 'inches';
        plt.Position = [0 0 8*nc 3*nr];
        index = reshape([1:fixedParams.nZP_size, ...
            zeros(1, nc * nr - fixedParams.nZP_size)], [nc nr])';
        for ii = 1:fixedParams.nZP_size
            subplot(nr,nc,index(ii))
            x = squeeze(Z(ii,:,fixedParams.ZP_C_index,:));
            F = griddedInterpolant(depth, time, x, smooth);
            Fsmooth = flip(F(depthGrid, timeGrid));
            contourf(Fsmooth)
            cb = colorbar;
            cb.Label.String = 'mmol C / m^3';
            title([num2str(round(fixedParams.ZPdia(ii),2,'significant')) ' \mum'])
            xlabel('year-day')
            ylabel('depth (m)')
            xticks(100:100:fixedParams.nt)
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'zooplankton carbon concentration given cell diameter'})
        colormap plasma

    case 'zooplankton_N'
        plt = figure;
        nc = floor(sqrt(fixedParams.nZP_size));
        nr = ceil(fixedParams.nZP_size / nc);
        plt.Units = 'inches';
        plt.Position = [0 0 8*nc 3*nr];
        index = reshape([1:fixedParams.nZP_size, ...
            zeros(1, nc * nr - fixedParams.nZP_size)], [nc nr])';
        for ii = 1:fixedParams.nZP_size
            subplot(nr,nc,index(ii))
            x = squeeze(Z(ii,:,fixedParams.ZP_N_index,:));
            F = griddedInterpolant(depth, time, x, smooth);
            Fsmooth = flip(F(depthGrid, timeGrid));
            contourf(Fsmooth)
            cb = colorbar;
            cb.Label.String = 'mmol N / m^3';
            title([num2str(round(fixedParams.ZPdia(ii),2,'significant')) ' \mum'])
            xlabel('year-day')
            ylabel('depth (m)')
            xticks(100:100:fixedParams.nt)
            xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
            yticks(linspace(0,abs(fixedParams.zw(end)),7))
            yticklabels(linspace(fixedParams.zw(end),0,7))
        end
        sgtitle({['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'], ...
            'zooplankton nitrogen concentration given cell diameter'})
        colormap plasma

%     case 'zooplankton_C'
%         plt = figure;
%         plt.Units = 'inches';
%         plt.Position = [0 0 8 3];
%         x = squeeze(Z(:,fixedParams.ZP_C_index,:));
%         F = griddedInterpolant(depth, time, x, smooth);
%         Fsmooth = flip(F(depthGrid, timeGrid));
%         contourf(Fsmooth)
%         cb = colorbar;
%         cb.Label.String = 'mmol C / m^3';
%         title('zooplankton abundance')
%         xlabel('year-day')
%         ylabel('depth (m)')
%         xticks(100:100:fixedParams.nt)
%         xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
%         yticks(linspace(0,abs(fixedParams.zw(end)),7))
%         yticklabels(linspace(fixedParams.zw(end),0,7))
%         title(['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'])
%         subtitle('zooplankton abundance')
%         colormap plasma
%         
%     case 'zooplankton_N'
%         plt = figure;
%         plt.Units = 'inches';
%         plt.Position = [0 0 8 3];
%         x = squeeze(Z(:,fixedParams.ZP_N_index,:));
%         F = griddedInterpolant(depth, time, x, smooth);
%         Fsmooth = flip(F(depthGrid, timeGrid));
%         contourf(Fsmooth)
%         cb = colorbar;
%         cb.Label.String = 'mmol N / m^3';
%         title('zooplankton abundance')
%         xlabel('year-day')
%         ylabel('depth (m)')
%         xticks(100:100:fixedParams.nt)
%         xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
%         yticks(linspace(0,abs(fixedParams.zw(end)),7))
%         yticklabels(linspace(fixedParams.zw(end),0,7))
%         title(['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'])
%         subtitle('zooplankton abundance')
%         colormap plasma

        
        
end
        
function plt = plot_timeSeries_barplot(var, event, traj, out, ...
    fixedParams, forcing, dat)

nt = fixedParams.nt;
nz = fixedParams.nz;
nPP_size = fixedParams.nPP_size;

% Extract outputs
N = squeeze(out.N(:,:,:,traj));
P = squeeze(out.P(:,:,:,:,traj));
Z = squeeze(out.Z(:,:,:,traj));
OM = squeeze(out.OM(:,:,:,:,traj));

switch var
    
    case 'phytoZooPlankton'
        
        plt = figure;
        plt.Units = 'inches';
        plt.Position = [0 0 8 6];
        % abundance summed over depth
        Y = P(:,:,fixedParams.PP_C_index,:,:);
        Y = repmat(reshape(fixedParams.zwidth, [1 nz]), [nPP_size 1]) .* Y; % conc -> quantity
        Y = squeeze(sum(Y, 2));
        Y = cat(1, Y, sum(fixedParams.zwidth .* Z));
        mY = median(Y, ndims(Y));
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





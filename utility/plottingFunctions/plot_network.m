function plt = plot_network(Type, nutrient, auxVars, FixedParams, Forc, traj, varargin)

extractVarargin(varargin)

if ~exist('axesTextSize', 'var')
    axesTextSize = 12;
end
if ~exist('titleTextSize', 'var')
    titleTextSize = 12;
end
if ~exist('legendTextSize', 'var')
    legendTextSize = 12;
end
if ~exist('legendPointSize', 'var')
    legendPointSize = 36;
end




% plt = figure;

switch Type
    case 'feedingFluxes'
    
%         plt.Units = 'inches';
%         plt.Position = [0 0 16 8];

        % Define network nodes and links
        nodes = table((1:FixedParams.nPP_size + FixedParams.nZP_size)');
        nodes.Properties.VariableNames = {'Label'};
        nodes.ESD = [FixedParams.PPdia; FixedParams.ZPdia];
        nodes.Volume = [FixedParams.PPsize; FixedParams.ZPsize];
        
        trajLabel = [];
        if ischar(traj)
            trajLabel = traj;
            traj = strcmp(Forc.waterMass, traj);
        end
        
        % primary production
        uptake = squeeze(auxVars.PP_uptake_tot(:,:,:,:,traj));
        uptake = mean(uptake, 3, 'omitnan'); % average over trajectories
        nodes.uptake_C = [uptake(:,FixedParams.PP_C_index); ...
            zeros(FixedParams.nZP_size, 1)]; % mmol C / m2 / year(simulation time)
        nodes.uptake_N = [uptake(:,FixedParams.PP_N_index); ...
            zeros(FixedParams.nZP_size, 1)]; % mmol N / m2 / year(simulation time)
        
        % grazing gains
        predation = squeeze(auxVars.predation_gains_tot(:,:,:,:,:,traj));
        predation = mean(predation, 3, 'omitnan'); % average over trajectories
        nodes.grazeGains_C = [zeros(FixedParams.nPP_size, 1); ...
            predation(:,FixedParams.ZP_C_index)];
        nodes.grazeGains_N = [zeros(FixedParams.nPP_size, 1); ...
            predation(:,FixedParams.ZP_N_index)];
                
        % Network links defined by feeding fluxes (predation losses)
        nutrients = {'C','N'};
        for i = 1:length(nutrients)
            Links.(['links' nutrients{i}]) = struct('start', [], 'end', [], 'grazeLosses', [], 'weight' ,[]);
        end
        
        pred_losses = squeeze(auxVars.predation_losses_tot(:,:,:,:,:,traj));
        pred_losses = mean(pred_losses, 4, 'omitnan'); % average over trajectories
        
        % For each prey class, set to zero any feeding flux that accounts for less
        % than x% of the total losses.
%         x = 1;
        x = 5;
        nearZero = pred_losses ./ sum(pred_losses) < (x / 100);
        pred_losses(nearZero) = 0;
        pred_losses = [zeros(FixedParams.nPP_size, FixedParams.nPP_size + FixedParams.nZP_size, size(pred_losses, 3)); ...
            pred_losses]; % extend to square matrix
        for j = 1:length(nutrients)
            for prey = 1:FixedParams.nPP_size + FixedParams.nZP_size
                for pred = 1:FixedParams.nPP_size + FixedParams.nZP_size
                    if pred_losses(pred, prey, j) > 0
                        Links.(['links' nutrients{j}]).start = [Links.(['links' nutrients{j}]).start; prey];
                        Links.(['links' nutrients{j}]).end = [Links.(['links' nutrients{j}]).end; pred];
                        Links.(['links' nutrients{j}]).grazeLosses = [Links.(['links' nutrients{j}]).grazeLosses; pred_losses(pred, prey, j)];
                    end
                end
            end
            Links.(['links' nutrients{j}]).weight = ...
                sqrt(Links.(['links' nutrients{j}]).grazeLosses ./ ...
                max(Links.(['links' nutrients{j}]).grazeLosses));
            Links.(['links' nutrients{j}]) = struct2table(Links.(['links' nutrients{j}]));
        end
        
        % Choose nutrient to plot
        switch nutrient
            case 'carbon'
                pnut = 'C';
            case 'nitrogen'
                pnut = 'N';
        end
        
        links = Links.(['links' pnut]);
        
        % Make network object
        gr = digraph(links.start, links.end);
        
        % Position the nodes
        x = log10(nodes.ESD);
        ang = linspace(0, 45 * pi / 180, FixedParams.nZP_size - 1)';
        y = [zeros(FixedParams.nPP_size, 1); ...
            [0.5; 0.5 + cumsum(tan(ang) .* diff(x(FixedParams.zooplankton)))]];
        
        plt = plot(gr, 'XData', x, 'YData', y, 'lineWidth', 8 * links.weight);
        
        plt.Parent.YTick = [];
        plt.Parent.YTickLabel = [];
        
        set(gca, 'FontSize', axesTextSize)
        
        xnat = [linspace(0.1, 1, 10), linspace(2, 10, 9), linspace(20, 100, 9), linspace(200, 1000, 9)];
        xlog = log10(xnat);
        plt.Parent.XTick = xlog;
        plt.Parent.XTickLabel = xnat;
        plt.Parent.XTickLabel(~ismember(xnat, [1, 10, 100]),:) = ' ';
        xlabel('ESD (\mum)')
        
        % Size nodes relative to production
        uptake = nodes.(['uptake_' pnut]);
        for j = 1:FixedParams.nPP_size + FixedParams.nZP_size
            if FixedParams.phytoplankton(j)
                plt.NodeLabel{j} = num2str(round(uptake(j)));
            end
        end
        uptake = ceil(uptake);
        plt.MarkerSize = [16 * sqrt(uptake(FixedParams.phytoplankton) ./ max(uptake)); ...
            4 * ones(FixedParams.nZP_size, 1)];

        graze = nodes.(['grazeGains_' pnut])';
        for j = 1:FixedParams.nPP_size + FixedParams.nZP_size
            if FixedParams.zooplankton(j)
                plt.NodeLabel{j} = num2str(round(graze(j)));
            end
        end
        graze = ceil(graze);
        plt.MarkerSize = [plt.MarkerSize(FixedParams.phytoplankton), ...
            16 * sqrt(graze(FixedParams.zooplankton) ./ max(graze))] ;
        
        plt.NodeFontSize = 14;
        
        % Add some white-space in front of node labels
        for j = 1:length(plt.NodeLabel)
            plt.NodeLabel{j} = [' ' plt.NodeLabel{j}];
        end

        % Different colours for autotrophs and heterotrophs
        green = [0.1 1 0.2];
        red = [1 0 0];
        cols = [repmat(green, [FixedParams.nPP_size, 1]);
            repmat(red, [FixedParams.nZP_size, 1])];
        plt.NodeColor = cols;
        plt.NodeLabelColor = cols;
        
        plt.EdgeColor = 0.5 .* [1 1 1];
        % Include edge labels
%         labeledge(plt, links.start, links.end, round(links.grazeLosses))
%         plt.EdgeFontSize = 12;
        
        xl = plt.Parent.XLim;
        yl = plt.Parent.YLim;
        
        text(xl(1) + 0.05 * diff(xl), yl(2) - 0.05 * diff(yl), 'autotrophs', 'FontSize', legendTextSize)
        text(xl(1) + 0.05 * diff(xl), yl(2) - 0.15 * diff(yl), 'heterotrophs', 'FontSize', legendTextSize)
        
        hold on
        scatter(xl(1) + 0.025 * diff(xl), yl(2) - 0.05 * diff(yl), ...
            'MarkerFaceColor', green, 'MarkerEdgeColor', green, 'SizeData', legendPointSize)
        scatter(xl(1) + 0.025 * diff(xl), yl(2) - 0.15 * diff(yl), ...
            'MarkerFaceColor', red, 'MarkerEdgeColor', red, 'SizeData', legendPointSize)
        hold off
        
        plt.ArrowSize = 16;

        switch pnut
            case 'C'
                titleText = 'Carbon production and feeding fluxes (mmol C / m^2/ yr)';
            case 'N'
                titleText = 'Nitrogen uptake and feeding fluxes (mmol N / m^2/ yr)';
        end
        if isempty(trajLabel)
            title(titleText, 'FontSize', titleTextSize)
        else
            title({titleText, [trajLabel ' waters']}, 'FontSize', titleTextSize)
        end
    
        
        
    case 'OMfluxes'
        
        % total fluxes to DOM & POM, from autotrophs & heterotrophs summed
        % over size classes
        
%         plt.Units = 'inches';
%         plt.Position = [0 0 8 8];
        
        % Define network nodes and links
        nodes = table({'P'; 'Z'; 'DOM'; 'POM'});
        nodes.Properties.VariableNames = {'Label'};
        
        trajLabel = [];
        if ischar(traj)
            trajLabel = traj;
            traj = strcmp(Forc.waterMass, traj);
        end
%         ntraj = sum(traj > 0);
        
        % OM production comprises mortality (from P and Z) & messy feeding (from Z)
        OM_mort_DOM_tot = squeeze(auxVars.OM_mort_DOM_tot(:,:,:,:,traj));
        OM_mort_POM_tot = squeeze(auxVars.OM_mort_POM_tot(:,:,:,:,traj));
        OM_mort_DOM_tot = mean(OM_mort_DOM_tot, 3, 'omitnan'); % average over trajectories
        OM_mort_POM_tot = mean(OM_mort_POM_tot, 3, 'omitnan'); % average over trajectories
        
        % Define mortality links and sum over size classes
        PtoDOM_mort = sum(OM_mort_DOM_tot(FixedParams.phytoplankton,:));
        PtoPOM_mort = sum(OM_mort_POM_tot(FixedParams.phytoplankton,:));
        ZtoDOM_mort = sum(OM_mort_DOM_tot(FixedParams.zooplankton,:));
        ZtoPOM_mort = sum(OM_mort_POM_tot(FixedParams.zooplankton,:));
        
        OM_mess_DOM_tot = squeeze(auxVars.OM_mess_DOM_tot(:,:,:,:,:,traj));
        OM_mess_POM_tot = squeeze(auxVars.OM_mess_POM_tot(:,:,:,:,:,traj));
        OM_mess_DOM_tot = mean(OM_mess_DOM_tot, 3, 'omitnan'); % average over trajectories
        OM_mess_POM_tot = mean(OM_mess_POM_tot, 3, 'omitnan'); % average over trajectories
        
        % Define messy feeding links and sum over size classes
        ZtoDOM_mess = sum(OM_mess_DOM_tot);
        ZtoPOM_mess = sum(OM_mess_POM_tot);

        % Network links
        nutrients = {'C','N'};
        for i = 1:length(nutrients)
            Links.(['links' nutrients{i}]) = struct('start', [], 'end', [], 'OMflux', [], 'weight' ,[]);
        end
        for j = 1:length(nutrients)
            jnut = find(FixedParams.(['PP_' nutrients{j}, '_index']));
            s = [find(strcmp(nodes.Label,'P')); ...
                find(strcmp(nodes.Label,'P')); ...
                find(strcmp(nodes.Label,'Z')); ...
                find(strcmp(nodes.Label,'Z')); ...
                find(strcmp(nodes.Label,'Z')); ...
                find(strcmp(nodes.Label,'Z'))];
            e = [find(strcmp(nodes.Label,'DOM')); ...
                find(strcmp(nodes.Label,'POM')); ...
                find(strcmp(nodes.Label,'DOM')); ...
                find(strcmp(nodes.Label,'DOM')); ...
                find(strcmp(nodes.Label,'POM')); ...
                find(strcmp(nodes.Label,'POM'))];
            Links.(['links' nutrients{j}]).start = s;
            Links.(['links' nutrients{j}]).end = e;
            Links.(['links' nutrients{j}]).OMflux = [PtoDOM_mort(jnut); PtoPOM_mort(jnut); ...
                ZtoDOM_mort(jnut); ZtoDOM_mess(jnut); ZtoPOM_mort(jnut); ZtoPOM_mess(jnut)];
            Links.(['links' nutrients{j}]).weight = ...
                sqrt(Links.(['links' nutrients{j}]).OMflux ./ ...
                max(Links.(['links' nutrients{j}]).OMflux));
            Links.(['links' nutrients{j}]) = struct2table(Links.(['links' nutrients{j}]));
        end
        
        % Define nodes as total OM production (sum the fluxes)
        Pnode = PtoDOM_mort + PtoPOM_mort;
        Znode = ZtoDOM_mort + ZtoPOM_mort + ZtoDOM_mess + ZtoPOM_mess;
        DOMnode = PtoDOM_mort + ZtoDOM_mort + ZtoDOM_mess;
        POMnode = PtoPOM_mort + ZtoPOM_mort + ZtoPOM_mess;
        nodes.production_C = [Pnode(FixedParams.PP_C_index); ...
            Znode(FixedParams.PP_C_index); DOMnode(FixedParams.PP_C_index); ...
            POMnode(FixedParams.PP_C_index)];
        nodes.production_N = [Pnode(FixedParams.PP_N_index); ...
            Znode(FixedParams.PP_N_index); DOMnode(FixedParams.PP_N_index); ...
            POMnode(FixedParams.PP_N_index)];
        
        
        % Choose nutrient to plot
        switch nutrient
            case 'carbon'
                pnut = 'C';
            case 'nitrogen'
                pnut = 'N';
        end
        
        links = Links.(['links' pnut]);
        
        % Make network object
        gr = digraph(links.start, links.end);
        
        % Position the nodes
        x = [0; 1; 0; 1];
        y = [1; 1; 0; 0];
        % Plot network
        plt = plot(gr, 'XData', x, 'YData', y, 'lineWidth', 8 * links.weight);

        plt.Parent.XTick = [];
        plt.Parent.XTickLabel = [];
        plt.Parent.YTick = [];
        plt.Parent.YTickLabel = [];
                
        % Size nodes relative to production
        production = nodes.(['production_' pnut]);
        for j = 1:length(production)
            plt.NodeLabel{j} = num2str(round(production(j)));
        end
        production = ceil(production);
        plt.MarkerSize = 16 * sqrt(production ./ max(production));
        
        plt.NodeColor = [0 0 0];
        % Different node shapes for plankton and OM
        label = nodes.Label;
        
        label{strcmp(label, 'DOM')} = ['DO' pnut];
        label{strcmp(label, 'POM')} = ['PO' pnut];

        plt.Marker = {'o','o','square','square'};
        text(x(1), y(1) + 0.1, label(1), 'HorizontalAlignment', 'center')
        text(x(2), y(2) + 0.1, label(2), 'HorizontalAlignment', 'center')
        text(x(3), y(3) - 0.1, label(3), 'HorizontalAlignment', 'center')
        text(x(4), y(4) - 0.1, label(4), 'HorizontalAlignment', 'center')
        
        plt.NodeFontSize = 10;
        % Add some white-space in front of node labels
        for j = 1:length(plt.NodeLabel)
            plt.NodeLabel{j} = ['  ' plt.NodeLabel{j}];
        end
        
        % Different colours for mortality and messy feeding
        green = [0 1 0];
        red = [1 0 0];
        plt.EdgeColor = [red; red; ...
            red; green; ...
            red; green];
        
        edgelabel = cell(1, height(links));
        for j = 1:height(links)
            edgelabel{j} = num2str(round(links.OMflux(j)));
        end
        plt.EdgeLabel = edgelabel;
        plt.EdgeFontSize = 12;

        plt.Parent.YLim(2) = 0.2 * diff(plt.Parent.YLim) + plt.Parent.YLim(2);
        
        xl = plt.Parent.XLim;
        yl = plt.Parent.YLim;
        
        text(xl(1) + 0.15 * diff(xl), yl(2) - 0.05 * diff(yl), 'mortality')
        text(xl(1) + 0.15 * diff(xl), yl(2) - 0.1 * diff(yl), 'messy feeding')
        
        hold on
        
        plot(xl(1) + [0.05, 0.1], repmat(yl(2) - 0.05 * diff(yl), [1 2]), 'Color', red);
        plot(xl(1) + [0.05, 0.1], repmat(yl(2) - 0.1 * diff(yl), [1 2]), 'Color', green);
        
        hold off
        
        plt.ArrowSize = 24;
        
        switch pnut
            case 'C'
%                 titleText = 'Organic carbon production and fluxes (mmol C / m^2/ yr)';
                titleText = 'Organic carbon production and fluxes';
            case 'N'
%                 titleText = 'Organic nitrogen production and fluxes (mmol N / m^2/ yr)';
                titleText = 'Organic nitrogen production and fluxes';
        end
        if isempty(trajLabel)
            title(titleText)
        else
            title({titleText, [trajLabel ' waters']})
        end

       
        
        
end


% % Define network nodes and links
% nodes = table((1:FixedParams.nPP_size + FixedParams.nZP_size)');
% nodes.Properties.VariableNames = {'Label'};
% nodes.ESD = [FixedParams.PPdia; FixedParams.ZPdia];
% nodes.Volume = [FixedParams.PPsize; FixedParams.ZPsize];
% 
% traj = 1:10; % choose trajectories (let's group by water origin later on)
% 
% % primary production
% uptake = squeeze(auxVars.PP_uptake_tot(:,:,:,:,traj));
% uptake = mean(uptake, 3); % average over trajectories
% nodes.uptake_C = [uptake(:,FixedParams.PP_C_index); ...
%     zeros(FixedParams.nZP_size, 1)]; % mmol C / m2 / year(simulation time)
% nodes.uptake_N = [uptake(:,FixedParams.PP_N_index); ...
%     zeros(FixedParams.nZP_size, 1)]; % mmol N / m2 / year(simulation time)
% 
% % grazing gains
% predation = squeeze(auxVars.predation_gains_tot(:,:,:,:,:,traj));
% predation = mean(predation, 3); % average over trajectories
% nodes.grazeGains_C = [zeros(FixedParams.nPP_size, 1); ...
%     predation(:,FixedParams.ZP_C_index)];
% nodes.grazeGains_N = [zeros(FixedParams.nPP_size, 1); ...
%     predation(:,FixedParams.ZP_N_index)];
% 
% disp(nodes)
% 
% % Network links defined by feeding fluxes (predation losses)
% nutrients = {'C','N'};
% for i = 1:length(nutrients)
%     Links.(['links' nutrients{i}]) = struct('start', [], 'end', [], 'grazeLosses', [], 'weight' ,[]);
% end
% 
% pred_losses = squeeze(auxVars.predation_losses_tot(:,:,:,:,:,traj));
% pred_losses = mean(pred_losses, 4); % average over trajectories
% 
% % For each prey class, set to zero any feeding flux that accounts for less
% % than x% of the total losses.
% x = 1;
% nearZero = pred_losses ./ sum(pred_losses) < (x / 100);
% pred_losses(nearZero) = 0;
% pred_losses = [zeros(FixedParams.nPP_size, FixedParams.nPP_size + FixedParams.nZP_size, length(nutrients)); ...
%     pred_losses]; % extend to square matrix
% for j = 1:length(nutrients)
%     for prey = 1:FixedParams.nPP_size + FixedParams.nZP_size
%         for pred = 1:FixedParams.nPP_size + FixedParams.nZP_size
%             if pred_losses(pred, prey, j) > 0
%                 Links.(['links' nutrients{j}]).start = [Links.(['links' nutrients{j}]).start; prey];
%                 Links.(['links' nutrients{j}]).end = [Links.(['links' nutrients{j}]).end; pred];
%                 Links.(['links' nutrients{j}]).grazeLosses = [Links.(['links' nutrients{j}]).grazeLosses; pred_losses(pred, prey, j)];
%             end
%         end
%     end
%     Links.(['links' nutrients{j}]).weight = ... 
%         sqrt(Links.(['links' nutrients{j}]).grazeLosses ./ ... 
%         max(Links.(['links' nutrients{j}]).grazeLosses));
%     Links.(['links' nutrients{j}]) = struct2table(Links.(['links' nutrients{j}]));
% end
% 
% % Choose nutrient to plot
% pnut = 'C';
% 
% links = Links.(['links' pnut]);
% 
% gr = digraph(links.start, links.end);
% 
% figure
% 
% plot(gr)
% 
% 
% % Position the nodes
% x = log10(nodes.ESD);
% % y = [repmat(0, [FixedParams.nPP_size, 1]); ... 
% %     repmat(1, [FixedParams.nZP_size, 1])];
% % y = [repmat(0, [FixedParams.nPP_size, 1]); ... 
% %     linspace(0.5, 1.5, FixedParams.nZP_size)'];
% 
% ang = linspace(0, 45 * pi / 180, FixedParams.nZP_size - 1)';
% 
% y = [repmat(0, [FixedParams.nPP_size, 1]); ... 
%     [0.5; 0.5 + cumsum(tan(ang) .* diff(x(FixedParams.zooplankton)))]];
% 
% % plt = plot(gr, 'XData', x, 'YData', y);
% plt = plot(gr, 'XData', x, 'YData', y, 'lineWidth', 8 * links.weight);
% 
% plt.Parent.YTick = [];
% plt.Parent.YTickLabel = [];
% 
% xnat = [linspace(0.1, 1, 10), linspace(2, 10, 9), linspace(20, 100, 9), linspace(200, 1000, 9)];
% xlog = log10(xnat);
% plt.Parent.XTick = xlog;
% plt.Parent.XTickLabel = xnat;
% plt.Parent.XTickLabel(~ismember(xnat, [1, 10, 100]),:) = ' ';
% xlabel('ESD')
% 
% % Size nodes relative to production
% uptake = round(nodes.(['uptake_' pnut]));
% plt.MarkerSize = [16 * sqrt(uptake(FixedParams.phytoplankton) ./ max(uptake)); ... 
%     4 * ones(FixedParams.nZP_size, 1)];
% for j = 1:FixedParams.nPP_size + FixedParams.nZP_size
%     if FixedParams.phytoplankton(j)
%         plt.NodeLabel{j} = num2str(uptake(j));
%     end
% end
% 
% graze = ceil(nodes.(['grazeGains_' pnut]))';
% plt.MarkerSize = [plt.MarkerSize(FixedParams.phytoplankton), ... 
%     16 * sqrt(graze(FixedParams.zooplankton) ./ max(graze))] ;
% 
% for j = 1:FixedParams.nPP_size + FixedParams.nZP_size
%     if FixedParams.zooplankton(j)
%         plt.NodeLabel{j} = num2str(graze(j));
%     end
% end
% 
% plt.NodeFontSize = 16;
% 
% % Different colours for autotrophs and heterotrophs
% green = [0 1 0];
% red = [1 0 0];
% cols = [repmat(green, [FixedParams.nPP_size, 1]);
%     repmat(red, [FixedParams.nZP_size, 1])];
% plt.NodeColor = cols;
% plt.NodeLabelColor = cols;
% 
% plt.EdgeColor = [0 0 0];
% labeledge(plt, links.start, links.end, round(links.grazeLosses))
% plt.EdgeFontSize = 12;
% 
% xl = plt.Parent.XLim;
% yl = plt.Parent.YLim;
% 
% text(xl(1) + 0.05 * diff(xl), yl(2) - 0.05 * diff(yl), 'autotrophs')
% text(xl(1) + 0.05 * diff(xl), yl(2) - 0.1 * diff(yl), 'heterotrophs')
% 
% hold on
% scatter(xl(1) + 0.025 * diff(xl), yl(2) - 0.05 * diff(yl), ...
%     'MarkerFaceColor', green, 'MarkerEdgeColor', green)
% scatter(xl(1) + 0.025 * diff(xl), yl(2) - 0.1 * diff(yl), ...
%     'MarkerFaceColor', red, 'MarkerEdgeColor', red)
% hold off
% 
% switch pnut
%     case 'C'
%         title('Carbon production and feeding fluxes (mmol C / m^2/ yr)')
%     case 'N'
%         title('Nitrogen uptake and feeding fluxes (mmol N / m^2/ yr)')
% end
% 

% Additional plots of model output
% --------------------------------


% run plots.m first to get model output... 

% % add directory with new plotting functions to path
% % later they should be stored in utility/plottingFunctions/...
% addpath('~/Documents/microARC model/newPlottingFunctions')
% 
% folder = '/Documents/microARC model/testnewfuncts'; % save plots here
% 
% 
% % save plots?
% save = true;

%% simple maps of trajectories
map_trajectories(Forc, FixedParams)
map_trajectories(Forc, FixedParams, 'highlightTrajs', [1897 1886 271 855 883])

%% modelled plankton histrograms and size spectra

plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'histogram',...
    'Cbiomass', 'phyto', 200, 'averaged', 'maxDepth', 20, 'event', 10)
plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'histogram',...
    'Cbiomass', 'phyto', 200, 'integrated')
plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'histogram',...
    'Chl', 'phyto', 200, 'integrated', 'waterOrigin', 'Atlantic')

plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'histogram',...
    'abundance', 'phyto', 200, 'integrated', 'waterOrigin', 'Atlantic')

plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'histogram',...
    'biovolume', 'zoo', 200, 'integrated', 'waterOrigin', 'Atlantic')

plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'spectrum',...
    'Cbiomass', 'phyto', 200, 'averaged', 'maxDepth', 20, 'event', 10)

plotModelSizeSpectra(out, auxVars, FixedParams, Data, Forc, 'spectrum',...
    'abundance', 'phyto', 200, 'averaged', 'maxDepth', 20, 'event', 10, 'normalised', true)

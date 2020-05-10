function [v0, initialState] = initialiseVariables_MonodModelSinglePredator(FixedParams,Forc)

% This state variable initialisation depends upon model parameter values,
% but the parameters are uncertain and also require some initialisation...
% Let's just use the parameter values given in Ward 2016 to avoid the state
% variable initialisation varying with our choices of initial parameter
% values.

% Q_C = 1.45e-11 .* [FixedParams.PPsize] .^ 0.8;             % cell carbon content (mmol C / cell)
% Q_N = 16/106 .* Q_C;                                     % cell nitrogen content (mmol N / cell)

% Qmin_N = 0.05;                                           % minimum nitrogen quotas (mmol N / mmol C)
% Qmax_N = 0.17;                                           % maximum nitrogen quotas (mmol N / mmol C)
% Q_N_init =  Qmin_N + 0.75 .* (Qmax_N - Qmin_N);          % initialise with 3/4 full quotas
% theta_N_max = 3;                                         % maximum Chl:N ratio (mg Chl / mmol N)
% Chl_init = 1;                                            % initialise Chl at maximum relative to nitrogen

for iy = 1:length(FixedParams.years)
    y_index = ['y' num2str(FixedParams.years(iy))];
    
    % Inorganic nutrients - dimension = [1, depth, trajectory]
    initialState.(y_index).NO3  = reshape(Forc.(y_index).NO3ic(1,:,:), ...
        [1 FixedParams.nz FixedParams.(y_index).nTraj]);                  % mmol N / m^3 (NO3 output from SINMOD)
    
    % Phytoplankton - dimension = [size, depth, trajectory]
    % Split the SINMOD (PS & PL) plankton estimates evenly across length classes
    sinmod_N_small = squeeze(Forc.(y_index).PSic(1,:,:)); % total N (mmol / m3) at depth and location
    sinmod_N_large = squeeze(Forc.(y_index).PLic(1,:,:));    
    nSmall = sum(~FixedParams.diatoms);
    nLarge = sum(FixedParams.diatoms);
    sinmod_N_small = (1/nSmall) .* repmat(reshape(sinmod_N_small, [1 size(sinmod_N_small)]), [nSmall 1 1]);
    sinmod_N_large = (1/nLarge) .* repmat(reshape(sinmod_N_large, [1 size(sinmod_N_large)]), [nLarge 1 1]);    
    sinmod_N = cat(1,sinmod_N_small,sinmod_N_large);                        % merge arrays of small and large plankton

    initialState.(y_index).PP = sinmod_N;
       
    % Zooplankton - dimension = [1, depth, trajectory]
    % Zooplankton are not modelled at any particularly size classes, so
    % just initialise the their total nitrogen content as some fraction of
    % the total phytoplanktonnitrogen content
    ZP_PP_frac = 0.25;
    initialState.(y_index).ZP = ZP_PP_frac .* sum(initialState.(y_index).PP);

    % Organic matter - dimension = [1, depth, trajectory]
    initialState.(y_index).OM = 0.05 * (sum(initialState.(y_index).ZP) ...
        + sum(initialState.(y_index).PP));

    % Store state variables in array:
    % 1st dimension = variable
    % 2nd dimension = location (trajectory)
    % Order of variables = inorganic nutrients [depth] ->
    %                      phytoplankton [size, depth] ->
    %                      zooplankton [depth] -> 
    %                      organic matter [depth]
    
    
    v0.(y_index) = [squeeze(initialState.(y_index).NO3); ... 
        reshape(initialState.(y_index).PP, ... 
        [FixedParams.nPP*FixedParams.nz FixedParams.(y_index).nTraj]); ...
        squeeze(initialState.(y_index).ZP); ... 
        squeeze(initialState.(y_index).OM)];
        
    
    
%     v0.(y_index) = nan(FixedParams.nz * FixedParams.nVar, ... 
%         FixedParams.(y_index).nTraj);
% 
%     j = true;
%     i = 1;
%     k = 1;
%     % Inorganic nutrients    
%     while j
%         x = initialState.(y_index).(FixedParams.INtype{k});
%         y = numel(x) / FixedParams.(y_index).nTraj;
%         v0.(y_index)(i:y+i-1,:) = reshape(x, [y FixedParams.(y_index).nTraj]);
%         i = i+y;
%         k = k+1;
%         if k > FixedParams.nIN, j = false; end
%     end
%     
%     % Phytoplankton
%     x = initialState.(y_index).PP;
%     y = numel(x) / FixedParams.(y_index).nTraj;
%     v0.(y_index)(i:y+i-1,:) = reshape(x, [y FixedParams.(y_index).nTraj]);
%     i = i+y;
%     
%     % Zooplankton
%     x = initialState.(y_index).ZP;
%     y = numel(x) / FixedParams.(y_index).nTraj;
%     v0.(y_index)(i:y+i-1,:) = reshape(x, [y FixedParams.(y_index).nTraj]);    
%     i = i+y;
%     
%     % Organic matter
%     x = initialState.(y_index).OM;
%     y = numel(x) / FixedParams.(y_index).nTraj;
%     v0.(y_index)(i:y+i-1,:) = reshape(x, [y FixedParams.(y_index).nTraj]);    
    
end

end


function [v0, initialState] = initialiseVariables(FixedParams,Forc)

for iy = 1:length(FixedParams.years)
    y_index = ['y' num2str(FixedParams.years(iy))];
    
    % Inorganic nutrients - dimension = [1, depth, trajectory]
    initialState.(y_index).NO3  = reshape(Forc.(y_index).NO3ic(1,:,:), ...
        [1 FixedParams.nz FixedParams.(y_index).nTraj]);                    % mmol N / m^3 (NO3 output from SINMOD)
    
    % Phytoplankton - dimension = [size, depth, trajectory]
    % Split the SINMOD (PS & PL) plankton estimates evenly across length classes
    sinmod_N_small = squeeze(Forc.(y_index).PSic(1,:,:));                   % total N (mmol / m3) at depth and location
    sinmod_N_large = squeeze(Forc.(y_index).PLic(1,:,:));
    nSmall = sum(~FixedParams.diatoms);
    nLarge = sum(FixedParams.diatoms);
    sinmod_N_small = (1/nSmall) .* repmat(reshape(sinmod_N_small, ...
        [1 size(sinmod_N_small)]), [nSmall 1 1]);
    sinmod_N_large = (1/nLarge) .* repmat(reshape(sinmod_N_large, ...
        [1 size(sinmod_N_large)]), [nLarge 1 1]);
    sinmod_N = cat(1,sinmod_N_small,sinmod_N_large);                        % merge arrays of small and large plankton
    
    initialState.(y_index).PP = sinmod_N;
    
    % Zooplankton - dimension = [1, depth, trajectory]
    % Zooplankton are not assigned any particular size class, so just
    % initialise their total nitrogen content as some fraction of the total
    % phytoplankton nitrogen content
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
    
end

end


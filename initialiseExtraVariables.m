function [namesExtra, nExtra, AUXVARS, AUXVARS_2d, AUXVARS_3d] = ...
    initialiseExtraVariables(v0, parameterList, Forc, returnExtra)

nt = parameterList.FixedParams.nt;
nz = parameterList.FixedParams.nz;
% nPP = parameterList.FixedParams.nPP;
nPP_size = parameterList.FixedParams.nPP_size;
nZP_size = parameterList.FixedParams.nZP_size;
nPP_nut = parameterList.FixedParams.nPP_nut;
nTraj = Forc.nTraj;

% Process 1st trajectory separately to find dimension of extra output
forcing.T = Forc.T(:,:,1);
forcing.K = Forc.K(:,:,1);
forcing.PARsurf = Forc.PARsurf(:,:,1);

[~, extraOutput_1d, extraOutput_2d, extraOutput_3d] = ...
    ODEs(0, v0(:,1), parameterList, forcing, 2, returnExtra);

namesExtra_1d = fieldnames(extraOutput_1d);
namesExtra_2d = fieldnames(extraOutput_2d);
namesExtra_3d = fieldnames(extraOutput_3d);

nExtra = [length(namesExtra_1d) length(namesExtra_2d) length(namesExtra_3d)];
namesExtra = cat(1, namesExtra_1d, namesExtra_2d, namesExtra_3d);

AUXVARS = nan(nExtra(1) * nz, nt, nTraj);
AUXVARS(:,1,1) = struct2array(extraOutput_1d);

AUXVARS_2d = nan(nExtra(2) * (nPP_size + nZP_size) * nz, nt, nTraj);
AUXVARS_2d(:,1,1) = struct2array(structfun(@(x)x(:)', ...
    extraOutput_2d, 'UniformOutput', false));

AUXVARS_3d = nan(nExtra(3) * (nPP_size + nZP_size) * nz * nPP_nut, nt, nTraj);
AUXVARS_3d(:,1,1) = struct2array(structfun(@(x)x(:)', ...
    extraOutput_3d, 'UniformOutput', false));

for i = 2:nTraj  % Loop through remaining trajectories
    forcing.T = Forc.T(:,:,i);
    forcing.K = Forc.K(:,:,i);
    forcing.PARsurf = Forc.PARsurf(:,:,i);
    [~, extraOutput_1d, extraOutput_2d, extraOutput_3d] = ...
        ODEs(0, v0(:,i), parameterList, forcing, 2, returnExtra);
    AUXVARS(:,1,i) = struct2array(extraOutput_1d);    
    AUXVARS_2d(:,1,i) = struct2array(structfun(@(x)x(:)', ...
        extraOutput_2d, 'UniformOutput', false));    
    AUXVARS_3d(:,1,i) = struct2array(structfun(@(x)x(:)', ...
        extraOutput_3d, 'UniformOutput', false));
end

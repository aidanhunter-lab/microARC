function [namesExtra, nExtra, AUXVARS, AUXVARS_2d] = ... 
    initialiseExtraVariables(v0, parameterList, Forc)

nt = parameterList.FixedParams.nt;
nz = parameterList.FixedParams.nz;
% nPP = parameterList.FixedParams.nPP;
nPP_size = parameterList.FixedParams.nPP_size;
nTraj = Forc.nTraj;

% Process 1st trajectory separately to find dimension of extra output
forcing.T = Forc.T(:,:,1);
forcing.K = Forc.K(:,:,1);
forcing.PARsurf = Forc.PARsurf(:,:,1);
[~, extraOutput, extraOutput_2d] = ODEs(0, v0(:,1), parameterList, forcing, 2, true);
namesExtra = fieldnames(extraOutput);
namesExtra_2d = fieldnames(extraOutput_2d);
nExtra = [length(namesExtra) length(namesExtra_2d)];
namesExtra = cat(1, namesExtra, namesExtra_2d);
AUXVARS = nan(nExtra(1) * nz, nt, nTraj);
AUXVARS(:,1,1) = struct2array(extraOutput);
AUXVARS_2d = nan(nExtra(2) * nPP_size * nz, nt, nTraj);
AUXVARS_2d(:,1,1) = struct2array(structfun(@(x)x(:)', ...
    extraOutput_2d, 'UniformOutput', false));
for i = 2:nTraj  % Loop through remaining trajectories
    forcing.T = Forc.T(:,:,i);
    forcing.K = Forc.K(:,:,i);
    forcing.PARsurf = Forc.PARsurf(:,:,i);
    [~, extraOutput, extraOutput_2d] = ODEs(0, v0(:,i), parameterList, forcing, 2, true);
    AUXVARS(:,1,i) = struct2array(extraOutput);
    AUXVARS_2d(:,1,i) = struct2array(structfun(@(x)x(:)', ...
        extraOutput_2d, 'UniformOutput', false));
end

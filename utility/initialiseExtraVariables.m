function [namesExtra, dimsExtra, indexExtra, AUXVARS] = ...
    initialiseExtraVariables(v0, parameterList, Forc, returnExtra)

nt = parameterList.FixedParams.nt;
nTraj = Forc.nTraj;

% Process 1st trajectory separately to find dimension of extra output
forcing.T = Forc.T(:,:,1);
forcing.K = Forc.K(:,:,1);
forcing.PARsurf = Forc.PARsurf(:,:,1);
forcing.lat = Forc.y(:,1)';
forcing.yd = yearday(Forc.t(:,1)');

[~, extraOutput] = ODEs(0, v0(:,1), parameterList, forcing, 2, returnExtra);

namesExtra = fieldnames(extraOutput);
nExtra = length(namesExtra);
dimsExtra = structfun(@(x) size(x), extraOutput, 'UniformOutput', false);

extraOutput = structfun(@(x) x(:)', extraOutput, 'UniformOutput', false);

x = struct2array(extraOutput);
n = length(x);

lab = cell(n,1);
jj = 0;
for i = 1:length(namesExtra)
    ni = prod(dimsExtra.(namesExtra{i}));
    lab(jj+1:jj+ni) = namesExtra(i);
    jj = jj + ni;
end

if nExtra > 0
    for i = 1:length(namesExtra)
        indexExtra.(namesExtra{i}) = strcmp(lab, namesExtra{i});
    end
else
    indexExtra = struct();
end


AUXVARS = nan(n, nt, nTraj);
AUXVARS(:,1,1) = x; clear x


for i = 2:nTraj  % Loop through remaining trajectories
    forcing.T = Forc.T(:,:,i);
    forcing.K = Forc.K(:,:,i);
    forcing.PARsurf = Forc.PARsurf(:,:,i);
    forcing.lat = Forc.y(:,i)';
    forcing.yd = yearday(Forc.t(:,i)');
    [~, extraOutput] = ODEs(0, v0(:,i), parameterList, forcing, 2, returnExtra);
    extraOutput = structfun(@(x) x(:)', extraOutput, 'UniformOutput', false);
    AUXVARS(:,1,i) = struct2array(structfun(@(x) x(:)', extraOutput, ...
        'UniformOutput', false));    
end

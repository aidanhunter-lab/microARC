function [namesExtra, nExtra, AUXVARS, AUXVARS_2d, RATES] = ... 
    initialiseExtraVariables(v0, parameterList, Forc)

returnExtras = parameterList.FixedParams.returnExtras;

nt = parameterList.FixedParams.nt;
nz = parameterList.FixedParams.nz;
nPP = parameterList.FixedParams.nPP;
nEquations = parameterList.FixedParams.nEquations;
nTraj = Forc.nTraj;

if strcmp(returnExtras, 'none')
    namesExtra = []; nExtra = []; % AUXVARS = []; RATES = [];
    AUXVARS = nan(1, nt, nTraj);
    AUXVARS_2d = nan(nPP, nt, nTraj);
    RATES = nan(1, nt, nTraj);            
else    
    % Process 1st trajectory separately to find dimension of extra output
    forcing.T = Forc.T(:,:,1);
    forcing.K = Forc.K(:,:,1);
    forcing.PARsurf = Forc.PARsurf(:,:,1);
    switch returnExtras
        case 'auxiliary'
            RATES = nan(1, nt, nTraj);
            [~, extraOutput, extraOutput_2d] = ODEs(0, v0(:,1), parameterList, forcing, 2);
            namesExtra = fieldnames(extraOutput);
            nExtra = zeros(1,2);
            nExtra(1) = length(namesExtra);
            AUXVARS = nan(nExtra(1) * nz, nt, nTraj);
            AUXVARS(:,1,1) = struct2array(extraOutput);            
            namesExtra_2d = fieldnames(extraOutput_2d);
            nExtra(2) = length(namesExtra_2d);
            namesExtra = cat(1, namesExtra, namesExtra_2d);            
            AUXVARS_2d = nan(nExtra(2) * nPP * nz, nt, nTraj);
            AUXVARS_2d(:,1,1) = struct2array(structfun(@(x)x(:), ... 
                extraOutput_2d, 'UniformOutput', false));
            for i = 2:nTraj
                forcing.T = Forc.T(:,:,i);
                forcing.K = Forc.K(:,:,i);
                forcing.PARsurf = Forc.PARsurf(:,:,i);
                [~, extraOutput] = ODEs(0, v0(:,i), parameterList, forcing, 2);
                AUXVARS(:,1,i) = struct2array(extraOutput);
                AUXVARS_2d(:,1,i) = struct2array(structfun(@(x)x(:), ... 
                    extraOutput_2d, 'UniformOutput', false));
            end
        case 'rates'
            AUXVARS = nan(1, nt, nTraj); namesExtra = []; nExtra = [];
            AUXVARS_2d = nan(nPP, nt, nTraj); namesExtra = []; nExtra = [];
            RATES = nan(nEquations, nt, nTraj);            
            RATES(:,1,1) = ODEs(0, v0(:,1), parameterList, forcing, 2);                        
            for i = 2:nTraj
                forcing.T = Forc.T(:,:,i);
                forcing.K = Forc.K(:,:,i);
                forcing.PARsurf = Forc.PARsurf(:,:,i);
                RATES(:,1,i) = ODEs(0, v0(:,i), parameterList, forcing, 2);
            end
        case 'auxiliaryAndRates'
            [dvdt, extraOutput, extraOutput_2d] = ODEs(0, v0(:,1), parameterList, forcing, 2);                        
            RATES = nan(nEquations, nt, nTraj);            
            RATES(:,1,1) = dvdt;
            namesExtra = fieldnames(extraOutput);
            nExtra = zeros(1,2);
            nExtra(1) = length(namesExtra);
            AUXVARS = nan(nExtra(1) * nz, nt, nTraj);
            AUXVARS(:,1,1) = struct2array(extraOutput);
            namesExtra_2d = fieldnames(extraOutput_2d);
            nExtra(2) = length(namesExtra_2d);
            namesExtra = cat(1, namesExtra, namesExtra_2d);
            AUXVARS_2d = nan(nExtra(2) * nPP * nz, nt, nTraj);
            AUXVARS_2d(:,1,1) = struct2array(structfun(@(x)x(:), ...
                extraOutput_2d, 'UniformOutput', false));
            for i = 2:nTraj
                forcing.T = Forc.T(:,:,i);
                forcing.K = Forc.K(:,:,i);
                forcing.PARsurf = Forc.PARsurf(:,:,i);
                [dvdt, extraOutput, extraOutput_2d] = ODEs(0, v0(:,i), parameterList, forcing, 2);
                RATES(:,1,i) = dvdt;
                AUXVARS(:,1,i) = struct2array(extraOutput);
                AUXVARS_2d(:,1,i) = struct2array(structfun(@(x)x(:), ...
                    extraOutput_2d, 'UniformOutput', false));
            end
    end    
end


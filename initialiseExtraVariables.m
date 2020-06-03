function [namesExtra, nExtra, AUXVARS, RATES] = ... 
    initialiseExtraVariables(v0, parameterList, Forc)

returnExtras = parameterList.FixedParams.returnExtras;

nt = parameterList.FixedParams.nt;
nz = parameterList.FixedParams.nz;
nEquations = parameterList.FixedParams.nEquations;
nTraj = Forc.nTraj;

if strcmp(returnExtras, 'none')
    namesExtra = []; nExtra = []; % AUXVARS = []; RATES = [];
    AUXVARS = nan(1, nt, nTraj);
    RATES = nan(1, nt, nTraj);            
else    
    % Process 1st trajectory separately to find dimension of extra output
    forcing.T = Forc.T(:,:,1);
    forcing.K = Forc.K(:,:,1);
    forcing.PARsurf = Forc.PARsurf(:,:,1);
    switch returnExtras
        case 'auxiliary'
            RATES = nan(1, nt, nTraj);            
            [~, extraOutput] = ODEs(0, v0(:,1), parameterList, forcing, 2);                        
            namesExtra = fieldnames(extraOutput);
            nExtra = length(namesExtra);            
            AUXVARS = nan(nExtra * nz, nt, nTraj);            
            AUXVARS(:,1,1) = struct2array(extraOutput);
            for i = 2:nTraj
                forcing.T = Forc.T(:,:,i);
                forcing.K = Forc.K(:,:,i);
                forcing.PARsurf = Forc.PARsurf(:,:,i);
                [~, extraOutput] = ODEs(0, v0(:,i), parameterList, forcing, 2);
                AUXVARS(:,1,i) = struct2array(extraOutput);
            end
        case 'rates'
            AUXVARS = nan(1, nt, nTraj); namesExtra = []; nExtra = [];
            RATES = nan(nEquations, nt, nTraj);            
            RATES(:,1,1) = ODEs(0, v0(:,1), parameterList, forcing, 2);                        
            for i = 2:nTraj
                forcing.T = Forc.T(:,:,i);
                forcing.K = Forc.K(:,:,i);
                forcing.PARsurf = Forc.PARsurf(:,:,i);
                RATES(:,1,i) = ODEs(0, v0(:,i), parameterList, forcing, 2);
            end
        case 'auxiliaryAndRates'
            [dvdt, extraOutput] = ODEs(0, v0(:,1), parameterList, forcing, 2);                        
            RATES = nan(nEquations, nt, nTraj);            
            RATES(:,1,1) = dvdt;
            namesExtra = fieldnames(extraOutput);
            nExtra = length(namesExtra);            
            AUXVARS = nan(nExtra * nz, nt, nTraj);            
            AUXVARS(:,1,1) = struct2array(extraOutput);
            for i = 2:nTraj
                forcing.T = Forc.T(:,:,i);
                forcing.K = Forc.K(:,:,i);
                forcing.PARsurf = Forc.PARsurf(:,:,i);
                [dvdt, extraOutput] = ODEs(0, v0(:,i), parameterList, forcing, 2);
                RATES(:,1,i) = dvdt;
                AUXVARS(:,1,i) = struct2array(extraOutput);
            end
    end    
end


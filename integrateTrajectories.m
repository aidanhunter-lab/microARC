function [OUT, AUXVARS, AUXVARS_2d, RATES, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options)

nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nOM = FixedParams.nOM;
nTraj = Forc.nTraj;
nEquations = FixedParams.nEquations;


N_index = FixedParams.IN_index;
P_index = FixedParams.PP_index;
Z_index = FixedParams.ZP_index;
OM_index = FixedParams.OM_index;

T = Forc.T;
K = Forc.K;
PARsurf = Forc.PARsurf;
deepens = Forc.deepens;
infillDepth = Forc.infillDepth;
replicateConc = Forc.replicateConc;
dry = ~Forc.wet;

returnExtras = FixedParams.returnExtras;

parameterList.FixedParams = FixedParams;
parameterList.Params = Params;


OUT = nan(nEquations, nt, nTraj); % stores state variable outputs
OUT(:,1,:) = reshape(v0, [nEquations 1 nTraj]);

% Create and initialise arrays for extra output if needed
[namesExtra, nExtra, AUXVARS, AUXVARS_2d, RATES] = ... 
    initialiseExtraVariables(v0, parameterList, Forc);


% Loop through trajectories and integrate
parfor i = 1:nTraj
    % Forcing data
    forcing = struct();
    forcing.T = T(:,:,i);
    forcing.K = K(:,:,i);
    forcing.PARsurf = PARsurf(:,:,i);
    % Initial state
    v_in = v0(:,i);
    % Integrate step-wise between successive data points
    
    for j = 2:nt
        if deepens(:,j,i)
            % Extract state variable types from input vector
            N = v_in(N_index);
            P = reshape(v_in(P_index), [nPP nz]);
            Z = v_in(Z_index);
            OM = reshape(v_in(OM_index), [nOM nz]);
            % Infill values
            nfill = sum(infillDepth(:,j,i));
            N(infillDepth(:,j,i)) = repmat(N(replicateConc(:,j,i)), [nfill 1]);
            P(:,infillDepth(:,j,i)) = repmat(P(:,replicateConc(:,j,i)), [1 nfill]);
            Z(infillDepth(:,j,i)) = repmat(Z(replicateConc(:,j,i)), [nfill 1]);
            if nOM > 1
                OM(:,infillDepth(:,j,i)) = repmat(OM(:,replicateConc(:,j,i)), [1 nfill]);
            else
                OM(infillDepth(:,j,i)) = repmat(OM(replicateConc(:,j,i)), [nfill 1]);
            end
            % Recombine the input vector
            v_in = [N; P(:); Z; OM(:)];            
        end
        % Integrate
        sol=ode45(@(t, v_in) ODEs(t, v_in, parameterList, forcing, j), [0 1], v_in, ode45options);
        % Store solutions each day (each forcing data time-step)
        OUT(:,j,i) = deval(sol, 1);
        % Update initials for next time step
        v_in = OUT(:,j,i);
        % Extract extra outputs
        if ~strcmp(returnExtras, 'none')
            switch returnExtras
                case 'auxiliary'
                    [~, extraOutput, extraOutput_2d] = ... 
                        ODEs(1, v_in, parameterList, forcing, j);
                    AUXVARS(:,j,i) = struct2array(extraOutput);
                    AUXVARS_2d(:,j,i) = struct2array(structfun(@(x)x(:)', ...
                        extraOutput_2d, 'UniformOutput', false));
                case 'auxiliaryAndRates'
                    [dvdt, extraOutput, extraOutput_2d] = ODEs(1, v_in, parameterList, forcing, j);
                    RATES(:,j,i) = dvdt;
                    AUXVARS(:,j,i) = struct2array(extraOutput);
                    AUXVARS_2d(:,j,i) = struct2array(struct2fun(@(x)x(:)', ...
                        extraOutput_2d, 'UniformOutput', false));
                case 'rates'
                    RATES(:,j,i) = ODEs(1, v_in, parameterList, forcing, j);
            end
        end            
    end
end


% Omit values deeper than sea floor
if any(dry(:))    
    N = OUT(N_index,:,:);
    P = reshape(OUT(P_index,:,:), [nPP nz nt nTraj]);
    Z = OUT(Z_index,:,:);    
    OM = reshape(OUT(OM_index,:,:), [nOM nz nt nTraj]);
%     OM = OUT(OM_index,:,:);
    N(dry) = nan;
    P(repmat(reshape(dry, [1 nz nt nTraj]), [nPP 1 1 1])) = nan;
    P = reshape(P, [nPP * nz nt nTraj]);
    Z(dry) = nan;
    OM(repmat(reshape(dry, [1 nz nt nTraj]), [nOM 1 1 1])) = nan;
    OM = reshape(OM, [nOM * nz nt nTraj]);
%     OM(dry) = nan;
    OUT = [N; P; Z; OM];
    if ~strcmp(returnExtras, 'none')
        switch returnExtras
            case 'auxiliary'
                AUXVARS = reshape(AUXVARS, [nz nExtra(1) nt nTraj]);
                AUXVARS(repmat(reshape(dry, [nz 1 nt nTraj]), [1 nExtra(1) 1 1])) = nan;
                AUXVARS_2d = reshape(AUXVARS_2d, [nPP nz nExtra(2) nt nTraj]);
                AUXVARS_2d(repmat(reshape(dry, [1 nz 1 nt nTraj]), [nPP 1 nExtra(2) 1 1])) = nan;
            case 'rates'
                N_rate = RATES(N_index,:,:);
                P_rate = reshape(RATES(P_index,:,:), [nPP nz nt nTraj]);
                Z_rate = RATES(Z_index,:,:);
                OM_rate = RATES(OM_index,:,:);
                N_rate(dry) = nan;
                P_rate(repmat(reshape(dry, [1 nz nt nTraj]), [nPP 1 1 1])) = nan;
                P_rate = reshape(P_rate, [nPP * nz nt nTraj]);
                Z_rate(dry) = nan;
                OM_rate(dry) = nan;                
                RATES = [N_rate; P_rate; Z_rate; OM_rate];
            case 'auxiliaryAndRates'
                AUXVARS = reshape(AUXVARS, [nz nExtra nt nTraj]);
                AUXVARS(repmat(reshape(dry, [nz 1 nt nTraj]), [1 nExtra 1 1])) = nan;
                AUXVARS_2d = reshape(AUXVARS_2d, [nPP nz nExtra(2) nt nTraj]);
                AUXVARS_2d(repmat(reshape(dry, [1 nz 1 nt nTraj]), [nPP 1 nExtra(2) 1 1])) = nan;
                N_rate = RATES(N_index,:,:);
                P_rate = reshape(RATES(P_index,:,:), [nPP nz nt nTraj]);
                Z_rate = RATES(Z_index,:,:);
                OM_rate = RATES(OM_index,:,:);
                N_rate(dry) = nan;
                P_rate(repmat(reshape(dry, [1 nz nt nTraj]), [nPP 1 1 1])) = nan;
                P_rate = reshape(P_rate, [nPP * nz nt nTraj]);
                Z_rate(dry) = nan;
                OM_rate(dry) = nan;                
                RATES = [N_rate; P_rate; Z_rate; OM_rate];
        end
    end
end


% Tidy up
switch returnExtras
    case 'none'
        AUXVARS = []; AUXVARS_2d = []; RATES = []; nExtra = 0; namesExtra = [];
    case 'auxiliary'
        AUXVARS = reshape(AUXVARS, [nz nExtra(1) nt nTraj]);
        AUXVARS_2d = reshape(AUXVARS_2d, [nPP nz nExtra(2) nt nTraj]);
        RATES = [];
    case 'rates'
        AUXVARS = []; AUXVARS_2d = []; nExtra = 0; namesExtra = [];
    case 'auxiliaryAndRates'
        AUXVARS = reshape(AUXVARS, [nz nExtra(1) nt nTraj]);
        AUXVARS_2d = reshape(AUXVARS_2d, [nPP nz nExtra(2) nt nTraj]);        
end

function [OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions)

nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nPP_size = FixedParams.nPP_size;
nPP_nut = FixedParams.nPP_nut;
nOM = FixedParams.nOM;
nOM_type = FixedParams.nOM_type;
nOM_nut = FixedParams.nOM_nut;
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

parameterList.FixedParams = FixedParams;
parameterList.Params = Params;

OUT = nan(nEquations, nt, nTraj); % stores state variable outputs
OUT(:,1,:) = reshape(v0, [nEquations 1 nTraj]);

% Create and initialise arrays for extra output if needed
[namesExtra, nExtra, AUXVARS, AUXVARS_2d] = ... 
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
    % Integrating method
    odeSolver = str2func(odeIntegrator);
    % Integrate step-wise between successive data points    
    for j = 2:nt
        if deepens(:,j,i)
            % Extract state variable types from input vector
            N = v_in(N_index);
            P = reshape(v_in(P_index), [nPP_size nz nPP_nut]);
            Z = v_in(Z_index);
            OM = reshape(v_in(OM_index), [nOM_type nz nOM_nut]);
            % Infill values
            nfill = sum(infillDepth(:,j,i));
            N(infillDepth(:,j,i)) = repmat(N(replicateConc(:,j,i)), [nfill 1]);
            P(:,infillDepth(:,j,i),:) = repmat(P(:,replicateConc(:,j,i),:), [1 nfill 1]);
            Z(infillDepth(:,j,i)) = repmat(Z(replicateConc(:,j,i)), [nfill 1]);
            OM(:,infillDepth(:,j,i),:) = repmat(OM(:,replicateConc(:,j,i),:), [1 nfill 1]);
            % Recombine the input vector
            v_in = [N; P(:); Z; OM(:)];            
        end
        % Integrate
%         sol=ode45(@(t, v_in) ODEs(t, v_in, parameterList, forcing, j, false), [0 1], v_in, odeOptions);        
%         sol=ode23(@(t, v_in) ODEs(t, v_in, parameterList, forcing, j, false), [0 1], v_in, odeOptions);        
        sol = odeSolver(@(t, v_in) ODEs(t, v_in, parameterList, forcing, j, false), [0 1], v_in, odeOptions);        
        % Store solutions each day (each forcing data time-step)
        OUT(:,j,i) = deval(sol, 1);
        % Update initials for next time step
        v_in = OUT(:,j,i);
        % Extract extra outputs
        [~, extraOutput, extraOutput_2d] = ...
            ODEs(1, v_in, parameterList, forcing, j, true);
        AUXVARS(:,j,i) = struct2array(extraOutput);
        AUXVARS_2d(:,j,i) = struct2array(structfun(@(x)x(:)', ...
            extraOutput_2d, 'UniformOutput', false));
    end
end


% Omit values deeper than sea floor
if any(dry(:))    
    % state variables
    N = OUT(N_index,:,:);
    P = reshape(OUT(P_index,:,:), [nPP_size nz nPP_nut nt nTraj]);
    Z = OUT(Z_index,:,:);    
    OM = reshape(OUT(OM_index,:,:), [nOM_type nz nOM_nut nt nTraj]);
    N(dry) = nan;
    P(repmat(reshape(dry, [1 nz 1 nt nTraj]), [nPP_size 1 nPP_nut 1 1])) = nan;
    P = reshape(P, [nPP * nz nt nTraj]);
    Z(dry) = nan;
    OM(repmat(reshape(dry, [1 nz 1 nt nTraj]), [nOM_type 1 nOM_nut 1 1])) = nan;
    OM = reshape(OM, [nOM * nz nt nTraj]);
    OUT = [N; P; Z; OM];
    % extra outputs
    AUXVARS = reshape(AUXVARS, [nz nExtra(1) nt nTraj]);
    AUXVARS(repmat(reshape(dry, [nz 1 nt nTraj]), [1 nExtra(1) 1 1])) = nan;
    AUXVARS_2d = reshape(AUXVARS_2d, [nPP_size nz nExtra(2) nt nTraj]);
    AUXVARS_2d(repmat(reshape(dry, [1 nz 1 nt nTraj]), [nPP_size 1 nExtra(2) 1 1])) = nan;
end


% Tidy up
AUXVARS = reshape(AUXVARS, [nz nExtra(1) nt nTraj]);
AUXVARS_2d = reshape(AUXVARS_2d, [nPP_size nz nExtra(2) nt nTraj]);


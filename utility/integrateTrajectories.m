function [out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, ...
    v0, odeIntegrator, odeOptions, varargin)

extractVarargin(varargin)

nt = FixedParams.nt;
nz = FixedParams.nz;
% nPP = FixedParams.nPP;
nPP_size = FixedParams.nPP_size;
nPP_nut = FixedParams.nPP_nut;
% nZP = FixedParams.nZP;
nZP_size = FixedParams.nZP_size;
nZP_nut = FixedParams.nZP_nut;
nsize = nPP_size + nZP_size;
% nOM = FixedParams.nOM;
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

if ~exist('returnExtra', 'var')
% By default return all auxiliary variables (specified at end of ODEs.m)
    returnExtra = true;
end

% Initialise arrays for extra output
[namesExtra, dimsExtra, indexExtra, AUXVARS] = ... 
    initialiseExtraVariables(v0, parameterList, Forc, returnExtra);

if Forc.integrateFullTrajectory
    nt_traj = repmat(nt, [1 nTraj]);
else
   x = cumsum(Forc.sampleTime);
   [nt_traj,~] = find([zeros(1, nTraj); diff(x == max(x))]);
   nt_traj = nt_traj';
end

% clear Forc

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
        if j <= nt_traj(i)
            if deepens(:,j,i)
                % Extract state variable types from input vector
                N = v_in(N_index);
                P = reshape(v_in(P_index), [nPP_size nz nPP_nut]);
                Z = reshape(v_in(Z_index), [nZP_size nz nZP_nut]);
                OM = reshape(v_in(OM_index), [nOM_type nz nOM_nut]);
                % Infill values
                nfill = sum(infillDepth(:,j,i));
                N(infillDepth(:,j,i)) = repmat(N(replicateConc(:,j,i)), [nfill 1]);
                P(:,infillDepth(:,j,i),:) = repmat(P(:,replicateConc(:,j,i),:), [1 nfill 1]);
                Z(:,infillDepth(:,j,i),:) = repmat(Z(:,replicateConc(:,j,i),:), [1 nfill 1]);
                OM(:,infillDepth(:,j,i),:) = repmat(OM(:,replicateConc(:,j,i),:), [1 nfill 1]);
                % Recombine the input vector
                v_in = [N; P(:); Z(:); OM(:)];
            end
            % Integrate
            sol = odeSolver(@(t, v_in) ODEs(t, v_in, parameterList, forcing, j, false), [0 1], v_in, odeOptions);
            % Store solutions each day (each forcing data time-step)
            OUT(:,j,i) = deval(sol, 1);
            % Update initials for next time step
            v_in = OUT(:,j,i);
            % Extract extra outputs
            [~, extraOutput] = ODEs(1, v_in, parameterList, forcing, j, returnExtra);
            AUXVARS(:,j,i) = struct2array(structfun(@(x) x(:)', ... 
                extraOutput, 'UniformOutput', false));
        else
            break;
        end
    end
end


%~~~~~~~~~~
% Tidy up

% Extract solutions from array, OUT, into a more readable struct, out.
out.N = reshape(OUT(N_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(P_index,:,:), [nPP_size nz nPP_nut nt nTraj]);
out.Z = reshape(OUT(Z_index,:,:), [nZP_size nz nZP_nut nt nTraj]);
out.OM = reshape(OUT(OM_index,:,:), [nOM_type nz nOM_nut nt nTraj]);
% Same for extra outputs.
nExtra = length(namesExtra);
if nExtra > 0
    for i = 1:length(namesExtra)
%         auxVars.(namesExtra{i}) = AUXVARS(indexExtra.(namesExtra{i}),:,:);
%         auxVars.(namesExtra{i}) = reshape(auxVars.(namesExtra{i}), ...
%             [dimsExtra.(namesExtra{i}), nt, nTraj]);
        auxVars.(namesExtra{i}) = reshape( ...
            AUXVARS(indexExtra.(namesExtra{i}),:,:), ...
            [dimsExtra.(namesExtra{i}), nt, nTraj]);
    end
else
    auxVars = struct();
end

clear OUT AUXVARS

% Omit values deeper than sea floor
if any(dry(:))
    
    % state variables
    out.N(dry) = nan;
    out.P(repmat(reshape(dry, [1, nz, 1, nt, nTraj]), ... 
        [nPP_size, 1, nPP_nut, 1, 1])) = nan;
    out.Z(repmat(reshape(dry, [1, nz, 1, nt, nTraj]), ... 
        [nZP_size, 1, nZP_nut, 1, 1])) = nan;
    out.OM(repmat(reshape(dry, [1, nz, 1, nt, nTraj]), ...
        [nOM_type, 1, nOM_nut, 1, 1])) = nan;
    
    % auxiliary outputs
    for i = 1:nExtra
        d = dimsExtra.(namesExtra{i});
        nd = length(d);
        % Include an if-statement for each separate matrix dimension
        % contained in auxVars
        if nd == 2 && all(d == [1, nz])
            auxVars.(namesExtra{i})(dry) = nan;
        end
        if nd == 3 && all(d == [nPP_size + nZP_size, nz, nPP_nut])
            auxVars.(namesExtra{i})(repmat(reshape(dry, ... 
                [1, nz, 1, nt, nTraj]), [nPP_size + nZP_size, 1, nPP_nut, 1, 1])) = nan;
        end
        if nd == 3 && all(d == [nZP_size, nPP_size + nZP_size, nz])
            auxVars.(namesExtra{i})(repmat(reshape(dry, ... 
                [1, 1, nz, nt, nTraj]), [nZP_size, nPP_size + nZP_size, 1, 1, 1])) = nan;
        end        
        if nd == 4 && all(all(d == [nZP_size, 1, nz, nPP_nut]))
            auxVars.(namesExtra{i})(repmat(reshape(dry, ... 
                [1, 1, nz, 1, nt, nTraj]), [nZP_size, 1, 1, nPP_nut, 1])) = nan;
        end
        if nd == 2 && all(d == [nPP_size + nZP_size, nz])
            auxVars.(namesExtra{i})(repmat(reshape(dry, ... 
                [1, nz,nt, nTraj]), [nPP_size + nZP_size, 1, 1, 1])) = nan;
        end
    end
end


%~~~~~~~~~~~~~~~~~~~
% Derived quantities
%~~~~~~~~~~~~~~~~~~~

% These could be returned from the ODEs.m script but it's more
% memory-efficient to calculate quantities here...

% This should be moved to a unique function called 'derivedQuantities.m'

if (islogical(returnExtra) && returnExtra) || any(strcmp(returnExtra, 'all'))
    
    dt = unique(diff(Forc.t(:,1)));
    
    % Primary production
    auxVars.PP_uptake = out.P(:,:,FixedParams.PP_C_index,:,:) .* ... 
        auxVars.V(FixedParams.phytoplankton,:,:,:,:); % mmol x / m3 / day
%     auxVars.PP_uptake = out.P .* auxVars.V(FixedParams.phytoplankton,:,:,:,:); % mmol x / m3 / day
    auxVars.PP_uptake_depthInt = sum(reshape(FixedParams.zwidth, [1 nz]) .* auxVars.PP_uptake , 2); % depth-integrated production, mmol x / m2 / day
    auxVars.PP_uptake_tot = sum(dt .* auxVars.PP_uptake_depthInt, ...
        find(size(auxVars.PP_uptake_depthInt) == nt)); % total production of each size class along each trajectory, mmol C / m2 / year (simulation time) 
    
    % Feeding fluxes (grazing loss) and secondary production (grazing gain)
    Z_C = out.Z(:,:,FixedParams.ZP_C_index,:,:); % zooplankton carbon    
    auxVars.predation_losses = reshape(auxVars.Q(:,:,~FixedParams.PP_Chl_index,:,:), ...
        [1, nPP_size + nZP_size, nz, nZP_nut, nt, nTraj]) .* ...
        reshape(Z_C, [nZP_size, 1, nz, 1, nt, nTraj]) .* ...
        reshape(auxVars.G, [nZP_size, nPP_size + nZP_size, nz, 1, nt, nTraj]);
    auxVars.predation_losses_depthInt = sum(reshape(FixedParams.zwidth, [1, 1, nz]) .* ... 
        auxVars.predation_losses, 3);
    auxVars.predation_losses_tot = sum(dt .* auxVars.predation_losses_depthInt, ...
        find(size(auxVars.predation_losses_depthInt) == nt));
    
    predation_gains_all = auxVars.lambda(:,:,:,~FixedParams.PP_Chl_index,:,:) .* auxVars.predation_losses;
    mess = auxVars.predation_losses - predation_gains_all;
    
    auxVars.predation_gains = sum(predation_gains_all, 2);    
    auxVars.predation_gains_depthInt = sum(reshape(FixedParams.zwidth, [1, 1, nz]) .* ...
        auxVars.predation_gains, 3);    
    auxVars.predation_gains_tot = sum(dt .* auxVars.predation_gains_depthInt, ...
        find(size(auxVars.predation_gains_depthInt) == nt));

    % OM production
    % Messy feeding    
    beta_mess = reshape(Params.beta, [1, nsize]) .* mess;
    auxVars.OM_mess_DOM = sum(beta_mess, 2);
    auxVars.OM_mess_POM = sum(mess - beta_mess, 2);
    auxVars.OM_mess_DOM_depthInt = sum(reshape(FixedParams.zwidth, [1, 1, nz]) .* ...
        auxVars.OM_mess_DOM, 3);
    auxVars.OM_mess_DOM_tot = sum(dt .* auxVars.OM_mess_DOM_depthInt, ...
        find(size(auxVars.OM_mess_DOM_depthInt) == nt));
    auxVars.OM_mess_POM_depthInt = sum(reshape(FixedParams.zwidth, [1, 1, nz]) .* ...
        auxVars.OM_mess_POM, 3);
    auxVars.OM_mess_POM_tot = sum(dt .* auxVars.OM_mess_POM_depthInt, ...
        find(size(auxVars.OM_mess_POM_depthInt) == nt));    
    % Mortality
    B = [out.P(:,:,~FixedParams.PP_Chl_index,:,:); out.Z];
    mortality = (Params.m + Params.m2 .* B) .* B;
%     mortality = Params.m .* B;
    auxVars.OM_mort_DOM = Params.beta .* mortality;
    auxVars.OM_mort_POM = mortality - auxVars.OM_mort_DOM;
    auxVars.OM_mort_DOM_depthInt = sum(reshape(FixedParams.zwidth, [1, nz]) .* ...
        auxVars.OM_mort_DOM, 2);
    auxVars.OM_mort_DOM_tot = sum(dt .* auxVars.OM_mort_DOM_depthInt, ...
        find(size(auxVars.OM_mort_DOM_depthInt) == nt));
    auxVars.OM_mort_POM_depthInt = sum(reshape(FixedParams.zwidth, [1, nz]) .* ...
        auxVars.OM_mort_POM, 2);
    auxVars.OM_mort_POM_tot = sum(dt .* auxVars.OM_mort_POM_depthInt, ...
        find(size(auxVars.OM_mort_POM_depthInt) == nt));

    
end

end

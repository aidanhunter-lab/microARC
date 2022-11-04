function [out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, ...
    v0, odeIntegrator, odeOptions, varargin)

extractVarargin(varargin)

nt = FixedParams.nt;
nz = FixedParams.nz;
nPP_size = FixedParams.nPP_size;
nPP_nut = FixedParams.nPP_nut;
nZP_size = FixedParams.nZP_size;
nZP_nut = FixedParams.nZP_nut;
nOM_type = FixedParams.nOM_type;
nOM_nut = FixedParams.nOM_nut;
nEquations = FixedParams.nEquations;

N_index = FixedParams.IN_index;
P_index = FixedParams.PP_index;
Z_index = FixedParams.ZP_index;
OM_index = FixedParams.OM_index;

nTraj = Forc.nTraj;
T = Forc.T;
K = Forc.K;
Y = Forc.y;
TimeYD = yearday(Forc.t);

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

% Loop through trajectories and integrate
parfor i = 1:nTraj        
    % Forcing data
    forcing = struct();
    forcing.T = T(:,:,i);
    forcing.K = K(:,:,i);
    forcing.PARsurf = PARsurf(:,:,i);
    
    % also transfer t and y for calculation of day length and improved light  
    forcing.lat = Y(:,i)';
    forcing.yd = TimeYD(:,i)';

    
    % Initial state
    v_in = v0(:,i);
    % Integrating method
    odeSolver = str2func(odeIntegrator);
    % Integrate step-wise between successive data points
    for j = 2:nt  % daily loop
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
    
    auxVars = derivedQuantities(Forc, FixedParams, Params, out, auxVars);
    
end

end

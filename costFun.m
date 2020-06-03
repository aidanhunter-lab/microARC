function cost = costFun(x, FixedParams, Params, Forc, Data, v0, ode45options)

% Returns a scalar describing model misfit to data given parameter
% updates, x

nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nTraj = Forc.nTraj;

% Set parameter values
Params = updateParameters2(Params, FixedParams, x);

% Integrate model
[OUT, AUXVARS, RATES, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options);

% Extract solutions
out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [1 nz nt nTraj]);

returnExtras = FixedParams.returnExtras;
if ~strcmp(returnExtras, 'none')
    switch returnExtras
        case 'auxiliary'
            for k = 1:nExtra
                auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
            end
        case 'auxiliaryAndRates'
            for k = 1:nExtra
                auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
            end
            rates.N = reshape(RATES(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
            rates.P = reshape(RATES(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
            rates.Z = reshape(RATES(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
            rates.OM = reshape(RATES(FixedParams.OM_index,:,:), [1 nz nt nTraj]);
        case 'rates'
            rates.N = reshape(RATES(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
            rates.P = reshape(RATES(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
            rates.Z = reshape(RATES(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
            rates.OM = reshape(RATES(FixedParams.OM_index,:,:), [1 nz nt nTraj]);
    end
end


% Compare model to observations
for i = 1:Data.nEvent
    ti = Data.eventTraj(i,:);
    N = out.N(:,:,:,ti);
    P = out.P(:,:,:,ti);
    Z = out.Z(:,:,:,ti);
    OM = out.OM(:,:,:,ti);
end










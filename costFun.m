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
[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options);

% Extract solutions
out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [1 nz nt nTraj]);

if sum(nExtra) > 0
    for k = 1:nExtra(1)
        auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
    end
    for k = 1:nExtra(2)
        auxVars.(namesExtra{k+nExtra(1)}) = squeeze(AUXVARS_2d(:,:,k,:,:));
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










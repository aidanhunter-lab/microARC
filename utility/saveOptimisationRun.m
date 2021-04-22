function m = saveOptimisationRun(fileName, fittingOutput, Data, Forc, FixedParams, Params, v0)

m = matfile(fileName, 'Writable', true);
m.fittingOutput = fittingOutput;
m.Data = Data;
m.Forc = Forc;
m.FixedParams = FixedParams;
m.Params = Params;
if isempty(v0)
    v0 = 'Initial values were recalculated for each parameter set because some values are initialised in relation to quota parameters.';
end
m.v0 = v0;
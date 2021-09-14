function [m, fittingOutput, parnames, optPar, lb, ub, Data, Forc, FixedParams, ...
    Params, v0] = loadOptimisationRun(fileName)

m = load(fileName);
fittingOutput = m.fittingOutput;
parnames = fittingOutput.parNames;
optPar = fittingOutput.optPar;
lb = fittingOutput.lowerBound;
ub = fittingOutput.upperBound;
Data = m.Data;
Forc = m.Forc;
FixedParams = m.FixedParams;
Params = m.Params;
v0 = m.v0;

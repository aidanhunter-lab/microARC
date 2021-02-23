function [m, gaOutput, parnames, optPar, lb, ub, Data, Forc, FixedParams, ...
    Params, v0] = loadOptimisationRun(fileName)

m = load(fileName);
% m = matfile(fileName, 'Writable', true);
gaOutput = m.gaOutput;
parnames = gaOutput.parNames;
optPar = gaOutput.optPar;
lb = gaOutput.lowerBound;
ub = gaOutput.upperBound;
Data = m.Data;
Forc = m.Forc;
FixedParams = m.FixedParams;
Params = m.Params;
v0 = m.v0;

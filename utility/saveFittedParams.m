function saveFittedParams(fileName, optPar, parNames)
optPar = optPar(:);
parNames = parNames(:);
optPars_table = table(parNames, optPar);
optPars_table.Properties.VariableNames = {'Param','Value'};
saveObj = matfile(fileName);
saveObj.Properties.Writable = true;
saveObj.pars = optPars_table;
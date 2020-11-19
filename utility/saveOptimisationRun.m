function m = saveOptimisationRun(fileName, gaOutput, Data, Forc, FixedParams, Params, v0)
        
m = matfile(fileName, 'Writable', true);
m.gaOutput = gaOutput;
m.Data = Data;
m.Forc = Forc;
m.FixedParams = FixedParams;
m.Params = Params;
m.v0 = v0;

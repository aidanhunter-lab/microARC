function FixedParams = createIndexes(fixedParams)

% Indexes used to extract specific state variables from a single vector

FixedParams = fixedParams;

nz  = FixedParams.nz;
nIN = FixedParams.nIN;
nPP = FixedParams.nPP;
nZP = FixedParams.nZP;
nOM = FixedParams.nOM;

% Inorganic nutrients
FixedParams.IN_index = [true(1, nIN * nz) ...
                        false(1, (nPP + nZP + nOM) * nz)]';
                    
% Phytoplankton
FixedParams.PP_index = [false(1, nIN * nz) ...
                        true(1, nPP * nz) ... 
                        false(1, (nZP + nOM) * nz)]';
                    
% Zooplankton
FixedParams.ZP_index = [false(1, (nIN + nPP) * nz) ...
                        true(1, nZP * nz) ...
                        false(1, nOM * nz)]';

% Organic matter
FixedParams.OM_index = [false(1, (nIN + nPP + nZP) * nz) ...
                        true(1, nOM * nz)]';

end
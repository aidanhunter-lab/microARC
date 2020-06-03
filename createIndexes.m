function [FixedParams, IN_index, PP_index, ZP_index, OM_index] = createIndexes(fixedParams)

% Indexes used to extract specific state variables from a single vector

FixedParams = fixedParams;

nz  = FixedParams.nz;
nIN = FixedParams.nIN;
nPP = FixedParams.nPP;
nZP = FixedParams.nZP;
nOM = FixedParams.nOM;

% Inorganic nutrients
IN_index = [true(1, nIN * nz) false(1, (nPP + nZP + nOM) * nz)]';
FixedParams.IN_index = IN_index;
                    
% Phytoplankton
PP_index = [false(1, nIN * nz)  true(1, nPP * nz) false(1, (nZP + nOM) * nz)]';
FixedParams.PP_index = PP_index;
                    
% Zooplankton
ZP_index = [false(1, (nIN + nPP) * nz) true(1, nZP * nz) false(1, nOM * nz)]';
FixedParams.ZP_index = ZP_index;

% Organic matter
OM_index = [false(1, (nIN + nPP + nZP) * nz) true(1, nOM * nz)]';
FixedParams.OM_index = OM_index;

end
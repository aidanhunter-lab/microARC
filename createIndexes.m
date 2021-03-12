function [FixedParams, IN_index, PP_index, ZP_index, OM_index] = createIndexes(fixedParams)

% Indexes used to extract specific state variables from a single vector

FixedParams = fixedParams;

nz  = FixedParams.nz;
nIN = FixedParams.nIN;
nPP = FixedParams.nPP;
nPP_nut = FixedParams.nPP_nut;
nZP = FixedParams.nZP;
nZP_nut = FixedParams.nZP_nut;
nOM = FixedParams.nOM;
nOM_type = FixedParams.nOM_type;
nOM_nut = FixedParams.nOM_nut;

% Inorganic nutrients
IN_index = [true(1, nIN * nz) false(1, (nPP + nZP + nOM) * nz)]';
FixedParams.IN_index = IN_index;
                    
% Phytoplankton
PP_index = [false(1, nIN * nz)  true(1, nPP * nz) false(1, (nZP + nOM) * nz)]';
FixedParams.PP_index = PP_index;
if nPP_nut > 1
    % index phytoplankton nutrients
    for i = 1:FixedParams.nPP_nut
        FixedParams.(['PP_' FixedParams.PP_nut{i} '_index']) = ...
            strcmp(FixedParams.PP_nut{i}, FixedParams.PP_nut);
    end
end

% Zooplankton
ZP_index = [false(1, (nIN + nPP) * nz) true(1, nZP * nz) false(1, nOM * nz)]';
FixedParams.ZP_index = ZP_index;
if nZP_nut > 1
    % index heterotroph nutrients
    for i = 1:FixedParams.nZP_nut
        FixedParams.(['ZP_' FixedParams.ZP_nut{i} '_index']) = ...
            strcmp(FixedParams.ZP_nut{i}, FixedParams.ZP_nut);
    end
end


% Organic matter
OM_index = [false(1, (nIN + nPP + nZP) * nz) true(1, nOM * nz)]';
FixedParams.OM_index = OM_index;

if nOM_type > 1
    % index DOM and POM
    for i =1:FixedParams.nOM_type
        FixedParams.([FixedParams.OM_type{i} '_index']) = ... 
            strcmp(FixedParams.OM_type{i}, FixedParams.OM_type);
    end
end
if nOM_nut > 1
    % index OM nutrients
    for i =1:FixedParams.nOM_nut
        FixedParams.(['OM_' FixedParams.OM_nut{i} '_index']) = ... 
            strcmp(FixedParams.OM_nut{i}, FixedParams.OM_nut);
    end
end







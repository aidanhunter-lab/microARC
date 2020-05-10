function FixedParams = createIndexes_MonodModelSinglePredator(fixedParams)

FixedParams = fixedParams;

nz  = FixedParams.nz;
nIN = FixedParams.nIN;
nPP = FixedParams.nPP;
nZP = FixedParams.nZP;
nOM = FixedParams.nOM;

%~~~~~~~~~~~~~~~~~~~~
% Inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~

FixedParams.IN_index = [true(1, nIN * nz) ...
                        false(1, (nPP + nZP + nOM) * nz)]';                 % all inorganic nutrients

% FixedParams.iDIC = find(strcmp(FixedParams.INtype, 'DIC'));                 % specific nutrients
% FixedParams.iNH4 = find(strcmp(FixedParams.INtype, 'NH4'));
% FixedParams.iNO2 = find(strcmp(FixedParams.INtype, 'NO2'));
% FixedParams.iNO3 = find(strcmp(FixedParams.INtype, 'NO3'));


%~~~~~~~~~~~~~~
% Phytoplankton
%~~~~~~~~~~~~~~

FixedParams.PP_index = [false(1, nIN * nz) ...
                        true(1, nPP * nz) ... 
                        false(1, (nZP + nOM) * nz)]';                       % all phytoplankton

% FixedParams.ipC = find(strcmp(FixedParams.PPnut, 'C'));                     % specific nutrients in phytoplankton
% FixedParams.ipN = find(strcmp(FixedParams.PPnut, 'N'));
% FixedParams.ipChl = find(strcmp(FixedParams.PPnut, 'Chl'));
% FixedParams.ipCN = [FixedParams.ipC FixedParams.ipN];


%~~~~~~~~~~~~
% Zooplankton
%~~~~~~~~~~~~

FixedParams.ZP_index = [false(1, (nIN + nPP) * nz) ...
                        true(1, nZP * nz) ...
                        false(1, nOM * nz)]';                               % all zooplankton

% FixedParams.izC = find(strcmp(FixedParams.ZPnut, 'C'));                     % specific nutrients in zooplankton
% FixedParams.izN = find(strcmp(FixedParams.ZPnut, 'N'));
% FixedParams.izCN = [FixedParams.izC FixedParams.izN];


%~~~~~~~~~~~~~~~
% Organic matter
%~~~~~~~~~~~~~~~

FixedParams.OM_index = [false(1, (nIN + nPP + nZP) * nz) ...
                        true(1, nOM * nz)]';                                % all organic matter

% FixedParams.ioC = find(strcmp(FixedParams.OMnut, 'C'));                     % nutrients in organic matter
% FixedParams.ioN = find(strcmp(FixedParams.OMnut, 'N'));

% FixedParams.kDOM = find(strcmp(FixedParams.OMtype, 'DOM'));                 % type (DOM or POM) of organic matter
% FixedParams.kPOM = find(strcmp(FixedParams.OMtype, 'POM'));

end
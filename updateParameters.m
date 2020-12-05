function Params = updateParameters(params,FixedParams,varargin)

Params = params;

if ~isempty(varargin)
    if nargin == 3 && isnumeric(varargin{1}) % if new param values are passed as a single vector (with labels stored in FixedParams)
        newvals = varargin{1};
        parnames = FixedParams.tunePars;
        npars = length(parnames);
        V_PP = FixedParams.PPsize;
        for i = 1:npars
            pn = parnames{i};
            Params.(pn) = newvals(i); % update scalars
            % update size-dependent parameters - constructed from scalars
            % suffixed with _a and _b
            if ~isempty(regexp(pn,'_a', 'once')) || ~isempty(regexp(pn,'_b', 'once'))
                Name = pn(1:end-2);
                param_a = Params.([Name '_a']);
                param_b = Params.([Name '_b']);
                if ~(any(strcmp(Name, {'Qmax_delQ'})) || any(strcmp(Name, {'wp'})))
                    Params.(Name) = powerFunction(param_a, param_b, V_PP);
                else
                    if any(strcmp(Name, {'wp'}))
                        Params.(Name) = powerFunction(param_a, param_b, V_PP, ...
                            'UpperBound', FixedParams.maxSinkSpeed, 'Transpose', true);
                    end
                    if any(strcmp(Name, {'Qmax_delQ'}))
                        Params.(Name) = 1 ./ (1 - powerFunction(param_a, param_b, V_PP));
                    end
                end
            end
            
            if strcmp(pn, 'beta1') || strcmp(pn, 'beta2') || strcmp(pn, 'beta3')
                Params.beta = doubleLogisticFunction(Params.beta1, Params.beta2, ...
                    Params.beta3, log10(V_PP));
                Params.beta(FixedParams.nPP_size+1) = Params.beta(FixedParams.nPP_size); % assume beta for zooplankton is equivalent to largest phytoplankton size class
            end
            
        end
        
        if any(strcmp('Qmin_QC_a', parnames)) || any(strcmp('Qmin_QC_b', parnames)) || ...
                any(strcmp('Qmax_delQ_a', parnames)) || any(strcmp('Qmax_delQ_b', parnames))
            Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
            Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;
        end
        
        if any(strcmp('Vmax_QC_a', parnames)) || any(strcmp('Vmax_QC_b', parnames)) || ...
                any(strcmp('aN_QC_a', parnames)) || any(strcmp('aN_QC_b', parnames))
            Params.kN = Params.Vmax_QC ./ Params.aN_QC;
        end
        
        if any(strcmp('rDOC', parnames)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_C_index) = Params.rDOC; end
        if any(strcmp('rDON', parnames)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_N_index) = Params.rDON; end
        if any(strcmp('rPOC', parnames)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_C_index) = Params.rPOC; end
        if any(strcmp('rPON', parnames)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_N_index) = Params.rPON; end
        
        if any(strcmp('wDOM', parnames)), Params.wk(FixedParams.DOM_index) = Params.wDOM; end
        if any(strcmp('wPOM', parnames)), Params.wk(FixedParams.POM_index) = Params.wPOM; end
        
    else % if new params are passed as name-values pairs
        
        % Scalar parameters
        parNames = Params.scalars;
        for i = 1:length(parNames)
            position = find(strcmpi(varargin, parNames{i}));
            if ~isempty(position)
                name = parNames{i};
                value = varargin{position+1};
                Params.(name) = value;
            end
        end
        
        % Size-dependent parameters
        V_PP = FixedParams.PPsize(:);
        parNames = Params.sizeDependent;
        suf = {'_a', '_b'};
        for i = 1:length(parNames)
            position = find(strcmpi(varargin, parNames{i}));
            if ~isempty(position)
                name = varargin{position};
                value = varargin{position+1};
                Params.(name) = value; % update the scalar
                % Most size-dependent terms are parameterised by power functions of
                % size using two parameters subscripted with '_a' or '_b', one
                % of which has just been updated... now recalculate the vector
                if any(strcmp(name(end-1:end), suf))
                    Name = name(1:length(name)-2);
                    param_a = Params.([Name '_a']);
                    param_b = Params.([Name '_b']);
                    if ~(any(strcmp(Name, {'Qmax_delQ'})) || any(strcmp(Name, {'wp'})))
                        Params.(Name) = powerFunction(param_a, param_b, V_PP);
                    else
                        if any(strcmp(Name, {'wp'}))
                            Params.(Name) = powerFunction(param_a, param_b, V_PP, ...
                                'UpperBound', FixedParams.maxSinkSpeed, 'Transpose', true);
                        end
                        if any(strcmp(Name, {'Qmax_delQ'}))
                            Params.(Name) = 1 ./ (1 - powerFunction(param_a, param_b, V_PP));
                        end
                    end
                    
                    if strcmp(name, 'beta1') || strcmp(name, 'beta2') || strcmp(name, 'beta3')
                        Params.beta = doubleLogisticFunction(Params.beta1, Params.beta2, ...
                            Params.beta3, log10(V_PP)); % flexible 3-parameter double logistic function
                        Params.beta(FixedParams.nPP_size+1) = Params.beta(FixedParams.nPP_size); % assume beta for zooplankton is equivalent to largest phytoplankton size class
                    end
                end
            end
        end
        
        
        if any(strcmp('Qmin_QC_a', varargin)) || any(strcmp('Qmin_QC_b', varargin)) || ...
                any(strcmp('Qmax_delQ_a', varargin)) || any(strcmp('Qmax_delQ_b', varargin))
            Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
            Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;
        end
        
        if any(strcmp('Vmax_QC_a', varargin)) || any(strcmp('Vmax_QC_b', varargin)) || ...
                any(strcmp('aN_QC_a', varargin)) || any(strcmp('aN_QC_b', varargin))
            Params.kN = Params.Vmax_QC ./ Params.aN_QC;
        end
        
        if any(strcmp('rDOC', varargin)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_C_index) = Params.rDOC; end
        if any(strcmp('rDON', varargin)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_N_index) = Params.rDON; end
        if any(strcmp('rPOC', varargin)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_C_index) = Params.rPOC; end
        if any(strcmp('rPON', varargin)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_N_index) = Params.rPON; end
        
        if any(strcmp('wDOM', varargin)), Params.wk(FixedParams.DOM_index) = Params.wDOM; end
        if any(strcmp('wPOM', varargin)), Params.wk(FixedParams.POM_index) = Params.wPOM; end
        
    end
end

end

function y = powerFunction(a, b, x, varargin)
y = a .* x .^ b;
if ~isempty(varargin)
    tr = strcmp(varargin, 'Transpose');
    if any(tr) && varargin{find(tr)+1}, y = y'; end
    ub = strcmp(varargin, 'UpperBound');
    if any(ub)
        cap = varargin{find(ub)+1};
        y(y>cap) = cap;
    end
end
end

function y = doubleLogisticFunction(a, b, c, x)
u = exp(x - c);
y = a ./ (1 + u) .* (1 + b .* u);
end

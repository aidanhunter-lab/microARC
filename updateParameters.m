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
                if ~any(strcmp(Name, {'Qmax_delQ'}))
                    Params.(Name) = volumeDependent(param_a, param_b, V_PP);
                else
                    Params.(Name) = 1 ./ (1 - volumeDependent(param_a, param_b, V_PP));
                end
            end
        end
        
        if any(strcmp('Qmin_QC_a', parnames)) || any(strcmp('Qmin_QC_b', parnames)) || ...
                any(strcmp('Qmax_delQ_a', parnames)) || any(strcmp('Qmax_delQ_b', parnames))
            Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
            Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;
        end
        
        if any(strcmp('rDOC', parnames)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_C_index) = Params.rDOC; end
        if any(strcmp('rDON', parnames)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_N_index) = Params.rDON; end
        if any(strcmp('rPOC', parnames)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_C_index) = Params.rPOC; end
        if any(strcmp('rPON', parnames)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_N_index) = Params.rPON; end

        if any(strcmp('wk', parnames))
            % wk probably shouldn't be numerically optimised because it's
            % related to the maximum permissible integration time step
            Params.wk = [0 Params.wk];
        end
        
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
                    if ~any(strcmp(Name, {'Qmax_delQ'}))
                        Params.(Name) = volumeDependent(param_a, param_b, V_PP);
                    else
                        Params.(Name) = 1 ./ (1 - volumeDependent(param_a, param_b, V_PP));
                    end
                end
            end
        end
        
        
        if any(strcmp('Qmin_QC_a', varargin)) || any(strcmp('Qmin_QC_b', varargin)) || ...
                any(strcmp('Qmax_delQ_a', varargin)) || any(strcmp('Qmax_delQ_b', varargin))
            Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
            Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;
        end
        
        if any(strcmp('rDOC', varargin)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_C_index) = Params.rDOC; end        
        if any(strcmp('rDON', varargin)), Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_N_index) = Params.rDON; end        
        if any(strcmp('rPOC', varargin)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_C_index) = Params.rPOC; end        
        if any(strcmp('rPON', varargin)), Params.rOM(FixedParams.POM_index,1,FixedParams.OM_N_index) = Params.rPON; end        
        
        if any(strcmp('wk', varargin))
            % wk probably shouldn't be numerically optimised because it's
            % related to the maximum permissible integration time step
            Params.wk = [0 Params.wk];
        end
        
    end
end


function p = volumeDependent(a,b,V)
p = a .* V .^ b;

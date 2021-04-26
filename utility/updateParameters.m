function Params = updateParameters(Params,FixedParams,varargin)

if ~isempty(varargin)
    if nargin == 3 && isnumeric(varargin{1}) % if new param values are passed as a single vector (with labels stored in FixedParams)
        newvals = varargin{1};
        parnames = FixedParams.tunePars;
        npars = length(parnames);
        Vol = FixedParams.sizeAll;
%         V_PP = FixedParams.PPsize;
        for i = 1:npars
            pn = parnames{i};
            Params.(pn) = newvals(i); % update scalars
            % update size-dependent parameters - constructed from scalars
            % suffixed with _a and _b
            if ~isempty(regexp(pn,'_a', 'once')) || ~isempty(regexp(pn,'_b', 'once'))
                Name = pn(1:end-2);
                param_a = Params.([Name '_a']);
                param_b = Params.([Name '_b']);
                if ~any(strcmp(Name, 'm'))
                    Params.(Name) = Params.([Name '_func'])(param_a, param_b, Vol);
                else
                    Params.(Name) = Params.([Name '_func'])(param_a, param_b, FixedParams.m_min, Vol);
                end
                
%                 if ~(any(strcmp(Name, {'Qmax_delQ'})) || any(strcmp(Name, {'wp'})))
%                     Params.(Name) = powerFunction(param_a, param_b, Vol);
%                 else
%                     if any(strcmp(Name, {'wp'}))
%                         Params.(Name) = powerFunction(param_a, param_b, Vol, ...
%                             'UpperBound', FixedParams.maxSinkSpeed_P, 'Transpose', true);
%                     end
%                     if any(strcmp(Name, {'Qmax_delQ'}))
%                         Params.(Name) = 1 ./ (1 - powerFunction(param_a, param_b, Vol));
%                     end
%                 end
            end
            
            if strcmp(pn, 'beta1') || strcmp(pn, 'beta2') || strcmp(pn, 'beta3')
                Name = pn(1:end-1);
                param1 = Params.([Name '1']);
                param2 = Params.([Name '2']);
                param3 = Params.([Name '3']);
                Params.(Name) = Params.([Name '_func'])(param1, param2, param3, Vol);
                Params.beta = doubleLogisticFunction(Params.beta1, Params.beta2, ...
                    Params.beta3, log10(Vol));
            end
            
        end
        
        % Functions of parameters
        
        if any(strcmp('Qmin_QC_a', parnames)) || any(strcmp('Qmin_QC_b', parnames)) || ...
                any(strcmp('Qmax_delQ_a', parnames)) || any(strcmp('Qmax_delQ_b', parnames))
            Params.Qmax_QC = Params.Qmax_QC_func(Params.Qmin_QC, Params.Qmax_delQ);
            Params.delQ_QC = Params.delQ_QC_func(Params.Qmin_QC, Params.Qmax_QC);
%             Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
%             Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;
        end
        
        if any(strcmp('Vmax_QC_a', parnames)) || any(strcmp('Vmax_QC_b', parnames)) || ...
                any(strcmp('aN_QC_a', parnames)) || any(strcmp('aN_QC_b', parnames))
            Params.kN = Params.kN_func(Params.Vmax_QC, Params.aN_QC);
%             Params.kN = Params.Vmax_QC ./ Params.aN_QC;
        end
        
        if any(strcmp('Gmax_a', parnames)) || any(strcmp('Gmax_b', parnames))
            Params.Gmax = Params.Gmax(FixedParams.zooplankton);
        end
        
        if any(strcmp('rDOC', parnames)) || any(strcmp('rDON', parnames)) || ...
                any(strcmp('rPOC', parnames)) || any(strcmp('rPON', parnames))
            Params.rOM = Params.rOM_func(Params.rDOC, Params.rDON, ...
                Params.rPOC, Params.rPON, FixedParams.DOM_index, ...
                FixedParams.POM_index, FixedParams.OM_C_index, FixedParams.OM_N_index);
        end
        
        if any(strcmp('wDOM', parnames)) || any(strcmp('wPOM', parnames))
            Params.wk = Params.wk_func(Params.wDOM, Params.wPOM, ... 
                FixedParams.DOM_index, FixedParams.POM_index);
        end
        
                
    else % if new params are passed as name-values pairs
        
        extractVarargin(varargin)
        
        % Scalar parameters
        parNames = Params.scalars;
        for i = 1:length(parNames)
            if exist(parNames{i}, 'var')
                name = parNames{i};
                value = eval(parNames{i});
                Params.(name) = value;
            end
        end
        
        % Size-dependent parameters
        Vol = FixedParams.sizeAll(:);
%         V_PP = FixedParams.PPsize(:);
        parNames = Params.sizeDependent;
        suf = {'_a', '_b'};
        
        for i = 1:length(parNames)
            
            if exist(parNames{i}, 'var')
            
                name = parNames{i};
                value = eval(parNames{i});
                Params.(name) = value; % update the scalar

                % Most size-dependent terms are parameterised by power functions of
                % size using two parameters subscripted with '_a' or '_b', one
                % of which has just been updated... now recalculate the vector
                
                if any(strcmp(name(end-1:end), suf))
                    Name = name(1:length(name)-2);
                    param_a = Params.([Name '_a']);
                    param_b = Params.([Name '_b']);
                    if ~strcmp(Name, 'm')
                        Params.(Name) = Params.([Name '_func'])(param_a, param_b, Vol);
                    else
                        Params.(Name) = Params.([Name '_func'])(param_a, param_b, FixedParams.m_min, Vol);
                    end
                end                
                if any(strcmp(name, {'beta1','beta2','beta3'}))
                    Name = name(1:length(name)-1);
                    b1 = Params.([Name '1']);
                    b2 = Params.([Name '2']);
                    b3 = Params.([Name '3']);
                    Params.(Name) = Params.([Name '_func'])(b1, b2, b3, Vol);
                end                
                
            end
        end
        
        
        % Functions of parameters

        if exist('Qmin_QC_a', 'var') || exist('Qmin_QC_b', 'var') || ...
                exist('Qmax_delQ_a', 'var') || exist('Qmax_delQ_b', 'var')
            Params.Qmax_QC = Params.Qmax_QC_func(Params.Qmin_QC, Params.Qmax_delQ);
            Params.delQ_QC = Params.delQ_QC_func(Params.Qmin_QC, Params.Qmax_QC);
        end
        
        if exist('Vmax_QC_a', 'var') || exist('Vmax_QC_b', 'var') || ... 
                exist('aN_QC_a', 'var') || exist('aN_QC_b', 'var')
            Params.kN = Params.kN_func(Params.Vmax_QC, Params.aN_QC);
        end
        
        if exist('Gmax_a', 'var') || exist('Gmax_b', 'var')
            Params.Gmax = Params.Gmax(FixedParams.zooplankton);
        end
        
        if exist('rDOC', 'var') || exist('rDON', 'var') || ... 
                exist('rPOC', 'var') || exist('rPON', 'var')
            Params.rOM = Params.rOM_func(rDOC, rDON, rPOC, rPON, ...
                FixedParams.DOM_index, FixedParams.POM_index, ... 
                FixedParams.OM_C_index, FixedParams.OM_N_index);
        end
        
        if exist('wDOM', 'var') || exist('wPOM', 'var')
            Params.wk = Params.wk_func(wDOM, wPOM, ... 
                FixedParams.DOM_index, FixedParams.POM_index);
        end
                    
    end

end

end

% function y = powerFunction(a, b, x, varargin)
% y = a .* x .^ b;
% if ~isempty(varargin)
%     tr = strcmp(varargin, 'Transpose');
%     if any(tr) && varargin{find(tr)+1}, y = y'; end
%     ub = strcmp(varargin, 'UpperBound');
%     if any(ub)
%         cap = varargin{find(ub)+1};
%         y(y>cap) = cap;
%     end
% end
% end
% 
% function y = doubleLogisticFunction(a, b, c, x)
% u = exp(x - c);
% y = a ./ (1 + u) .* (1 + b .* u);
% end

function Params = updateParameters(params,FixedParams,varargin)

Params = params;

if ~isempty(varargin)
    if nargin == 3 && isnumeric(varargin{1}) % if new param values are passed as a single vector (with labels stored in FixedParams)
        newvals = varargin{1};
        parnames = FixedParams.tunePars;
        npars = length(parnames);
        V_PP = FixedParams.PPsize;
        for i = 1:npars
            Params.(parnames{i}) = newvals(i); % update scalars
            % update size-dependent parameters - constructed from scalars
            % suffixed with _a and _b
            pn = parnames{i};            
            if ~isempty(regexp(pn,'_a', 'once')) || ~isempty(regexp(pn,'_b', 'once'))
                Name = pn(1:end-2);
                param_a = Params.([Name '_a']);
                param_b = Params.([Name '_b']);
                if ~any(strcmp(Name, {'Qmax_over_delQ'}))
                    Params.(Name) = volumeDependent(param_a, param_b, V_PP);
                else
                    Params.(Name) = 1 ./ (1 - volumeDependent(param_a, param_b, V_PP));
                end
            end
        end
        
        
        if any(strcmp('Qmin_a', parnames)) || any(strcmp('Qmin_b', parnames)) || ...
                any(strcmp('Qmax_over_delQ_a', parnames)) || any(strcmp('Qmax_over_delQ_b', parnames))
            Params.Qmax = Params.Qmin .* (Params.Qmax_over_delQ ./ (Params.Qmax_over_delQ - 1));
            Params.delQ = Params.Qmax - Params.Qmin;
        end
        
        if any(strcmp('wk', parnames)) || any(strcmp('rPOM', parnames))
            % POM sinking and remineralisation matrix
            sinkTime = FixedParams.delz ./ Params.wk; % sinking time of particles moving between centers of adjacent depth layers
            r_x_sinkTime = Params.rPOM .* sinkTime;
            nz = FixedParams.nz;
            dzm = FixedParams.zwidth(2:nz) ./ ...
                (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
            dzp = FixedParams.zwidth(1:nz-1) ./ ...
                (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
            
            POM_to_IN_array = zeros(nz, nz); % (i,j) lower-tri matrix of proportions of POM remineralised while sinking from layer j to i
            POM_to_IN = []; % create a block-diagonal matrix - only needed when modelling multiple nutrients
            for i_nut = 1:FixedParams.nOM
                p = r_x_sinkTime .* tril(cumprod(tril(1 - repmat(r_x_sinkTime, [1 nz-1]), -1) + triu(ones(nz-1,nz-1))));
                POM_to_IN_array(1:nz-1,1:nz-1) = dzp .* p;
                POM_to_IN_array(2:nz,1:nz-1) = POM_to_IN_array(2:nz,1:nz-1) + dzm .* p;
                if ~FixedParams.POM_is_lost
                    % If POM does not sink below bottom depth layer then all that
                    % remains after sinking is remineralised on bottom layer
                    POM_to_IN_array(nz,:) = POM_to_IN_array(nz,:) + (1 - sum(POM_to_IN_array));
                end
                POM_to_IN = blkdiag(POM_to_IN, POM_to_IN_array);
            end
            Params.POM_to_IN = sparse(POM_to_IN); % using sparse matrix is increasingly useful the more nutrient s are modelled...
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
                    if ~any(strcmp(Name, {'Qmax_over_delQ'}))
                        Params.(Name) = volumeDependent(param_a, param_b, V_PP);
                    else
                        Params.(Name) = 1 ./ (1 - volumeDependent(param_a, param_b, V_PP));
                    end
                end
            end
        end
        
        
        if any(strcmp('Qmin_a', varargin)) || any(strcmp('Qmin_b', varargin)) || ...
                any(strcmp('Qmax_over_delQ_a', varargin)) || any(strcmp('Qmax_over_delQ_b', varargin))
            Params.Qmax = Params.Qmin .* (Params.Qmax_over_delQ ./ (Params.Qmax_over_delQ - 1));
            Params.delQ = Params.Qmax - Params.Qmin;
        end
        
        if any(strcmp('wk', varargin)) || any(strcmp('rPOM', varargin))
            % POM sinking and remineralisation matrix
            sinkTime = FixedParams.delz ./ Params.wk; % sinking time of particles moving between centers of adjacent depth layers
            r_x_sinkTime = Params.rPOM .* sinkTime;
            nz = FixedParams.nz;
            dzm = FixedParams.zwidth(2:nz) ./ ...
                (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
            dzp = FixedParams.zwidth(1:nz-1) ./ ...
                (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
            
            POM_to_IN_array = zeros(nz, nz); % (i,j) lower-tri matrix of proportions of POM remineralised while sinking from layer j to i
            POM_to_IN = []; % create a block-diagonal matrix - only needed when modelling multiple nutrients
            for i_nut = 1:FixedParams.nOM
                p = r_x_sinkTime .* tril(cumprod(tril(1 - repmat(r_x_sinkTime, [1 nz-1]), -1) + triu(ones(nz-1,nz-1))));
                POM_to_IN_array(1:nz-1,1:nz-1) = dzp .* p;
                POM_to_IN_array(2:nz,1:nz-1) = POM_to_IN_array(2:nz,1:nz-1) + dzm .* p;
                if ~FixedParams.POM_is_lost
                    % If POM does not sink below bottom depth layer then all that
                    % remains after sinking is remineralised on bottom layer
                    POM_to_IN_array(nz,:) = POM_to_IN_array(nz,:) + (1 - sum(POM_to_IN_array));
                end
                POM_to_IN = blkdiag(POM_to_IN, POM_to_IN_array);
            end
            Params.POM_to_IN = sparse(POM_to_IN); % using sparse matrix is increasingly useful the more nutrient s are modelled...
        end
        
    end
end


function p = volumeDependent(a,b,V)
p = a .* V .^ b;
% end

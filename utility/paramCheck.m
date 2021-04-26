function [Params, Bounds] = paramCheck(params, bounds)
% Make sure Param struct is consistent with the bounds and that every
% parameter has been assigned a value.

Params = params;
Bounds = bounds;

% Have all parameters listed in scalars and sizeDependent been assigned
% values? If not then display warning
allPars = [Params.scalars; Params.sizeDependent];
for i = 1:length(allPars)
    n = allPars{i};
    try p = Params.(n);
    catch, p = nan;
    end
    if any(isnan(p))
        warning(['Value of ' n ' has not been assigned. See defaultParameters.m'])
        Params.(n) = [];
    end
end

% Flag and correct any out-of-bounds parameters. This just catches any
% incompatabilities between values selected for parameters and their bounds.
% Also, if bounds have not been selected then they are set 
for i = 1:length(allPars)
    n = allPars{i};
    p = Params.(n);
    if ~isempty(p)
        try b = Bounds.(n);
        catch, b = nan;
        end
        if ~any(isnan(b)) % if bounds have been set
            if p < b(1)
                warning(['Default ' n ' value is less than its lower bound, so ' ...
                    n ' has been reset to the midpoint between its bounds. Choose consistent values in defaultParameters.m'])
                Params.(n) = 0.5 .* (b(1)+b(2));
            end
            if p > b(2)
                warning(['Default ' n ' value is greater than its upper bound, so ' ...
                    n ' has been reset to the midpoint between its bounds. Choose consistent values in defaultParameters.m'])
                Params.(n) = 0.5 .* (b(1)+b(2));
            end
        else % if bounds have not been set then parameter is deemed to be fixed
            warning(['Parameter ' n ' declared without bounding values. Bounds now set equal to ' n ' value. See defaultParameters.m'])
            Bounds.(n) = [p, p];
        end
    end
end

function [gaOutput, optPar] = storeGAoutput(gapopulationhistory, ... 
    gacosthistory, gaOptions, FixedParams, varargin)

% Store output from genetic algorithm within a singe struct

[I,J] = find(gacosthistory == min(gacosthistory(:))); % extract best parameter set
if length(I) > 1 || length(J) > 1
    % if tied then choose the first best from most recent iteration
    i = find(J == max(J), 1);
    I = I(i);
    J = J(i);
end

optPar = gapopulationhistory(I,:,J);

gaOutput.parNames = FixedParams.tunePars;
gaOutput.lowerBound = FixedParams.tunePars_lb;
gaOutput.upperBound = FixedParams.tunePars_ub;
gaOutput.optPar = optPar;
gaOutput.populationHistory = gapopulationhistory;
gaOutput.scoreHistory = gacosthistory;
gaOutput.gaOptions = gaOptions;

if ~isempty(varargin)
    i = strcmp(varargin, 'stoppedEarly');
    if any(i)
        stoppedEarly = varargin{i+1};
        if ~stoppedEarly
            % If algorithm was terminated early then  some outputs are not available
            gaOutput.fval = fval;
            gaOutput.exitflag = exitflag;
            gaOutput.output = output;
            gaOutput.population = population;
            gaOutput.scores = scores;
        end
    end
end


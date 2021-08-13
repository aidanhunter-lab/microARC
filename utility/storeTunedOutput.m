function [output, optPar] = storeTunedOutput(optimiser, populationhistory, ...
    costhistory, optimiserOptions, FixedParams, varargin)
% Store output from genetic algorithm within a single struct

extractVarargin(varargin)

switch optimiser
    % output may depend upon choice of optimising algorithm
    case 'ga'
        [I,J] = find(costhistory == min(costhistory(:))); % extract best parameter set
        if length(I) > 1 || length(J) > 1
            % if tied then choose the first best from most recent iteration
            i = find(J == max(J), 1);
            I = I(i);
            J = J(i);
        end
        
        for i = 1:length(FixedParams.tunePars)
            % Convert from parameter search space to natural space
            pn = FixedParams.tunePars{i};
            func = FixedParams.tuneParsInvTransform.(pn);
            populationhistory(:,i,:) = func(populationhistory(:,i,:));
        end
        
        optPar = populationhistory(I,:,J);
        
        output.parNames = FixedParams.tunePars;
        output.lowerBound = FixedParams.tunePars_lb;
        output.upperBound = FixedParams.tunePars_ub;
        output.optPar = optPar;
        output.optPar_summary = table(output.parNames', output.lowerBound', ...
            output.optPar', output.upperBound');
        output.optPar_summary.Properties.VariableNames = {'par','lower','opt','upper'};
        output.populationHistory = populationhistory;
        output.scoreHistory = costhistory;
        output.optimiserOptions = optimiserOptions;
        
        % Extra output options may be included if the genetic algorithm was not
        % terminated early
        if ~exist('stoppedEarly', 'var')
            % By default assume that tuning algorithm was halted before completion
            stoppedEarly = true;
        end
        if ~stoppedEarly
            output.fval = evalin('caller', 'fval');
            output.exitflag = evalin('caller', 'exitflag');
            output.output = evalin('caller', 'output');
            output.population = evalin('caller', 'population');
            output.scores = evalin('caller', 'scores');
        end

end


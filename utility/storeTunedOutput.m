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
        
%         for i = 1:length(FixedParams.tunePars)
%             % Convert from parameter search space to natural space
%             pn = FixedParams.tunePars{i};
%             func = FixedParams.tuneParsInvTransform.(pn);
%             populationhistory(:,i,:) = func(populationhistory(:,i,:));
%         end
        

        optPar = populationhistory(I,:,J);
        optPar_raw = nan(size(optPar)); % parameter output on natural scale (rather than search-space scale)
        
        % The 'output' struct stores parameter info on search-space scale
        % so that it may easily be used to restart an optimisation from
        % prior run.
        % The only fields where parameter are stored in their natural scale
        % are optPar and optPar_summary.
        for i = 1:length(optPar)
            % Convert from parameter search space to natural space
            pn = FixedParams.tunePars{i};
            func = FixedParams.tuneParsInvTransform.(pn);
            optPar_raw(i) = func(optPar(i));
        end
        
        output.parNames = FixedParams.tunePars;
        output.lowerBound = optimiserOptions.InitialPopulationRange(1,:);
        output.upperBound = optimiserOptions.InitialPopulationRange(2,:);
        output.optPar_searchSpace = optPar;
        output.optPar = optPar_raw;
        output.optPar_summary = table(output.parNames', FixedParams.tunePars_lb', ...
            output.optPar', FixedParams.tunePars_ub');
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


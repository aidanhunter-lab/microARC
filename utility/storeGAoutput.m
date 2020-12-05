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
gaOutput.optPar_summary = table(gaOutput.parNames', gaOutput.lowerBound', ... 
    gaOutput.optPar', gaOutput.upperBound');
gaOutput.optPar_summary.Properties.VariableNames = {'par','lower','opt','upper'};
gaOutput.populationHistory = gapopulationhistory;
gaOutput.scoreHistory = gacosthistory;
gaOutput.gaOptions = gaOptions;


% Extra output options may be included if the genetic algorithm was not
% terminated early

if ~isempty(varargin)
    i = strcmp(varargin, 'stoppedEarly');
    if any(i)
        stoppedEarly = varargin{find(i)+1};
    else
        stoppedEarly = false;
    end
    if ~stoppedEarly
        gaOutput.fval = evalin('base', 'fval');
        gaOutput.exitflag = evalin('base', 'exitflag');
        gaOutput.output = evalin('base', 'output');
        gaOutput.population = evalin('base', 'population');
        gaOutput.scores = evalin('base', 'scores');
    end
end


% if ~isempty(varargin)
%     i = strcmp(varargin, 'stoppedEarly');
%     if any(i)
%         stoppedEarly = varargin{i+1};
%         if ~stoppedEarly
%             % If algorithm was terminated early then  some outputs are not available
%             gaOutput.fval = fval;
%             gaOutput.exitflag = exitflag;
%             gaOutput.output = output;
%             gaOutput.population = population;
%             gaOutput.scores = scores;
%         end
%     end
% end
% 

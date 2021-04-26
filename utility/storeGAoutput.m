function [gaOutput, optPar] = storeGAoutput(gapopulationhistory, ... 
    gacosthistory, gaOptions, FixedParams, varargin)
% Store output from genetic algorithm within a single struct

extractVarargin(varargin)

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
if ~exist('stoppedEarly', 'var')
    % By default assume that tuning algorithm was halted before completion
    stoppedEarly = true;
end

if ~stoppedEarly
    gaOutput.fval = evalin('caller', 'fval');
    gaOutput.exitflag = evalin('caller', 'exitflag');
    gaOutput.output = evalin('caller', 'output');
    gaOutput.population = evalin('caller', 'population');
    gaOutput.scores = evalin('caller', 'scores');
end


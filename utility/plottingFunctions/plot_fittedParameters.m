function plt = plot_fittedParameters(outputSummary)

% Display the fitted parameters relative to their bounds

plt = figure;
plt.Units = 'inches';
plt.Position = [0 0 12 12];

%%

subplot(1,2,1)

outputSummary = flip(outputSummary);

% fitted parameter values and bounds
parName = outputSummary.par;
boundLower = outputSummary.lower;
boundUpper = outputSummary.upper;
parValue = outputSummary.opt;
% rescale within (0,1)
parValueScaled = (parValue - boundLower) ./ (boundUpper - boundLower);

npars = length(parValue);

% Adjust some parameter names for plotting
x = regexp(parName, '_');
for i = 1:npars
    if length(x{i}) == 2        
        xx = x{i};
        pn = parName{i};
        pn(xx(1)) = '/';
        pn_ = ['*' pn(1:xx(2)-1) '*' pn(xx(2):end)];
        x_ = regexp(pn_, '*');
        pn_(x_(1)) = '(';
        pn_(x_(2)) = ')';
        parName{i} = pn_;
    end
end

y = 1:npars;


colOK = [0 1 0];
colBad = [1 0 0];
colBorderline = [1 0.5 0];

Cols = repmat(colOK, [npars, 1]);

for i = 1:npars
    if parValueScaled(i) < 0.05 || parValueScaled(i) > 0.95
       Cols(i,:) = colBorderline;
    end
    if parValueScaled(i) < 0.01 || parValueScaled(i) > 0.99
       Cols(i,:) = colBad;
    end
end

for i = 1:npars
    scatter(parValueScaled(i), y(i), 'MarkerFaceColor', Cols(i,:), 'MarkerEdgeColor', Cols(i,:))
    if i == 1, hold on; end
end

for i = 1:npars
    plot([0 1], [y(i) y(i)], ':', 'Color', [0.5 0.5 0.5])
end

%%

gc = gca;
gc.YLim = [0 npars+1];
gc.YTick = y;
gc.YTickLabel = parName;
gc.TickLength = [0 0];
gc.XTick = [0.05 0.25 0.5 0.75 0.95];

plot([0.05 0.05], gc.YLim, '--', 'Color', [0.5 0.5 0.5])
plot([0.95 0.95], gc.YLim, '--', 'Color', [0.5 0.5 0.5])

title('Fitted parameter values relative to their bounds')

hold off

%%
sp = subplot(1,2,2);
pos = get(sp,'Position');
un = get(sp,'Units');
delete(sp)

outputSummary_ = outputSummary;

outputSummary_.Properties.VariableNames = {'Parameter','Lower','Value','Upper'};
outputSummary_.Properties.RowNames = outputSummary_.Parameter;
outputSummary_.Parameter = [];
outputSummary_ = flip(outputSummary_);

uitable('Data',outputSummary_{:,:},'ColumnName',outputSummary_.Properties.VariableNames,...
    'RowName',outputSummary_.Properties.RowNames,'Units', un', 'Position', pos);



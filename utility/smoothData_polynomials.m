function out = smoothData_polynomials(Data, maxDegree_scalar, maxDegree_size, varargin)

% Scalar data
obs = Data.scalar.obsInCostFunction;
% loop through data types
for i = 1:length(obs)
    ind = strcmp(Data.scalar.Variable, obs{i});
    y = Data.scalar.scaled_Value(ind);
    n = length(y);
    x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
    x = x ./ max(x); % evaluate polynomial over [-1,1]
    [y, yo] = sort(y);
    alpha = polyfit(x, y, maxDegree_scalar); % fit polynomial coefficients
    Data.scalar.(['polyXvals_' obs{i}]) = x;
    Data.scalar.(['polyCoefs_' obs{i}]) = alpha;
    Data.scalar.(['polyWeights_' obs{i}]) = coefWeights(maxDegree_scalar);
    Data.scalar.(['sortOrder_' obs{i}]) = yo;
    ip = ~isempty(varargin) && any(strcmp(varargin, ['plotPoly_' obs{i}]));
    if ip
        makePlot = varargin{find(strcmp(varargin, ['plotPoly_' obs{i}]))+1};
        if makePlot
            plotSmooth(x, y, alpha, obs{i});
        end
    end
end


% Size spectra
obs = unique(Data.size.dataBinned.Variable);
for i = 1:length(obs)
    ind = strcmp(Data.size.dataBinned.Variable, obs{i});
    y = Data.size.dataBinned.scaled_Value(ind);    
    n = length(y);
    x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
    x = x ./ max(x); % evaluate polynomial over [-1,1]
    [y, yo] = sort(y);
    alpha = polyfit(x, y, maxDegree_size);
    Data.size.(['polyXvals_' obs{i}]) = x;
    Data.size.(['polyCoefs_' obs{i}]) = alpha;
    Data.size.(['polyWeights_' obs{i}]) = coefWeights(maxDegree_size); % default weights for polynomial coefs in cost function
    Data.size.(['sortOrder_' obs{i}]) = yo;    
    ip = ~isempty(varargin) && any(strcmp(varargin, ['plotPolySize_' obs{i}]));
    if ip
        makePlot = varargin{find(strcmp(varargin, ['plotPolySize_' obs{i}]))+1};
        if makePlot
            plotSmooth(x, y, alpha, obs{i});
        end
    end
end


out = Data;

end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Each polynomial coefficient is ascribed a weighting in the cost function.
% Use coefWeights to specify these weights. Current function gives equal
% weight to order 0 and order 1 terms (mean and variance), and downweights
% higher order terms by scaling by the factorial of their order (the
% pattern follows a Taylor series).
function w = coefWeights(degree)
    w = 1 ./ factorial(degree:-1:0);
end

function fig = plotSmooth(x, y, alpha, var)
    fig = figure;
    ymod = polyval(alpha, x);
    plot(x, y, 'o', x, ymod, '-')
    scatter(x, y, 'MarkerEdgeColor', [0 0 1])
    hold on    
    plot(x, ymod, 'Color', [1 0 0])
    xlabel('x')
    ylabel(var)
    gc = gca;
    xl = gc.XLim; yl = gc.YLim;
    text(xl(1)+0.1*diff(xl), yl(2)-0.05*diff(yl), 'ordered standardised data')
    text(xl(1)+0.1*diff(xl), yl(2)-0.1*diff(yl), 'polynomial smoother')    
    scatter(xl(1)+0.05*diff(xl), yl(2)-0.05*diff(yl), 'MarkerEdgeColor', [0 0 1])
    line([xl(1)+0.025*diff(xl), xl(1)+0.075*diff(xl)], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], 'Color', [1 0 0])
    coefs = flip(round(alpha, 2, 'significant'));
    eqPoly = ['y = ' num2str(coefs(1))];
    if length(coefs) > 1
        for i = 2:length(coefs)
            o = i-1;
            if o == 1
                if coefs(i) == abs(coefs(i))
                    eqPoly = [eqPoly ' + ' num2str(coefs(i)) 'x'];
                else
                    eqPoly = [eqPoly ' - ' num2str(abs(coefs(i))) 'x'];
                end
            else
                if coefs(i) == abs(coefs(i))
                    eqPoly = [eqPoly ' + ' num2str(coefs(i)) 'x^' num2str(o)];
                else
                    eqPoly = [eqPoly ' - ' num2str(abs(coefs(i))) 'x^' num2str(o)];
                end
            end
        end
    end
    text(xl(1)+0.1*diff(xl), yl(2)-0.15*diff(yl), eqPoly);
    hold off
end


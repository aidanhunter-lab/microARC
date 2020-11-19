function out = smoothData_polynomials(Data, maxDegree_scalar, maxDegree_size, varargin)

v = reshape(varargin, [2 0.5*length(varargin)]);


% N
y = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, 'N'));
n = length(y);
x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
x = x ./ max(x); % evaluate polynomial over [-1,1]
[y, yo] = sort(y);
alpha = polyfit(x, y, maxDegree_scalar);
Data.scalar.polyXvals_N = x;
Data.scalar.polyCoefs_N = alpha;
Data.scalar.sortOrder_N = yo;
if ~isempty(v) && any(contains(v(1,:),'plotPolyN')) && v{2,strcmp(v(1,:), 'plotPolyN')}
    plotSmooth(x, y, alpha, 'N');
end

% PON
y = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, 'PON'));
n = length(y);
x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
x = x ./ max(x); % evaluate polynomial over [-1,1]
[y, yo] = sort(y);
alpha = polyfit(x, y, maxDegree_scalar);
Data.scalar.polyXvals_PON = x;
Data.scalar.polyCoefs_PON = alpha;
Data.scalar.sortOrder_PON = yo;
if ~isempty(v) && any(contains(v(1,:),'plotPolyPON')) && v{2,strcmp(v(1,:), 'plotPolyPON')}
    plotSmooth(x, y, alpha, 'PON');
end

% POC
y = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, 'POC'));
n = length(y);
x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
x = x ./ max(x); % evaluate polynomial over [-1,1]
[y, yo] = sort(y);
alpha = polyfit(x, y, maxDegree_scalar);
Data.scalar.polyXvals_POC = x;
Data.scalar.polyCoefs_POC = alpha;
Data.scalar.sortOrder_POC = yo;
if ~isempty(v) && any(contains(v(1,:),'plotPolyPOC')) && v{2,strcmp(v(1,:), 'plotPolyPOC')}
    plotSmooth(x, y, alpha, 'POC');
end


% Chl
y = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, 'chl_a'));
n = length(y);
x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
x = x ./ max(x); % evaluate polynomial over [-1,1]
[y, yo] = sort(y);
alpha = polyfit(x, y, maxDegree_scalar);
Data.scalar.polyXvals_chl_a = x;
Data.scalar.polyCoefs_chl_a = alpha;
Data.scalar.sortOrder_chl_a = yo;
if ~isempty(v) && any(contains(v(1,:),'plotPolyChl')) && v{2,strcmp(v(1,:), 'plotPolyChl')}
    plotSmooth(x, y, alpha, 'Chl');
end



% Size spectra
y = Data.size.dataBinned.scaled_Ntot;
n = length(y);
x = (-0.5*(n-1):1:0.5*(n-1))'; % center curve over origin for intercept near 0
x = x ./ max(x); % evaluate polynomial over [-1,1]
[y, yo] = sort(y);
alpha = polyfit(x, y, maxDegree_size);
Data.size.polyXvals = x;
Data.size.polyCoefs = alpha;
Data.size.sortOrder = yo;
if ~isempty(v) && any(contains(v(1,:),'plotPolySize')) && v{2,strcmp(v(1,:), 'plotPolySize')}
    plotSmooth(x, y, alpha, 'N-at-size');
end


out = Data;

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
    if strcmp(var, 'N-at-size')
        text(xl(1)+0.1*diff(xl), yl(2)-0.05*diff(yl), 'ordered data')
    else
        text(xl(1)+0.1*diff(xl), yl(2)-0.05*diff(yl), 'ordered standardised data')
    end    
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

% function fig = plotSmooth_size(x, y, ymod, var, Data)
%     fig = figure;
%     plot(x, y, 'Color', [0 0 0]);
%     xlabel('x')
%     ylabel(var)
%     hold on
%     plot(x, ymod, 'Color', [1 0 1])
%     gc = gca;
%     xl = gc.XLim; yl = gc.YLim;
%     text(xl(1)+0.1*diff(xl), yl(2)-0.05*diff(yl), 'ordered data')
%     text(xl(1)+0.1*diff(xl), yl(2)-0.1*diff(yl), 'polynomial smoother')
%     line([xl(1)+0.025*diff(xl), xl(1)+0.075*diff(xl)], [yl(2)-0.05*diff(yl), yl(2)-0.05*diff(yl)], 'Color', [0 0 0])
%     line([xl(1)+0.025*diff(xl), xl(1)+0.075*diff(xl)], [yl(2)-0.1*diff(yl), yl(2)-0.1*diff(yl)], 'Color', [1 0 1])
%     coefs = round(Data.size.polyCoefs, 2, 'significant');
%     eqPoly = ['y = ' num2str(coefs(1))];
%     if length(coefs) > 1
%         for i = 2:length(coefs)
%             o = i-1;
%             if o == 1
%                 if coefs(i) == abs(coefs(i))
%                     eqPoly = [eqPoly ' + ' num2str(coefs(i)) 'x'];
%                 else
%                     eqPoly = [eqPoly ' - ' num2str(abs(coefs(i))) 'x'];
%                 end
%             else
%                 if coefs(i) == abs(coefs(i))
%                     eqPoly = [eqPoly ' + ' num2str(coefs(i)) 'x^' num2str(o)];
%                 else
%                     eqPoly = [eqPoly ' - ' num2str(abs(coefs(i))) 'x^' num2str(o)];
%                 end
%             end
%         end
%     end
%     text(xl(1)+0.1*diff(xl), yl(2)-0.15*diff(yl), eqPoly);
%     hold off
% end

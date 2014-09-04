function fit = twoexponentialfit(d)
%fit = twoexponentialfit(varargin)
%param = [Atrip D0 A0 K k1 k2 c]
%
%problem = twoexponentialfit(problem) - Adds general p0, ub, and lb fields to
%                                       structure input, problem.
[nRows, nCols] = size(d.tData);
nTraces = nCols/2;

if ~isfield(d, 'p0')

    p0i = [1 0.1 1000 100 0];
    ubi = [10 1 inf inf 0.01];
    lbi = [0 0 0 0 -0.01];
    
    d.p0 = zeros(nTraces, 5);
    d.lb = zeros(nTraces, 5);
    d.ub = zeros(nTraces, 5);
    
    for i = 1:nTraces
        d.p0(i, :) = p0i;
        d.ub(i, :) = ubi;
        d.lb(i, :) = lbi;
        
        fit = d;
    end %for
    return
end %if

%Verify that matrices are the correct sizes
lb_size   = size(d.lb);
ub_size   = size(d.ub);
p_size    = size(d.p0);

if nTraces ~= lb_size(1)
    msgbox('lb wrong size')
elseif nTraces ~= ub_size(1)
    msgbox('ub wrong size')
elseif nTraces ~= p_size(1)
    msgbox('p0 wrong size')
end %if

%Divide traces into time and y

xData = zeros(nRows, nTraces);
yData = zeros(nRows, nTraces);

for i = 1:nTraces
    xData(:,i) = d.tData(:, i*2-1);
    yData(:,i) = d.tData(:, i*2);
end %for

fitfun  = @(p, t) twoexponentialfun(p, t, d.g);
%options = optimset('MaxFunEvals', 90000);
[pFit, ~, resi, ~, ~, ~, J] = lsqcurvefit(fitfun, d.p0, xData, yData, d.lb, d.ub);

%Error analysis
ciTemp = nlparci(pFit, resi, 'jacobian', J);
for i = 1:p_size(2);
    ci(:, i*2-1:i*2) = ciTemp(i*nTraces-(nTraces-1):i*nTraces, 1:2);
end %for

yFit    = fitfun(pFit, xData);

fit.yFit = yFit;
fit.pFit = pFit;
fit.ci = ci;

% Create figure
figure
hold on

for i = 1:nTraces
    plot(xData(:, i), yData(:, i), 'kx')
    plot(xData(:, i), yFit(:, i), 'r-')
end %for

set(gca, 'XScale', 'log')
end %twoexponentialfit
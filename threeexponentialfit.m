function fit = threeexponentialfit(d)
%fit = threeexponentialfit(varargin)
%param = 
%
%problem = twoexponentialfit(problem) - Adds general p0, ub, and lb fields to
%                                       structure input, problem.

[nRows, nCols] = size(d.tData); %Triplet data format alternating t and y cols.
nTraces = nCols/2;

if ~isfield(d, 'p0')

    p0i = [0.999 0.3 0.3 0]; %p = [aTrip f1 f2 c]
    ubi = [1 1 1 0.01];
    lbi = [0.95 0 0 -0.01];
    
    d.p0 = zeros(nTraces, 4);
    d.lb = zeros(nTraces, 4);
    d.ub = zeros(nTraces, 4);
    
    for i = 1:nTraces
        d.p0(i, :) = p0i;
        d.ub(i, :) = ubi;
        d.lb(i, :) = lbi;
        
    end %for
    d.fixedParam = [85 220 1300];
    d.fixLb      = [80 200 400];
    d.fixUb      = [90 300 10000];
    
    fit = d; %Breaks and returns the input with starting parameters.
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

%Vectorize parameters and constraints

[p, nFixedP, nCols] = xformparam(d.p0, d.fixedParam);
[lb, ~, ~] = xformparam(d.lb, d.fixLb);
[ub, ~, ~] = xformparam(d.ub, d.fixUb);

extras = [nFixedP nCols];

fitfun  = @(p, t) threeexponentialfun(p, t, extras);
options = optimset('MaxFunEvals', 90000);
[pFit, ~, resi, ~, ~, ~, J] = lsqcurvefit(fitfun, p, xData, yData, lb, ub, options);

%Error analysis
ciTemp = nlparci(pFit, resi, 'jacobian', J);

%Reconfigure ci
fixCi = ciTemp(1:3, :);
ciTemp = ciTemp(4:end, :);

ci = zeros(nTraces, p_size*2);

for i = 1:p_size(2);
    for j = 1:nTraces
        ci(j, i*2-1:i*2) = ciTemp((j-1)*p_size(2)+i, 1:2);
    end %for
end %for

yFit    = fitfun(pFit, xData);

%Devectorize parameters

[pFit, fixedParamFit] = xformparam(pFit, nFixedP, nCols);
fit.yFit = yFit;
fit.pFit = pFit;
fit.fixedParamFit = fixedParamFit;
fit.ci = ci;
fit.fixCi = fixCi;

% Create figure
figure
hold on

for i = 1:nTraces
    plot(xData(:, i), yData(:, i), 'kx')
    plot(xData(:, i), yFit(:, i), 'r-')
end %for

set(gca, 'XScale', 'log')
end %twoexponentialfit
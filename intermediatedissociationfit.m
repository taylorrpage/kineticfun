function fit = intermediatedissociationfit(d)

if ~isstruct(d)
    error('Nothing:atall', ['Input must be a structure with the following fields:\n' ...
           'xData\nyData\nconstants\nconcentrations\np0\nlb\nub\n' ...
           'globalAttr'])
end %if

%Verify inputs
%--------------------------
equalTraces = checkRows(d.xData', d.yData', d.constants, d.concentrations, ...
    d.p0, d.lb, d.ub, d.globalAttr);
if ~equalTraces
    error(['Mismatched inputs sizes. Check that number of all parameter rows ' ...
           'equals number of all data columns'])
end %if
equalPoints = checkColumns(d.xData, d.yData);
if ~equalTraces
    error('xData and yData have different number of points')
end %if



%Assign gloal paramters NaN
p0 = d.p0;
p0(d.globalAttr == 1) = NaN;

%Global fit
%--------------------------
fun = @(p0, xData) intermediatedissociationfun(p0, xData, d.constants, ...
    d.concentrations, d.globalAttr);
[pFit, ~, resi, ~, ~, ~, J] = lsqcurvefit(fun, d.p0, d.xData, d.yData, d.lb, ...
    d.ub);
yFit = intermediatedissociationfun(pFit, d.xData, ...
    d.constants, d.concentrations, d.globalAttr);

%Error analysis
ci = nlparci(pFit, resi, 'jacobian', J);

%Graph with fits
%--------------------------
figure
hold on

[nTraces, ~] = size(d.p0);
for i = 1:nTraces
    plot(d.xData(:,i), d.yData(:,1), 'ko')
end

for i = 1:nTraces
    plot(d.xData(:,i), yFit(:,1), 'r-')
end

set(gca, 'XScale', 'log')

%Plot non-shared parameters with 95% confidence intervals.
%Report shared parameters with 95% confidence intervales.

%Assemble struct
%--------------------------
fit.pFit = pFit;
fit.yFit = yFit;
fit.ci   = ci;
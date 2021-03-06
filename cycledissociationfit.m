function fit = cycledissociationfit(d)
%fit = cycledissociationfit(d)
%
%Fits problem, d to cylcedissociationfun and returns plot of triplets and
%intermediates with fits.
%
%d     Structure containing the following fields
%   tData - Matrix, n columns wide of triplet data with paired time data and
%           signal data.
%   iData - Matrix, n columns wide of intermediate data with paired time data
%           and signal data.
%   constants - [log(Ka) phi]; Matrix with n*2-1 rows.
%   concentrations - [zntAdded(uM) fe3Added(uM) fe2Added(uM)]; Matrix with n*2-1
%                    row.
%   globalAttr - Matrix, n*2-1 rows, indicating which parameters in p0 are
%                global. A value of 0 indicates a unique paramter. A value of 1
%                will copy the parameter value from row 1 of the same column.
%
%See cycledissociationfun for more information about the kinetic model.

if ~isstruct(d)
    error('Nothing:atall', ['Input must be a structure with the following fields:\n' ...
           'tData\niData\nconstants\nconcentrations\np0\nlb\nub\n' ...
           'globalAttr'])
end %if

%Stack triplet and intermediate data
%--------------------------
tSize = size(d.tData); %Collect triplet data size to manipulate the stacked data.
iSize = size(d.iData);

stackData = [d.tData; d.iData]; %Stack data into single matrix for lsqcurvefit.
stackXData = zeros(tSize(1) + iSize(1), tSize(2)/2);
stackYdata = zeros(tSize(1) + iSize(1), tSize(2)/2);

for i = 1:tSize(2)/2
    stackXData(:, i) = stackData(:, i*2-1); 
    stackYData(:, i) = stackData(:, i*2);
end

s = tSize(1); %Pass to cycledissociationfun to split triplet and intermediate.

%Verify inputs
%--------------------------
equalTraces = checkRows(stackXData', stackYData', d.constants, d.concentrations, ...
    d.p0, d.lb, d.ub, d.globalAttr);
if ~equalTraces
    error(['Mismatched inputs sizes. Check that number of all parameter rows ' ...
           'equals number of all data columns'])
end %if

%Assign gloal paramters NaN
p0 = d.p0;
p0(d.globalAttr == 1) = NaN;

%Global fit
%--------------------------
fun = @(p0, xData) cycledissociationfun(p0, xData, d.constants, ...
    d.concentrations, d.globalAttr, s);
[pFit, ~, resi, ~, ~, ~, J] = lsqcurvefit(fun, p0, stackXData, stackYData, ...
    d.lb, d.ub);
yFit = cycledissociationfun(pFit, stackXData, ...
    d.constants, d.concentrations, d.globalAttr, s);

tYFit = yFit(1:s, :);
iYFit = yFit(s+1:end, :);

%Error analysis
ci = nlparci(pFit, resi, 'jacobian', J);

%Graph with fits
%--------------------------

[nTraces, ~] = size(d.p0);
subplot(1,2,1);
hold on
for i = 1:nTraces
    plot(d.tData(:,i*2-1), d.tData(:,i*2), 'ko')
end

subplot(1,2,1);
for i = 1:nTraces
    plot(d.tData(:,i*2-1), tYFit(:,i), 'r-')
end

set(gca, 'XScale', 'log')

subplot(1,2,2);
hold on
for i = 1:nTraces
    plot(d.iData(:,i*2-1), d.iData(:,i*2), 'ko')
end

subplot(1,2,2);
for i = 1:nTraces
    plot(d.iData(:,i*2-1), iYFit(:,i), 'r-')
end

set(gca, 'XScale', 'log')

%Plot non-shared parameters with 95% confidence intervals.
%Report shared parameters with 95% confidence intervales.

%Assemble struct
%--------------------------
fit.pFit = pFit;
fit.tYFit = tYFit;
fit.iYFit = iYFit;
fit.ci   = ci;
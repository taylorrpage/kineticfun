function fit = rdfit(d)
%fit = rdfit(d)
%
%Fits problem, d to rdfun and returns plot of triplets and
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
%   p0 = {localParam, globalParam}
%       localParam - Matrix of initial local paramters.
%                    [aTrip aInt aJump kD kQ c], each row corresponds to one 
%                    trace.
%       globalParam - Vector of global parameters.
%                     [kbA; kbB; kbC; kbD; konA; konB; konC; konD; koffM; koffB].
%                     kB, kOn, and kOff are calculated using the following
%                     equations.
%                     kB   = kbA/(kbB*(x - kbC)) + kbD, x = fe3Added
%                     kOn  = konA/(konB*(x - konC)) + konD, x = fe3Added
%                     kOff = koffM*x + koffB, x = fe3Added
%   lb = {lbLocal, lbGlobal}
%   ub = {ubLocal, ubGlobal}
%   globalAttr - Matrix, n*2-1 rows, indicating which parameters in p0 are
%                global. A value of 0 indicates a unique paramter. A value of 1
%                will copy the parameter value from row 1 of the same column.
%
%See rdfun for more information about the kinetic model.

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
    d.p0{1}, d.lb{1}, d.ub{1}, d.globalAttr);
if ~equalTraces
    error(['Mismatched inputs sizes. Check that number of all parameter rows ' ...
           'equals number of all data columns'])
end %if

%Vectorize parameters and bounds
%--------------------------
localParam0  = d.p0{1};
globalParam0 = d.p0{2};

localParam0(d.globalAttr == 1) = NaN;

lbLocal = d.lb{1};
lbGlobal = d.lb{2};
ubLocal = d.ub{1};
ubGlobal = d.ub{2};

[p0Vector, nFixedP, nCols] = xformparam(localParam0, globalParam0);
[lbVector, ~, ~] = xformparam(lbLocal, lbGlobal);
[ubVector, ~, ~] = xformparam(ubLocal, ubGlobal);

%Assemble extras
%--------------------------
extras = [s nFixedP nCols]; %Should come from vectorization of localParam0.

%Global fit
%--------------------------
fun = @(p0, xData) rdfun(p0, xData, d.constants, ...
    d.concentrations, d.globalAttr, extras);

options = optimset('MaxFunEvals', 9000);
[pFit, ~, resi, ~, ~, ~, J] = lsqcurvefit(fun, p0Vector, stackXData, stackYData, ...
    lbVector, ubVector, options);
yFit = rdfun(pFit, stackXData, ...
    d.constants, d.concentrations, d.globalAttr, extras);

tYFit = yFit(1:s, :);
iYFit = yFit(s+1:end, :);

%Error analysis
%--------------------------
ci = nlparci(pFit, resi, 'jacobian', J);

%Reform pFit and to match input
%--------------------------
[pFitLocal, pFitGlobal] = xformparam(pFit, nFixedP, nCols);
pFit = {pFitLocal, pFitGlobal};

%Graph with fits
%--------------------------

[nTraces, ~] = size(d.p0{1});

figure; hold on;

subplot(1,2,1);
hold on
for i = 1:nTraces
    plot(d.tData(:,i*2-1), d.tData(:,i*2), 'ko')
end

for i = 1:nTraces
    plot(d.tData(:,i*2-1), tYFit(:,i), 'r-')
end

set(gca, 'XScale', 'log')

subplot(1,2,2);
hold on
for i = 1:nTraces
    plot(d.iData(:,i*2-1), d.iData(:,i*2), 'ko')
end

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
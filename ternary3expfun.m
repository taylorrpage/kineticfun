function Y = ternary3expfun(param, t, concentrations, g, extras)
%Y = ternary3expfun(param, t, g)

nFixedP = extras(1); %Fixed constants, K1, K2, k1, k2, k3
nCols = extras(2);

[param, fixedParam] = xformparam(param, nFixedP, nCols);

validInputs = checkRows(param, t', g);
if ~validInputs
    error('Param rows and t cols not equal')
end %if

param(g == 1) = NaN; %Uses g to determine which parameters are global.
param = populateconstants(param); %Replace NaN for global parameters.

[~, nTraces] = size(t);

K1 = fixedParam(1);
K2 = fixedParam(2);
k1 = fixedParam(3);
k2 = fixedParam(4);
k3 = fixedParam(5);

d0 = concentrations(1);
a0 = concentrations(2:end);

aTrip = param(:,1);
c     = param(:,2);

[d, ~, da, daa] = bindingtwosite(K1, K2, d0, a0);

totalD = d + da + daa;
f1 = d./totalD;
f2 = da./totalD;
f3 = 1 - f1 - f2;

Y = zeros(size(t));
for i = 1:nTraces
Y(:, i) = aTrip(i) .* (f1(i).*exp(-k1.*t(:, i)) + ...
                       f2(i).*exp(-k2.*t(:, i)) + ...
                       f3(i).*exp(-k3.*t(:, i))) + c(i);
end %for
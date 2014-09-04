function Y = twoexponentialfun(param, t, g)
%Y = twoexponentialfun(param, t, g)

[nRows, nCols] = size(t);

validInputs = checkRows(param, t');
if ~validInputs
    error('Param rows and t cols not equal')
end %if

param(g == 1) = NaN;
param = populateconstants(param); %Replace NaN for global parameters.

aTrip = param(:, 1);
f     = param(:, 2);
k1    = param(:, 3);
k2    = param(:, 4);
c     = param(:, 5);

Y = zeros(nRows, nCols);

for i = 1:nCols
    Y(:, i) = aTrip(i) .* ((1-f(i)).*exp(-k1(i).*t(:, i)) + ...
                            f(i).*exp(-(k2(i)).*t(:, i))) + c(i);
end %for

end %twoexponentialfun
function Y = threeexponentialfun(param, t, extras)
%Y = threeexponentialfun(param, t g, extras)

nFixedP = extras(1);
nCols   = extras(2);

[param, fixedParam] = xformparam(param, nFixedP, nCols);

validInputs = checkRows(param, t');
if ~validInputs
    error('Param rows and t cols not equal')
end %if

[~, nTraces] = size(t);

aTrip = param(:,1);
f1    = param(:,2);
f2    = param(:,3);
c     = param(:,4);

k1 = fixedParam(1);
k2 = fixedParam(2);
k3 = fixedParam(3);

Y = zeros(size(t));
for i = 1:nTraces
Y(:, i) = aTrip(i) .* (f1(i).*exp(-k1.*t(:, i)) + ...
                       f2(i).*exp(-k2.*t(:, i)) + ...
                       (1-f1(i)-f2(i)).*exp(-k3.*t(:, i))) + c(i);
end %for
function Y = rdfun(param, t, constants, concentrations, g, extras)
%Y = rdfun(param, t, constants, concentrations, g)
%
%REDUCED DISSOCIATION
%
%Simulates intermediate data for a simple dissociation model including
%dissociation of the intermediate, but with strictly slow exchange triplet
%decay. kB, kOn, and kOff are constrained.
%
%param          = [globalParam; localParam]. Vectorized parameters from
%                 xformparam. Pre-vectorized parameters:
%                 globalParam = [kbA; kbB; kbC; kbD; konA; konB; konC; konD; ...
%                               koffM; koffB].
%                 localParam  = [aTrip; aInt; aJump; kD; kQ; c]
%t              - Vertical vector or array of time data where each column
%                 corresponds to one trace.
%constants      = [log(Ka) phi]; Vector or array where each row corresponds to
%                 one trace. Ka is the association constant in M-1. phi is the
%                 quantum yield of triplet. This is necessary for determining
%                 initial concentrations of the excited state species for second
%                 order reactions.
%concentrations = [zntAdded(uM) fe3Added(uM) fe2Added(uM)]; Vector or array
%                 where each row corresponds to one trace. Concentrations must
%                 be in micromolar.
%g              - Vector or matrix with the same dimensions as param which
%                 designates global parameters. A value of 0 indicates a unique
%                 paramter. A value of 1 will copy the parameter value from row
%                 1 of the same column.
%extras         - [s nFixedP, nCols] Vector.

s = extras(1);
nFixedP = extras(2);
nCols = extras(3);

tT = t(1:s, :); %Triplet time vector
tI = t(s+1:end, :); %Intermediate time vector

tTSize = size(tT);
tISize = size(tI);

nTraces = tTSize(2);

[param, fixedParam] = xformparam(param, nFixedP, nCols);

%Error if mismatched arrays
validInputs = checkRows(param, constants, concentrations, tT', tI');

if ~validInputs
    error('Param rows, t cols, constants rows, and concentrations row must be equal')
end %if

[kbA, kbB, kbC, kbD, konA, konB, konC, konD, koffM, koffB] = deal(...
    fixedParam(1), fixedParam(2), fixedParam(3), fixedParam(4), fixedParam(5),...
    fixedParam(6), fixedParam(7), fixedParam(8), fixedParam(9), fixedParam(10));

Y = zeros(tTSize(1) + tISize(1), tTSize(2));

param(g == 1) = NaN;
param = populateconstants(param); %Replace NaN for global parameters.
    
for i = 1:nTraces
    
    %Parse parameters
    %----------------

    [aTrip, aInt, aJump, kD, kQ, c] = ...
    deal(param(i, 1), param(i, 2), param(i, 3), param(i, 4), param(i, 5), ...
         param(i, 6));
     
     kB = kbA/(kbB*(concentrations(i, 2) - kbC)) + kbD;
     kOn = konA/(konB*(concentrations(i, 2) - konC)) + konD;
     kOff = koffM*concentrations(i, 2) + koffB;

    zntAdded = concentrations(i, 1);
    fe3Added = concentrations(i, 2);
    fe2Added = concentrations(i, 3);
    Ka       = constants(i, 1);
    phi      = constants(i, 2);

    Ka       = 10^Ka;
    
    %Calculate preequibrium concentrations
    %----------------
    [znt, ~, zntFe3] = calculatebinding_onesite(zntAdded, fe3Added, Ka);
    
    %Account for quantum yield
    %----------------
    zng    = (1 - phi) * znt; %zng calculations before znt or phi is squared.
    zngFe3 = (1 - phi) * zntFe3;
    znt    = phi * znt;
    zntFe3 = phi * zntFe3;
    zni    = 0;
    zniFe2 = 0;

    %Check if t(1) is good estimate for t(0), add new t(0) if needed.
    %----------------
    [tTi, tripletHasAddedT] = checkTi(kQ+kD, tT(:, i));
    [tIi, intermediateHasAddedT] = checkTi(kQ+kD, tI(:, i));
    
    initialValues = [znt zntFe3 zni zniFe2 zng zngFe3]; %zni and zniFe2 are 0.
    [~, yT] = ode45(@fun, tTi, initialValues);
    [~, yI] = ode45(@fun, tIi, initialValues);
    
    triplet       = aTrip .* (yT(:, 1) + yT(:, 2)) + c;
    intermediate  = aInt .* (yI(:, 3) + yI(:,4)) + aJump .* (yI(:,1) + yI(:,2));
    
    
    if tripletHasAddedT
        triplet = triplet(2:end); %Remove added ti
    end %if
    if intermediateHasAddedT  
        intermediate = intermediate(2:end); %Remove added ti
    end %if
    
    Y(:,i) = [triplet; intermediate];
end %for    

    function dydt = fun(~, y)
        fe3 = fe3Added - y(2) - y(3) - y(4) - y(6);
        fe2 = fe2Added + y(3); 
        dydt = [-(kD)*y(1);     ... %znt
                -(kQ + kD)*y(2); ... %zntFe3
                - kOn*fe2*y(3) + kOff*y(4);     ... %zni
                -(kB + kOff)*y(4) + kOn*fe2*y(3) + kQ*y(2); ... %zniFe2
                + kD*y(1); ...  %zng
                + (kD+kQ)*y(2) + kB*y(4)];    %zngFe3
    end %fun
end %intermediatedissociationfun
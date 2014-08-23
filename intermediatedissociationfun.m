function Y = intermediatedissociationfun(param, t, constants, concentrations)
%Y = intermediatedissociationfun(param, t, constants, concentrations)
%
%Simulates intermediate data for a simple dissociation model including
%dissociation of the intermediate, but with strictly slow exchange triplet decay.
%
%
%param          = [aInt aJump kD kQ kB log(kOn) kOff]; Vector or array where
%                 each row corresponds to one trace.
%constants      = [log(Ka) phi]; Vector or array where each row corresponds to
%                 one trace. Ka is the association constant in M-1. phi is the
%                 quantum yield of triplet. This is necessary for determining
%                 initial concentrations of the excited state species for second
%                 order reactions.
%concentrations = [zntAdded(uM) fe3Added(uM) fe2Added(uM)]; Vector or array
%                 where each row corresponds to one trace. Concentrations must
%                 be in micromolar.
%t              - Vertical vector or array where each column corresponds to one 
%                 trace.
%
%Inputs must all be vectors or all be arrays.

paramSize          = size(param);
tSize              = size(t);
constantsSize      = size(constants);
concentrationsSize = size(concentrations);

nTraces = tSize(2);

%Error if mismatched arrays
validInputs = paramSize(1) == nTraces & constantsSize == nTraces & ...
    concentrationsSize == nTraces;

if ~validInputs
    error('Param rows, t cols, constants rows, and concentrations row must be equal')
end %if


Y = zeros(tSize);

for i = 1:nTraces
    
    %Parse parameters
    %----------------
    [aInt, aJump, kD, kQ, kB, kOn, kOff] = ...
    deal(param(i, 1), param(i, 2), param(i, 3), param(i, 4), param(i, 5), ...
         param(i, 6), param(i, 7));


    zntAdded = concentrations(i, 1);
    fe3Added = concentrations(i, 2);
    fe2Added = concentrations(i, 3);
    Ka       = constants(i, 1);
    phi      = constants(i, 2);
    
    kOn = 10^kOn * 10^-6;
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

    %Check if t(1) is good estimate for t(0), add ti if needed.
    %----------------
    [iT, hasAddedT] = checkTi(kQ+kD, t(:, i));
    
    
    initialValues = [znt zntFe3 zni zniFe2 zng zngFe3]; %zni and zniFe2 are 0.
    [~, yi] = ode45(@fun, iT, initialValues);
    y       = aInt .* (yi(:, 3) + yi(:,4)) + aJump .* (yi(:,1) + yi(:,2));
    
    
    if hasAddedT
        y = y(2:end); %Remove added ti
    end %if
    
    Y(:, i) = y;
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
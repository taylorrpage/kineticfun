function [d, a, da, daa] = bindingtwosite(K1, K2, d0, a0)
% ==============================================================================
% Calculates the equilibrium concentrations of species using a sequential
% binding model:
%
% [D, A, DA, DAA] = bindingtwosite(K1, K2, D0, A0) - Returns the equilibrium
% concentrations of D, A, DA, and DAA.
%
%   D + A <-> DA    K1
%   DA + A <-> DAA  K2
%
%   NOTES:
%       K1-     Binding constant for the first association
%       K2-     Binding constant for the second association
%       D0-     Initial donor concentration in micromolar
%       A0-     Initial acceptor concentration in micromolar
% ==============================================================================

%Convert binding constants to dissociation constants.
Kd1 = 1/K1 * 10^6;
Kd2 = 1/K2 * 10^6;

%Polynomial coefficients
c1 = 4/(Kd2^2) - 1/(Kd1*Kd2);
c2 = 4/Kd2 + 2*d0/(Kd1*Kd2) - 1/Kd1;
c3 = 1 - a0.*(2*d0/Kd2 - 1 - a0./Kd2)./Kd1 + d0/Kd1;
c4 = d0.*a0./Kd1;

%Cubic solution equations
a  = (27*c1^2.*c3 - 9*c1*c2^2)/(3*c1)^3;
b  = (3*c2^3 - 9*c1*c2.*c3 - 27*c1^2.*c4 - c2^3)/(3*c1)^3;
t  = 3.*b./(2.*a.*sqrt(-a./3));
r3 = (acosd(t))./3 + 4*180/3;    %This solution requires degrees
x3 = 2.*sqrt(-a/3).*cosd(r3) - c2/(3*c1);   %This solution requires degrees

%Calculating the final species
daa = (a0.*x3 - x3.^2)./(Kd2 + 2.*x3);
w1  = (x3 + sqrt(x3.^2 - 4.*daa.*(d0 - x3 - daa)))./2;
w2  = (x3 - sqrt(x3.^2 - 4.*daa.*(d0 - x3 - daa)))./2;
da  = w1 + w2;
d   = d0 - da - daa;
a   = a0 - da - 2.*daa;



function mu = phi2mu(phi);
% PHI2MU    Convert phi units to microns.
%
%   Syntax: mu = phi2mu(phi)
%
%   John Newgard, April 1, 2010

mu = 2.^(-phi).*1000;

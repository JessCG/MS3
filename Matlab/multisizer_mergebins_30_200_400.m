% MULTISIZER_MERGEBINS_30_200_400   Get the diameter bin centres.
%
%   D = MULTISIZER_MERGEBINS_30_200_400 returns the diameters (in microns)
%       corresponding to the central size bins of the merged Multisizer III 
%       distributions for the 30-, 200-, and 400-micron tubes.
%
%   [D,D_L,D_U] = ... returns the central bin diameters as well as the
%       values corresponding to the lower and upper bin edges.
%
%   See also MU2PHI, PHI2MU
%
%   By John Newgard, 2 April 2010

function [centd,varargout] = multisizer_mergebins_30_200_400

delta_phi = -.2; % bin-centre interval, in phi units
a = mu2phi(1); % diameter, in phi units, corresponding to 1 micron
b = 1.5864059; % this is the lowest diameter, in phi units, of the 3 distributions
centd_phi = [a:delta_phi:b]';
centd_mu = phi2mu(centd_phi);

% Prepare output
centd = centd_mu;
if nargout==3
    lowerd_phi = centd_phi-0.5*delta_phi;
    upperd_phi = centd_phi+0.5*delta_phi;
    lowerd_mu = phi2mu(lowerd_phi);
    upperd_mu = phi2mu(upperd_phi);
    
    varargout = {lowerd_mu,upperd_mu};
elseif nargout >1
    disp('ERROR: You must use either 1 or 3 output arguments.')
    return
end


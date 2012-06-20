function phi=mu2phi(mu);
% 
% This function converts diameters micrometres to phi.  It requires a
% value in microns as the input and the output is diameter in phi units.

phi = -log2(mu/1000);

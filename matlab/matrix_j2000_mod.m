function P = matrix_j2000_mod(JD)
% MATRIX_J2000_MOD - Compute the rotation matrix between the J2000 and MoD
% frames.
%
% The rotation matrix between the J2000 and MOD frames is the Precession
% Matrix.
%
% INPUTS:
%   JD         The Julian Date.
%
% OUTPUTS:
%   P          The Precession Matrix.
%
% References:
% [1] E. Suirana, J. Zoronoza, M. Hernandez-Pajares - GNSS Data Processing -
% Volume I: Fundamentals and Algorithms, ESA 2013.

% Number of Julian centuries from J2000 of TDT (A.19) in [1]: 
T = (JD - 2451545.0)/36525.0;

% IAU 1976 Precession Model (A.23):
coeffs = [0.6406161388,  3.0407777777e-04,  5.0563888888e-06; ...
          0.5567530277, -1.1851388888e-04, -1.1620277777e-05; ...
          0.6406161388,  8.3855555555e-05,  4.9994444444e-06] * [T; T*T; T*T*T];

z = coeffs(1);
theta = coeffs(2);
zeta = coeffs(3);

% (A.22)
P = matrix_rot3d(-z) * matrix_rot2d(theta) * matrix_rot3d(-zeta);

end
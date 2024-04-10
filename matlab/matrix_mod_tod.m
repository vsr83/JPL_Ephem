function N = matrix_mod_tod(JT, eps, deps, dpsi)
% MATRIX_MOD_TOD - Compute the rotation matrix between the MoD and ToD
% frames.
%
% The rotation matrix between the MoD and ToD frames is the Nutation
% Matrix.
%
% INPUTS:
%   JD         The Julian Time.
%   eps        Nutation parameter eps (optional, if parameter computation
%                                      is skipped)
%   deps       Nutation parameter deps (optional, if parameter computation
%                                       is skipped)
%   dpsi       Nutation parameter dpsi (optional, if parameter computation
%                                       is skipped)
%
% OUTPUTS:
%   N          The Nutation Matrix.
%
% References:
% [1] E. Suirana, J. Zoronoza, M. Hernandez-Pajares - GNSS Data Processing -
% Volume I: Fundamentals and Algorithms, ESA 2013.

if nargin < 4 
    [eps, deps, dpsi] = nutation_iau1980(JT);
end

% (A.24)
N = matrix_rot1d(-eps - deps) * matrix_rot3d(-dpsi) * matrix_rot1d(eps);

end
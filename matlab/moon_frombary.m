function [R, V] = moon_frombary(R, V, mu, ind_earth, ind_moon)
% MOON_FROMBARY - Convert the coordinates of the Moon and the Earth from 
% solar system barycentric coordinates.
%
% For numerical reasons, the Moon position is integrated relative to the
% Earth rather than the solar system barycenter. Thus, the Earth and Moon
% coordinates are transformed to the Earth-Moon barycenter and the position 
% vector of the Moon relative to the Earth-Moon barycenter.
%
% This implements equations (6)-(7) in [1] and (8.16)-(8.17) in [2].
%
% INPUTS:
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)
%   mu         The standard gravitational parameter for each point-mass
%              (au^3/d^2, num_targets)
%   ind_earth  The row corresponding to Earth in the matrices R and V.
%   ind_moon   The row corresponding to Moon in the matrices R and V.
%
% OUTPUTS:
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.

mu_earth = mu(ind_earth);
mu_moon = mu(ind_moon);

% We need temporary variables since the vectors are used more than once.
r_B = (mu_earth * R(ind_earth, :) + mu_moon * R(ind_moon), :) ...
    / (mu_earth + mu_moon);
r_M = R(ind_moon, :) - R(ind_earth, :);
v_B = (mu_earth * V(ind_earth, :) + mu_moon * V(ind_moon, :)) ...
    / (mu_earth + mu_moon);
v_M = V(ind_moon, :) - V(ind_earth, :);

R(ind_earth, :) = r_B;
R(ind_moon, :) = r_M;
V(ind_earth, :) = v_B;
V(ind_moon, :) = v_M;

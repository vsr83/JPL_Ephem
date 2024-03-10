function [R, V] = adjust_sun(R, V, mu, ind_sun)
% ADJUST_SUN - Compute the position and velocity of the Sun so that the
% relativistic barycenter of the Solar System is at origin.
%
% The position and velocity of Sun is not integrated in the actual 
% numerical solution. Instead, the position is obtained from the equations
% for the relativistic barycenter so that the the barycenter of the Solar
% System is always at the origin [1], [2].
%
% This implements equations (4)-(5) in [1] and (8.2) in [2]. The method
% corresponds to the method fixsun in [3].
%
% INPUTS:
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)
%   mu         The standard gravitational parameter for each point-mass
%              (au^3/d^2, num_targets)
%   ind_sun    The row corresponding to Sun in the matrices R and V.
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
%  [3] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip
%

% In [1], [3], the number of iterations is 2. [2] does not specify any
% number.
num_iterations = 2;

[num_targets, tmp] = size(R);

% Speed of light (au/d).
c = 173.144632720536344565;
c2 = c * c;

% Distance between each pair of point-masses.
Rij = zeros(num_targets, num_targets);

% Relativistic standard gravitational parameter for the evaluation of the
% barycenter.
mu_star = zeros(num_targets, 1);

for ind_iteration = 1:num_iterations
    % Compute the relativistic standard gravitational parameters:
    for ind_target = 1:num_targets
        tmp = norm(V(ind_target, :))^2;

        for ind_source = 1:num_targets
            if ind_source == ind_target
                continue;
            end

            tmp = tmp - mu(ind_source) / norm(R(ind_target, :) - R(ind_source, :));
        end
        tmp = 1 - tmp * 0.5 / c2;
        mu_star(ind_target) = mu(ind_target) * tmp;
    end

    % Compute the part of the equation for barycenter not including the
    % contribution from the Sun.
    r = [0, 0, 0];
    v = [0, 0, 0];
    for ind_target = 1:num_targets
        if ind_target == ind_sun
            continue;
        end
        r = r + mu_star(ind_target) * R(ind_target, :);
        v = v + mu_star(ind_target) * V(ind_target, :);
    end

    % Position and velocity of the Sun.
    R(ind_sun, :) = -(1 / mu_star(ind_sun)) * r;
    V(ind_sun, :) = -(1 / mu_star(ind_sun)) * v;
end
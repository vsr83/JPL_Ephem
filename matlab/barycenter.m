function [r_bary, v_bary] = barycenter(R, V, mu, with_rel)
% BARYCENTER - Compute the classical or relativistic barycenter for the 
% given point masses.
%
% INPUTS:
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)
%   mu         The standard gravitational parameter for each point-mass
%              (au^3/d^2, num_targets)
%   with_rel   Flag indicating whether to include relativistic effects.
%
% OUTPUTS:
%   r_bary     The barycenter position (au, 1x3)
%   v_bary     The barycenter velocity (au/d, 1x3)
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.

[num_targets, tmp] = size(R);

% Relativistic or classical standard gravitational parameter for the 
% evaluation of the barycenter. 
mu_star = zeros(num_targets, 1);

if with_rel
    % Speed of light (au/d).
    c = 173.144632720536344565;
    c2 = c * c;
    
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
else 
    mu_star = mu;
end

% Compute the part of the equation for barycenter not including the
% contribution from the Sun.
r_bary = [0, 0, 0];
v_bary = [0, 0, 0];
mu_sum = sum(mu_star);

for ind_target = 1:num_targets
    r_bary = r_bary + mu_star(ind_target) * R(ind_target, :);
    v_bary = v_bary + mu_star(ind_target) * V(ind_target, :);
end

r_bary = r_bary / mu_sum;
v_bary = v_bary / mu_sum;
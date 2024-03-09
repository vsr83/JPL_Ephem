function [A_rel, A_newt] = acc_pointmass(R, V, mu)
% ACC_POINTMASS - Compute the relativistic and Newtonian parts of 
% acceleration.
%
% This routine computes the relativistic and Newtonian parts of the
% acceleration for arbitrary number of point-masses. The total
% accelerations are obtained as the sum of the two matrices. The
% computation is based on the equation (1) in [1] and (8.1) in [2].
%
% The method has been extensively tested by comparisons to DE118/DE200 
% ephemeris implemented by Steve Moshier [3].
%
% Important: This method has not been implemented for performance. Rather,
% it is supposed to serve as a basis of implementation.
%
% INPUTS:
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)
%   mu         The standard gravitational parameter for each point-mass
%              (au^3/d^2, num_targets)
%
% OUTPUTS:
%   A_rel      The relativistic acceleration   (au/d^2, num_targets x 3)
%   A_newt     The Newtonian acceleration      (au/d^2, num_targets x 3)
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning fort-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [3] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip
%

[num_targets, tmp] = size(R);

A_rel = zeros(num_targets, 3);
A_newt = zeros(num_targets, 3);

% Distance between each pair of point-masses.
Rij = zeros(num_targets, num_targets);

for ind_target = 1:num_targets
    for ind_source = 1:num_targets
        Rij(ind_source, ind_target) = norm(R(ind_target, :) - R(ind_source, :));
    end
end

% Third power of the distance between each pair of point-masses.
Rij3 = Rij.^3;

% Speed of light (au/d).
c = 173.144632720536344565;
c2 = c * c;

% For numerical accuracy, it is very important to compute the relativistic 
% accceleration separately. Otherwise, one has to do large amount of 
% floating computations that involve adding small numbers to larger ones.

% Compute the Newtonian accelerations first.
for ind_target = 1:num_targets
    for ind_source = 1:num_targets
        if ind_source == ind_target
            continue;
        end
        r_target_source = R(ind_source, :) - R(ind_target, :);
        acc_newton = (mu(ind_source) / Rij3(ind_source, ind_target)) ...
                   * r_target_source;

        A_newt(ind_target, :) = A_newt(ind_target, :) + acc_newton;
    end
end

% Compute the relativistic accelerations.
for ind_target = 1:num_targets
    acc_target = zeros(1, 3);
    r_target = R(ind_target, :);
    v_target = V(ind_target, :);

    for ind_source = 1:num_targets
        if ind_source == ind_target
            continue;
        end
        r_source = R(ind_source, :);
        v_source = V(ind_source, :);

        r_target_source = r_source - r_target;
        v_target_source = v_source - v_target;

        % Newtonian part of the acceleration.
        acc_newton = (mu(ind_source) / Rij3(ind_source, ind_target)) ...
                   * r_target_source;

        % The first part of the acceleration formula involves
        % multiplication of the Newtonian acceleration with an expression.
        newton_mult = 0;
        for ind_target2 = 1:num_targets
            if ind_target ~= ind_target2
                newton_mult = newton_mult ...
                            - (4 / c2) * mu(ind_target2) / Rij(ind_target, ind_target2);
            end

            if ind_source ~= ind_target2
                newton_mult = newton_mult ...
                            - (1 / c2) * mu(ind_target2) / Rij(ind_source, ind_target2);
            end
        end

        newton_mult = newton_mult ...
                    +     norm(V(ind_target, :))^2 / c2 ...
                    + 2 * norm(V(ind_source, :))^2 / c2 ...
                    - (4/c2) * dot(V(ind_target, :), V(ind_source, :)) ...
                    - (1.5/c2) * (dot(r_target_source, V(ind_source, :)) / Rij(ind_source, ind_target))^2 ...
                    + (0.5/c2) * dot(r_target_source, A_newt(ind_source, :));

        acc_target = acc_target ...
                   + newton_mult * acc_newton ...
                   + (1/c2) * mu(ind_source) / Rij3(ind_source, ind_target) ...
                   * dot(r_target_source, 4 * v_target - 3 * v_source) * v_target_source ...
                   + (3.5/c2) * mu(ind_source) * A_newt(ind_source, :) ...
                   / Rij(ind_source, ind_target);
    end
    
    A_rel(ind_target, :) = acc_target;
end

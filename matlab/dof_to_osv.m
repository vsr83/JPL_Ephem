function [R, V] = dof_to_osv(y);
% DOF_TO_OSV - Transform degrees of freedom to OSVs for n point masses.
%   
% INPUTS: 
%   y          Degrees of freedom (num_targets * 6 x 1).
% OUTPUTS:
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)

[num_equations, tmp] = size(y);
num_targets = num_equations / 6;

% Convert degrees of freedom to OSVs.
R = zeros(num_targets, 3);
V = zeros(num_targets, 3);

for ind_target = 1:num_targets
    ind_start = (ind_target - 1) * 6;
    R(ind_target, :) = y(ind_start + (1:3));
    V(ind_target, :) = y(ind_start + (4:6));
end

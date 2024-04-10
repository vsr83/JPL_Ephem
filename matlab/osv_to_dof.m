function y = osv_to_dof(R, V);
% OSV_TO_DOF - Transform OSVs into degrees of freedom for n point masses.
%   
% INPUTS: 
%   R          The position of each point-mass (au, num_targets x 3)
%   V          The velocity of each point-mass (au/d, num_targets x 3)
% OUTPUTS:
%   y          Degrees of freedom (num_targets * 6 x 1).

[num_targets, tmp] = size(R);
num_equations = num_targets * 6;

y = zeros(num_equations, 1);

for ind_target = 1:num_targets
    ind_start = (ind_target - 1) * 6;
    y(ind_start + (1:3)) = R(ind_target, :)';
    y(ind_start + (4:6)) = V(ind_target, :)';
end

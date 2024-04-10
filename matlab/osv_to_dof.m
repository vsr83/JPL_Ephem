function y = osv_to_dof(R, V);

[num_targets, tmp] = size(R);
num_equations = num_targets * 6;

y = zeros(num_equations, 1);

for ind_target = 1:num_targets
    ind_start = (ind_target - 1) * 6;
    y(ind_start + (1:3)) = R(ind_target, :)';
    y(ind_start + (4:6)) = V(ind_target, :)';
end

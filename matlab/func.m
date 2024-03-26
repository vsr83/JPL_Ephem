function k = func(t, y, mu)

[num_equations, tmp] = size(y);
num_targets = num_equations / 6;
k = zeros(num_equations, 1);

[R, V] = dof_to_osv(y);
[A_rel, A_newt] = acc_pointmass(R, V, mu);
A = A_newt + A_rel;

for ind_target = 1:num_targets
    ind_start = (ind_target - 1) * 6;
    k(ind_start + (1:3)) = V(ind_target, :)';
    k(ind_start + (4:6)) = A(ind_target, :)';
end

function k = func_lib(t, y, mu)
% FUNC_LIB - Compute the time derivative of the degrees of freedom for 
% numerical integration including librations.
%
% INPUTS:
%   t          The time after epoch at evaluation.
%   y          The degrees of freedom at evaluation 
%              (6 * (num_targets + 1) x 1).
%   mu         The standard gravitational parameter for each point-mass
%              (au^3/d^2, num_targets)
%
% OUTPUTS:
%   k          The time derivative of the DoFs for numerical integration.
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [3] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip

% Ugly: This should be replaced, if developed further.
JT = t + 2440400.50;

% Determine the number of targets and equations.
[num_equations, tmp] = size(y);
num_targets = num_equations / 6 - 1;
k = zeros(num_equations, 1);

index_sun   = 1;
index_earth = 4;
index_moon  = 5;

% Ugly: This should be replaced, if developed further.
mu_sun   = 2.959122082855911e-04;
mu_earth = 8.887692445123496e-10;
mu_moon  = 1.093189450705846e-11;
MU = [mu_sun, mu_earth, mu_moon];

% The coefficients are apart from S_22 and C_22 are taken from Tables 8.4 
% and 8.5 in [2]. The values in a proper implementation should function 
% only as "rigid" inputs for more detailed computation of the harmonic 
% coefficients. The coefficients S_22 and C_22 are from [3].

% Zonal harmonics J_2, J_3 and J_4 in Earth's potential.
Je = [0.001082626, -0.000002533, -0.000001616];

% Zonal harmonics J_2, J_3 and J_4 in Moon's potential.
Jm = [2.04312007e-4; 0.000008785470; -0.000000145383];

% Tesseral harmonics in Moon's potential.
CSnm = [
    2, 2,  2.230351309e-5,    0.0;; ...
    3, 1,  3.07082741328e-5,  5.61066891941e-6; ...
    3, 2,  4.88840471683e-6,  1.68743052295e-6; ...
    3, 3,  1.43603108489e-6, -3.343544677e-7; ...
    4, 1, -7.17780149806e-6,  2.94743374914e-6; ...
    4, 2, -1.43951838385e-6, -2.8843721272e-6; ...
    4, 3, -8.54788154819e-8, -7.88967312839e-7; ...
    4, 4, -1.5490389313e-7,   5.6404155572e-8
];

CSnm = [
    2, 2,  2.230351309e-5,    0.0;; ...
    3, 1,  0.000030803810,    0.000004259329; ...
    3, 2,  0.000004879807,    0.000001695516; ...
    3, 3,  0.000001770176,   -0.000000270970; ...
    4, 1, -0.000007177801,    0.000002947434; ...
    4, 2, -0.000001439518,   -0.000002884372; ...
    4, 3, -0.000000085479,   -0.000000788967; ...
    4, 4, -0.000000154904,    0.000000056404
];

[R, V] = dof_to_osv(y(7:end));
[A_rel, A_newt] = acc_pointmass(R, V, mu);
A = A_newt + A_rel;

OSV = [
    R(index_sun, :)', R(index_earth, :)', R(index_moon, :)'; ...
    V(index_sun, :)', V(index_earth, :)', V(index_moon, :)'
];

LIB = y(1:6);

[A_oblate, acc_moon] = acc_oblateness(OSV, MU, LIB, JT, Je, Jm, CSnm, false);
 
A(index_sun,   :) = A(index_sun,   :) + A_oblate(:, 1)';
A(index_earth, :) = A(index_earth, :) + A_oblate(:, 2)';
A(index_moon,  :) = A(index_moon,  :) + A_oblate(:, 3)';

k(1) = y(2);
k(2) = acc_moon(1);
k(3) = y(4);
k(4) = acc_moon(2);
k(5) = y(6);
k(6) = acc_moon(3);

for ind_target = 1:num_targets
    ind_start = ind_target * 6;
    k(ind_start + (1:3)) = V(ind_target, :)';
    k(ind_start + (4:6)) = A(ind_target, :)';
end

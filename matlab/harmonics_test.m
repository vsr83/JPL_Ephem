% This scripts tests the computation of the point mass in body frame with
% inputs and outputs from the DE118 implementation by Steve Moshier.

au = 149597870.691;

% Zonal harmonics J_2, J_3 and J_4 in Earth's potential.
Je = [1.08263e-3; -2.54e-6; -1.61e-6];

% Zonal harmonics J_2, J_3 and J_4 in Moon's potential.
Jm = [2.02150907893e-4; 1.21260448837e-5; -1.45383007072e-7];

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

r_body_earth_sun = [
    -0.10851582874626557185; ...
     0.92735373714107316445; ...
     0.40210907185031696809 ...
 ];

% Equatorial radius of the Earth (km).
a_earth = 6.378139999385382e+03;

mu = 1;
a_1 = acc_body(r_body_earth_sun, a_earth/au, -mu, Je, []);
a_1_exp = [
    0.0000000000000642325883410871061839355210; ...
    -0.0000000000005489183608746217128182116910; ...
    -0.0000000000024245505027412273649070483729];
rel_error_1 = norm(a_1_exp - a_1) / norm(a_1_exp)

r_body_moon_sun = [
    0.94374576705943424848; ...
    0.38268794216078455550;
    0.02710296020739921818    
];

% Equatorial radius of the Moon (km).
a_moon = 1.737999999832521e+03;
mu = 1;
a_2 = acc_body(r_body_moon_sun, a_moon/au, -mu, Jm, CSnm);
a_2_exp = [
    -0.0000000000000473892866153501820441870181; ...
    -0.0000000000000318149662536360739254924260; ...
    -0.0000000000000038287252192564183969255991 ...
];

rel_error_2 = norm(a_2_exp - a_2) / norm(a_2_exp)

r_body_moon_earth = [
    0.00239026490334162461; ...
    -0.00017227335425532372; ...
    0.00026607374895177280 ...
];
a_3 = acc_body(r_body_moon_earth, a_moon/au, -mu, Jm, CSnm);
a_3_exp = [
    -0.00189164806069694151; ...
    0.00021261987821214140; ...
    -0.00053425131651080283
];
rel_error_3 = norm(a_3_exp - a_3) / norm(a_3_exp)

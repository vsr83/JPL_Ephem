function [A, acc_moon] = acc_oblateness(OSV, MU, lib_moon, JT, Je, Jm, CSnm, use_de118)
% ACC_OBLATENESS - Compute the accelerations due to Earth and Moon figure 
%                  and tides.
% 
% This method is heavily based on the oblate method in [2].
%
% Important: This method has not been implemented for performance but for 
% simplicity and helping with an implementation. The large number of
% unnecessary intermediate variables is for testing.
%
% INPUTS:
%   OSV        6x3 matrix [osv_sun, osv_earth, osv_moon], where
%     osv_sun    The orbit state vector for the Sun.
%     osv_earth  The orbit state vector for the Earth.
%     osv_moon   The orbit state vector for the Moon.
%   MU         1x3 vector [mu_sun, mu_earth, mu_moon] with standard 
%              gravitational parameters (au^3/d^2).
%   lib_moon   The libration state of the moon [phi; phi1; theta; theta1; 
%              psi; psi1] (rad, rad/d).
%   JT         Julian time.
%   Je         Zonal harmonics for the extended body of Earth starting 
%              from n = 2 (num_harmonics_earth x 1).
%   Jm         Zonal harmonics for the extended body of Moon starting from 
%              n = 2 (num_harmonics x 1).
%   CSnm       Tesseral harmonics in the (n, m, C_nm, Snm) row format 
%              (num_harmonics x 4).
%   use_de118  Flag to indicate whether to use DE118 frame instead of J2000
%              for the computations.
%
%  The orbit state vectors have the format [x; y; z; v_x; v_y; v_z] in the
%  DE118/J2000 frame with barycentric origin and units (au, au/d).
%
% OUTPUTS:
%   A          6x3 matrix [acc_sun, acc_earth, acc_moon], where
%     acc_sun    The acceleration of the Sun due to figure effects.
%     acc_earth  The acceleration of the Earth due to figure effects.
%     acc_moon   The acceleration of the Moon due to figure effects.
%   acc_moon   The angular accelerations [phi2; theta2; psi2] of the Moon
%              (rad/d^2).
%  
%  The accelerations are 3x1 vectors (au/d^2).
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip

mu_s = MU(1);
mu_e = MU(2);
mu_m = MU(3);

% Parse OSVs.
r_s = OSV(1:3, 1);
r_e = OSV(1:3, 2);
r_m = OSV(1:3, 3);
v_s = OSV(4:6, 1);
v_e = OSV(4:6, 2);
v_m = OSV(4:6, 3);

% Parse libration angles.
phi    = lib_moon(1);
phi1   = lib_moon(2);
theta  = lib_moon(3);
theta1 = lib_moon(4);
psi    = lib_moon(5);
psi1   = lib_moon(6);

% Matrix from DE118/J2000 coordinates to body coordinates.
matrix_body = matrix_to_body(phi, theta, psi);

% Precession matrix for the transformation from J2000 to MoD coordinates. 
% DE118 TBD.
matrix_precession = matrix_j2000_mod(JT);

% Nutation matrix for the transformation from MoD to ToD coordinates.
matrix_nutation = matrix_mod_tod(JT);

% Astronomical unit (km).
au = 149597870.691;

% Equatorial radius of the Moon (au).
a_moon = 1.737999999832521e+03 / au;

% Equatorial radius of the Earth (au).
a_earth = 6378.14 / au;

% The position of the Earth w.r.t. Moon body center in DE118/J2000 and 
% body coordinates.
r_em_j2000 = r_e - r_m;
r_em_body = matrix_body * r_em_j2000;
% Acceleration/mu and
[acc_em_body_tmp, T_earth] = acc_body(r_em_body, a_moon, 1, Jm, CSnm);

% 1. Accelerations from the interaction between the Moon figure and Earth.
acc_em_body      = -mu_m * acc_em_body_tmp;
acc_em_j2000_fig = matrix_body' * acc_em_body;
acc_me_body      = mu_e * acc_em_body_tmp;
acc_me_j2000_fig = matrix_body' * acc_me_body;

r_sm_j2000 = r_s - r_m;
r_sm_body  = matrix_body * r_sm_j2000;
[acc_sm_body_tmp, T_sun] = acc_body(r_sm_body, a_moon, 1, Jm, CSnm);

% 2. Accelerations from the interaction between the Moon figure and Sun.
acc_sm_body      = -mu_m * acc_sm_body_tmp;
acc_sm_j2000_fig = matrix_body' * acc_sm_body
acc_ms_body      = mu_s * acc_sm_body_tmp;
acc_ms_j2000_fig = matrix_body' * acc_ms_body

% 3. Libration of the Moon.

% Compute the total torque on the Moon and the angular accelerations.
T = mu_e * T_earth + mu_s * T_sun;

[phi2, theta2, psi2] = libration_moon(phi, theta, psi, phi1, theta1, psi1, T);
acc_moon = [phi2; theta2; psi2];

% 4. Oblateness of the Earth.

% The position of the Moon w.r.t. Earth body center in DE118/J2000.
%r_me_j2000 = r_m - r_e;
r_me_j2000 = matrix_de118_j2000 * (r_m - r_e);

% The position of the Sun w.r.t. Earth body center in DE118/J2000.
%r_se_j2000 = r_s - r_e;
r_se_j2000 = matrix_de118_j2000 * (r_s - r_e);

% Transform the relative position of the Moon to the True-of-Date frame.
r_me_mod = matrix_precession * r_me_j2000;
r_me_tod = matrix_nutation * r_me_mod;

% Transform the relative position of the Sun to the True-of-Date frame.
r_se_mod = matrix_precession * r_se_j2000;
r_se_tod = matrix_nutation * r_se_mod;

[acc_me_tod_tmp, T_sun] = acc_body(r_me_tod, a_earth, 1, Je, []);
[acc_se_tod_tmp, T_sun] = acc_body(r_se_tod, a_earth, 1, Je, []);

acc_me_tod = -mu_e * acc_me_tod_tmp;
acc_se_tod = -mu_e * acc_se_tod_tmp;
acc_em_tod = mu_m * acc_me_tod_tmp;
acc_es_tod = mu_s * acc_se_tod_tmp;

% 5. Accelerations from the interaction between Earth tides and the Moon.
[acc_me_tod_tides, acc_em_tod_tides] = acc_tides(r_me_tod, mu_e, mu_m);

% Convert accelerations from Earth oblateness and tides to J2000 frame.
%matrix_tod_j2000 = matrix_precession' * matrix_nutation';
matrix_tod_j2000 = matrix_de118_j2000' * matrix_precession' * matrix_nutation';
acc_me_j2000_obl = matrix_tod_j2000 * acc_me_tod;
acc_se_j2000_obl = matrix_tod_j2000 * acc_se_tod;
acc_em_j2000_obl = matrix_tod_j2000 * acc_em_tod;
acc_es_j2000_obl = matrix_tod_j2000 * acc_es_tod;
acc_me_j2000_tides = matrix_tod_j2000 * acc_me_tod_tides;
acc_em_j2000_tides = matrix_tod_j2000 * acc_em_tod_tides;

acc_s_j2000 = acc_sm_j2000_fig + acc_se_j2000_obl ;
acc_e_j2000 = acc_es_j2000_obl + acc_em_j2000_fig + acc_em_j2000_obl + acc_em_j2000_tides;
acc_m_j2000 = acc_ms_j2000_fig + acc_me_j2000_fig + acc_me_j2000_obl + acc_me_j2000_tides;

% Assemble outputs to an acceleration matrix.
A = [acc_s_j2000, acc_e_j2000, acc_m_j2000];


disp('Moon Figure <-> Earth : Earth Acceleration');
disp(acc_em_j2000_fig');
disp('Moon Figure <-> Earth : Moon Acceleration');
disp(acc_me_j2000_fig');
disp('Moon Figure <-> Sun : Sun Acceleration');
disp(acc_sm_j2000_fig');
disp('Moon Figure <-> Sun : Moon Acceleration');
disp(acc_ms_j2000_fig');
disp('Earth Oblateness <-> Moon : Earth Acceleration');
disp(acc_em_j2000_obl');
disp('Earth Oblateness <-> Moon : Moon Acceleration');
disp(acc_me_j2000_obl');
disp('Earth Oblateness <-> Sun : Earth Acceleration');
disp(acc_es_j2000_obl');
disp('Earth Oblateness <-> Sun : Sun Acceleration');
disp(acc_se_j2000_obl');
disp('Earth Tides <-> Moon : Earth Acceleration');
disp(acc_em_j2000_tides');
disp('Earth Tides <-> Moon : Moon Acceleration');
disp(acc_me_j2000_tides');

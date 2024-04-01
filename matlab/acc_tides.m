function acc_point = acc_tides(r_me_tod, mu_e, mu_m)
% ACC_TIDES - Compute the acceleration 
%
% This method computes the expression in equation (2) of [1] or (8.3) in
% [2] and transforms the acceleration to body coordinates.
%
% Important: This method has not been implemented for performance but for 
% simplicity and helping with an implementation. 
%
% INPUTS:
%   r_me_tod   The position of the Moon w.r.t. Earth in the ToD frame 
%   (au, 3 x 1)
%   mu_e       Standard gravitational parameter (au^3/d^2) for Earth.
%   mu_m       Standard gravitational parameter (au^3/d^2) for Moon.
%
% OUTPUTS:
%   acc_point  The acceleration of the point mass in body coordinates 
%              (au/d^2, num_targets x 3).
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip

% Potential love number for the Earth [2]. In [1], the value is 0.29.
love = 0.30;

% Astronomical unit in km.
au = 1.495978706910000e+08;

% Equatorial radius of Earth [2]. 6378.156 km in [1].
a = 6378.14 / au;

% Distance between Earth and the Moon.
r_em = norm(r_me_tod);

% Angle between the bulge and the Earth-Moon line. 0.04635 rad in [1].
phase = 4.0700012e-2;

acc_point = -(3 * love * mu_m) * (1 + mu_m/mu_e) * (a^5 / r_em^8) * ...
    [r_me_tod(1) + phase * r_me_tod(2) ; ...
     r_me_tod(2) - phase * r_me_tod(1) ; ...
     r_me_tod(3)];

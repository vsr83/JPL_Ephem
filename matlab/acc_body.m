function [acc_point, T_body] = acc_body(r_point, a, mu, Jn, CSnm)
% ACC_BODY - Compute the acceleration and torque due to zonal and tesseral 
% harmonics from an extended body.
%
% This method computes the expression in equation (2) of [1] or (8.3) in
% [2] and transforms the acceleration to body coordinates.
%
% Important: This method has not been implemented for performance but for 
% simplicity and helping with an implementation. 
%
% INPUTS:
%   r_point    The position of the point-mass w.r.t. the body center 
%              in body coordinates (au, 3 x 1)
%   mu         Standard gravitational parameter (au^3/d^2) or 1 if
%              the results are multiplied with -mu afterwards.
%   a          The equatorial radius of the extended body (au).
%   Jn         Zonal harmonics for the extended body starting from n = 2 
%              (num_harmonics x 1).
%   CSnm       Tesseral harmonics in the (n, m, C_nm, Snm) row format 
%              (num_harmonics x 4).
%
% OUTPUTS:
%   acc_point  The acceleration of the point mass in body coordinates 
%              (au/d^2, num_targets x 3).
%   T_body     Torque per unit mass on the body in body coordinates 
%              (au^2/d^2, 3 x 1). If mu is set to 1
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [3] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip

% Distance between center of the extended body and the point mass.
r = norm(r_point);

% Latitude and longitude of the point mass w.r.t. body coordinates (rad).
sin_lat = r_point(3) / norm(r_point);
lat_point = asin(sin_lat);
lon_point = atan2(r_point(2), r_point(1));
cos_lat = cos(lat_point);

% Number of zonal harmonics starting from n=2.
number_zonal = length(Jn);
number_tesseral = length(CSnm);

acc_point_zonal = [0, 0, 0]';
acc_point_tesseral = [0, 0, 0]';
acc_point = [0, 0, 0]';

% Evaluate zonal harmonics.
for ind_zonal = 1:number_zonal
    n = ind_zonal + 1;

    % Legendre value and derivative terms.
    Pn = legendre_value(n, sin_lat);
    Pn_dot = legendre_deriv(n, sin_lat);

    acc_point_zonal = acc_point_zonal + Jn(ind_zonal) * ((a/r)^n) ...
              * [(n + 1) * Pn; 0; -cos_lat * Pn_dot];
end

acc_point_zonal = acc_point_zonal * (-mu / (r * r));

% Evaluate tesseral harmonics.
for ind_tesseral = 1:number_tesseral
    n    = CSnm(ind_tesseral, 1);
    m    = CSnm(ind_tesseral, 2);
    C_nm = CSnm(ind_tesseral, 3);
    S_nm = CSnm(ind_tesseral, 4);

    cos_mlon = cos(m * lon_point);
    sin_mlon = sin(m * lon_point);

    % Associated Legendre function of degree n and order m and its 
    % derivative.
    
    % TBD: The Moshier implementation [3] seems to contain (-1)^m term for 
    % associated Legendre functions and their derivatives. These are not
    % present in [1] and [2] and the associated Legendre functions are 
    % computed correctly when compared to the Matlab implementation.
    P_nm = ((-1)^m) * legendre_assoc(n, m, sin_lat);
    P_nm_dot = ((-1)^m) * legendre_assocd(n, m, sin_lat);

    acc_point_tesseral = acc_point_tesseral + ((a / r)^n) * ...
        [   -(n + 1) * P_nm     * ( C_nm * cos_mlon + S_nm * sin_mlon); ...
         (m/cos_lat) * P_nm     * (-C_nm * sin_mlon + S_nm * cos_mlon); ...
             cos_lat * P_nm_dot * ( C_nm * cos_mlon + S_nm * sin_mlon)];
end
acc_point_tesseral = acc_point_tesseral * (-mu / (r * r));

acc_point = acc_point + acc_point_zonal;
acc_point = acc_point + acc_point_tesseral;

% Rotate back to the body coordinates.
acc_point = matrix_rot3(-lon_point) * matrix_rot2(lat_point) * acc_point;

% Torque per unit mass is r x acc_point.
T_body = cross(r_point, acc_point);
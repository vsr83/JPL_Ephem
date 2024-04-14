function [phi2, theta2, psi2] = libration_moon(phi, theta, psi, phi1, theta1, psi1, N)
% LIBRATION_MOON - Compute the second derivatives for the Moon libration
% angles.
%
% This method computes the expression in equation (3) of [1] or (8.6)-(8.8) 
% in [2].
%
% Important: This method has not been implemented for performance but for 
% simplicity and helping with an implementation. 
%
% INPUTS: 
%   phi            The clockwise angle along the xy-plane to the line of  
%                  nodes from the x axis (rad).
%   theta          The clockwise inclination of the body equator (rad).
%   psi            The clockwise angle from the node to the prime meridian
%                  along the body equator (rad).
%   phi1           The time derivative of phi (rad/d).
%   theta1         The time derivative of theta (rad/d).
%   psi1           The time derivative of psi (rad/d).
%   N              Torque per unit mass in body coordinates (3x1, au^2/d^2)
%
% OUTPUTS:
%   phi2           The second derivative of phi.
%   theta2         The second derivative of theta.
%   psi2           The second derivative of psi.
%
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [3] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip
%  [4] Ferrari et. al. - Geophysical Parameters of the Earth-Moon System,
%  Journal of Geophysical Research, 1980.

% The moment of inertia factors C/MR^2 are dimensionless. For example, [4]
% gives C/MR^2 = 0.3905 +- 0.0023. The lunar radius in Table 4 of [4] is
% 1738 km = 1.1617812e-05 au. Thereafter, thee principal moment of inertia
% per unit mass C/M = 0.3905 * (1.1617812e-05^2) au^2 = 5.2707e-11 au^2.

% Three principal moments of inertia per unit mass (au^2) [3].
A = 5.2699461151199766018387453E-11;
B = 5.2711485389894113965222358E-11;
C = 5.2732758299330413853212309E-11;

%A = 5.297704171694732e-11;
%B = 5.298912929015983e-11;
%C = 5.301051424905749e-11;

% Section 2D of [1].
beta_L = (C - A) / B;
gamma_L = (B - A) / C;
%beta_L = 6.3102131934887270E-04;
%gamma_L = 2.2773171480091860E-04;

% Angular velocity vector.
omega_x = phi1 * sin(theta) * sin(psi) + theta1 * cos(psi);
omega_y = phi1 * sin(theta) * cos(psi) - theta1 * sin(psi);
omega_z = phi1 * cos(theta) + psi1;

% Differential equations for the angular velocity from Euler's equations.
omega_x1 = omega_y * omega_z * (gamma_L - beta_L) / (1 - beta_L * gamma_L) ...
         + N(1) / A;
omega_y1 = omega_z * omega_x * beta_L + N(2) / B;
omega_z1 = -omega_x * omega_y * gamma_L + N(3) / C;

% Differential equations for the three Euler angles.
phi2 = (omega_x1 * sin(psi) + omega_y1 * cos(psi) + ...
        theta1 * (psi1 - phi1 * cos(theta))) / sin(theta);
theta2 = omega_x1 * cos(psi) - omega_y1 * sin(psi) - phi1 * psi1 * sin(theta);
psi2 = omega_z1 - phi2 * cos(theta) + phi1 * theta1 * sin(theta);

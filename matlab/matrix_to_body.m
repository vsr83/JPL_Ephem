function M = matrix_to_body(phi, theta, psi)
% MATRIX_TO_BODY - Compute Transformation to Body Coordinates.
%   
% INPUTS: 
%   phi            The clockwise angle along the xy-plane to the line of  
%                  nodes from the x axis (radians).
%   theta          The clockwise inclination of the body equator (radians).
%   psi            The clockwise angle from the node to the prime meridian
%                  along the body equator.
% OUTPUTS:
%   M              Rotation matrix to body coordinates (3x3).
% REFERENCES: 
%  [1] Newhall, Standish, Williams - DE 102: a numerically integrated
%  ephemeris of the Moon and planets spanning forty-four centuries,
%  Astronomy and Astrophysics, 125, 150-167, 1983.
%  [2] Urban, Seidelmann - Explanatory Supplement to the Astronomical
%  Almanac, 3rd edition, University Science Books, 2013.
%  [3] Steve Moshier, DE118i available at 
%  http://www.moshier.net/de118i-2.zip

M = matrix_rot3(psi) * matrix_rot1(theta) * matrix_rot3(phi);

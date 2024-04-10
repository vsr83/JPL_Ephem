function M = matrix_rot3d(angled)
% MATRIX_ROT3 - Compute matrix for rotation w.r.t. third coordinate axis in
% clockwise direction.
%
% INPUTS:
%   angled     The clockwise rotation angle (in degrees).
%
% OUTPUTS:
%   M          The rotation matrix.
%

M =  [cosd(angled), sind(angled), 0; -sind(angled), cosd(angled), 0; 0, 0, 1];

end
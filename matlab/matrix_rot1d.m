function M = matrix_rot1d(angled)
% MATRIX_ROT1 - Compute matrix for rotation w.r.t. first coordinate axis in
% clockwise direction.
%
% INPUTS:
%   angled     The clockwise rotation angle (in degrees).
%
% OUTPUTS:
%   M          The rotation matrix.
%

M = [1, 0, 0; 0, cosd(angled), sind(angled); 0, -sind(angled), cosd(angled)];

end
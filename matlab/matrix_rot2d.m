function M = matrix_rot2(angled)
% MATRIX_ROT2 - Compute matrix for rotation w.r.t. second coordinate axis in
% clockwise direction.
%
% INPUTS:
%   angled     The clockwise rotation angle (in degrees).
%
% OUTPUTS:
%   M          The rotation matrix.
%

M = [cosd(angled), 0, -sind(angled); 0, 1, 0; sind(angled), 0, cosd(angled)];

end
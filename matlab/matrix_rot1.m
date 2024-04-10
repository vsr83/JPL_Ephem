function M = matrix_rot1(angle)
% MATRIX_ROT1 - Compute matrix for rotation w.r.t. first coordinate axis in
% clockwise direction.
%
% INPUTS:
%   angle      The clockwise rotation angle (in radians).
%
% OUTPUTS:
%   M          The rotation matrix.
%

M = [1, 0, 0; 0, cos(angle), sin(angle); 0, -sin(angle), cos(angle)];

end
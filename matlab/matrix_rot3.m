function M = matrix_rot3(angle)
% MATRIX_ROT3 - Compute matrix for rotation w.r.t. third coordinate axis in
% clockwise direction.
%
% INPUTS:
%   angle      The clockwise rotation angle (in radians).
%
% OUTPUTS:
%   M          The rotation matrix.
%

M =  [cos(angle), sin(angle), 0; -sin(angle), cos(angle), 0; 0, 0, 1];

end
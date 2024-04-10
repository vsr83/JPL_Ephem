function M = matrix_rot2(angle)
% MATRIX_ROT2 - Compute matrix for rotation w.r.t. second coordinate axis in
% clockwise direction.
%
% INPUTS:
%   angle      The clockwise rotation angle (in radians).
%
% OUTPUTS:
%   M          The rotation matrix.
%

M = [cos(angle), 0, -sin(angle); 0, 1, 0; sin(angle), 0, cos(angle)];

end
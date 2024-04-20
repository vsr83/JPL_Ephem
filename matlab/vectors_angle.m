function a = vectors_angle(u, v);
% VECTORS_ANGLE - Compute the angle between two vectors.
%
% This method computes the angle between two vectors with more accurate
% method than the typical combination of dot product and acos.
%
% INPUTS:
%   u          The first vector (3 x 1 or 1 x 3).
%   v          The second vector (3 x 1 or 1 x 3).
%
% OUTPUTS:
%   a          The angle in degrees.
%
% REFERENCES: 
%  [1] https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors

a = norm(u);
b = norm(v);

if a < b
    d = a;
    a = b;
    b = d;
end

c = norm(u - v);

if b >= c && c >= 0
    mu = c - (a - b);
elseif c > b && b >= 0
    mu = b - (a - c);
else
   % error 'foo';
    a = 0;
    return;
end

a = 2 * atand(sqrt((((a - b) + c) * mu) / ((a + (b + c)) * ((a - c) + b))));
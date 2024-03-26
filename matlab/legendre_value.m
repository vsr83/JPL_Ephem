function v = legendre_value(n, x)
% LEGENDRE_VALUE - Compute value of a Legendre polynomial.
%   
% INPUTS: 
%   n              The degree of the polynomial.
%   x              The parameter to the polynomial.
% OUTPUTS:
%   v              The value.
% REFERENCES: 
%  [1] https://en.wikipedia.org/wiki/Legendre_polynomials, Referenced
%  26.3.2024.

switch n
    case 0
        v = 1;
    case 1
        v = x;
    case 2
        v = 0.5 * (3 * x.^2 - 1);
    case 3
        v = 0.5 * (5 * x.^3 - 3 * x);
    case 4
        v = 0.125 * (35 * x.^4 - 30 * x.^2 + 3);
    case 5
        v = 0.125 * (63 * x.^5 - 70 * x.^3 + 15 * x);
    case 6
        v = 0.0625 * (231 * x.^6 - 315 * x.^4 + 105 * x.^2 - 5);
end
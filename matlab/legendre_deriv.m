function v = legendre_deriv(n, x)
% LEGENDRE_DERIV - Compute derivative of a Legendre polynomial.
%   
% INPUTS: 
%   n              The degree of the polynomial.
%   x              The parameter to the polynomial.
% OUTPUTS:
%   v              The derivative.
% REFERENCES: 
%  [1] https://en.wikipedia.org/wiki/Legendre_polynomials, Referenced
%  26.3.2024.

switch n
    case 0
        v = 0;
    case 1
        v = 1;
    case 2
        v = 3 * x;
    case 3
        v = 7.5 * x.^2 - 1.5;
    case 4
        v = 17.5 * x.^3 - 7.5 * x;
    case 5
        v = 39.375 * x.^4 - 26.25 * x.^2 + 1.875;
    case 6
        v = 86.625 * x.^5 - 78.75 * x.^3 + 13.125 * x;
end
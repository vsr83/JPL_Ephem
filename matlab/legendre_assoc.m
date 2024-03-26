function v = legendre_assoc(l, m, x)
% LEGENDRE_ASSOC - Compute value of associated Legendre polynomial.
%   
% INPUTS: 
%   l              The degree of the polynomial.
%   m              The order of the polynomial.
%   x              The parameter to the polynomial.
% OUTPUTS:
%   v              The value.
% REFERENCES: 
%  [1] https://en.wikipedia.org/wiki/Associated_Legendre_polynomials, 
%  Referenced 26.3.2024.

switch l
    case 1
        switch m
            case 0
                v = x;
            case 1
                v = -sqrt(1 - x.^2);
        end
    case 2
        switch m
            case 0
                v = 0.5 * (3 * x.^2 - 1);
            case 1
                v = -3 * x .* sqrt(1 - x.^2);
            case 2
                v = 3 * (1 - x.^2);
        end
    case 3
        switch m
            case 0
                v = 0.5 * (5 * x.^3 - 3 * x);
            case 1
                v = 1.5 * (1 - 5 * x.^2) .* sqrt(1 - x.^2);
            case 2
                v = 15 * x .* (1 - x.^2);
            case 3
                v = -15 * sqrt(1 - x.^2).^3;
        end
    case 4
        switch m
            case 0
                v = 0.125 * (35 * x.^4 - 30 * x.^2 + 3);
            case 1
                v = -2.5 * (7 * x.^3 - 3 * x) .* sqrt(1 - x.^2);
            case 2
                v = 7.5 * (7 * x.^2 - 1) .* (1 - x.^2);
            case 3
                v = -105 * x .* sqrt(1 - x.^2).^3;
            case 4
                v = 105 * (1 - x.^2).^2;
        end
end
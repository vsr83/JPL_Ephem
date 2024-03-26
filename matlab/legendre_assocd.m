function v = legendre_assocd(l, m, x)
% LEGENDRE_ASSOCD - Compute derivative of associated Legendre polynomial.
%   
% INPUTS: 
%   l              The degree of the polynomial.
%   m              The order of the polynomial.
%   x              The parameter to the polynomial.
% OUTPUTS:
%   v              The derivative.
% REFERENCES: 
%  [1] https://en.wikipedia.org/wiki/Associated_Legendre_polynomials, 
%  Referenced 26.3.2024.

switch l
    case 1
        switch m
            case 0
                v = 1;
            case 1
                v = x ./ sqrt(1 - x.^2);
        end
    case 2
        switch m
            case 0
                v = 3 * x;
            case 1
                v = -3 * sqrt(1 - x.^2) + 3 * x.^2 ./ sqrt(1 - x.^2);
            case 2
                v = -6 * x;
        end
    case 3
        switch m
            case 0
                v = 7.5 * x.^2 - 1.5;
            case 1
                v = (-1.5 + 7.5 * x.^2) .* x ./ sqrt(1 - x.^2) - 15 * x .* sqrt(1 - x.^2);
            case 2
                v = 15 - 45 * x.^2;
            case 3
                v = 45 * x .* sqrt(1 - x.^2);
        end
    case 4
        switch m
            case 0
                v = 17.5 * x.^3 - 7.5 * x;
            case 1
                v = (17.5 * x.^4 - 7.5 * x.^2) ./ sqrt(1 - x.^2) ...
                  - (52.5 * x.^2 - 7.5) .* sqrt(1 - x.^2);
            case 2
                v = 120 * x - 210 * x.^3;
            case 3
                v = -105 * sqrt(1 - x.^2).^3 + 315 * x.^2 .* sqrt(1 - x.^2);
            case 4
                v = -420 * (x - x.^3);
        end
end
function [y_out, t_out, F_out] = adams8(func_in, t_in, y_in, F_in, h, param)
% RUNGE4 - Perform numerical integration step with 8th order 
% Adams-Bashforth-Moulton method.
%
% This method performs a single integration step for the initial value
% problem
%    dy/dt = f(t, y) 
%    y(t_in) = y_in
%
% INPUTS:
%   func_in     Function handle for f(t, y).
%   t_in        Time before the step.
%   y_in        Previous solution.
%   F_in        Eight states f(t, y) before the step.
%   h           Time step size.
%   param       Parameter passed to func_in.
%
% OUTPUTS:
%   Y_out       Four states after the step.
%   t_out       Time after the step.
%   F_out       The eight states for the next step.

f_0  = F_in(:, 1);
f_m1 = F_in(:, 2);
f_m2 = F_in(:, 3);
f_m3 = F_in(:, 4);
f_m4 = F_in(:, 5);
f_m5 = F_in(:, 6);
f_m6 = F_in(:, 7);
f_m7 = F_in(:, 8);

% Predictor
y_new = y_in + (h / 120960) * (434241 * f_0 - 1152169 * f_m1 + 2183877 * f_m2 ... 
    - 2664477 * f_m3 + 2102243 * f_m4 - 1041723 * f_m5 + 295767 * f_m6 - 36799 * f_m7); 

f_1 = func_in(t_in + h, y_new, param);
% Corrector
y_out = y_in + (h / 120960) * (36799 * f_1 + 139849 * f_0 - 121797 * f_m1 ...
    + 123133 * f_m2 - 88547 * f_m3 + 41499 * f_m4 - 11351 * f_m5 + 1375 * f_m6);

f_1 = func_in(t_in + h, y_out, param);

F_out = [f_1, f_0, f_m1, f_m2, f_m3, f_m4, f_m5, f_m6];
t_out = t_in + h;

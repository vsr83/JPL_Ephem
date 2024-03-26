function [y_out, t_out] = runge4(func_in, t_in, y_in, h, param)
% RUNGE4 - Perform numerical integration step with 4th order Runge-Kutta.
%
% This method performs a single integration step for the initial value
% problem
%    dy/dt = f(t, y) 
%    y(t_in) = y_in
%
% INPUTS:
%   func_in     Function handle for f(t, y).
%   t_in        Time before the step.
%   y_in        State before the step.
%   h           Time step size.
%   param       Parameter passed to func_in.
%
% OUTPUTS:
%   y_out       State after the step.
%   t_out       Time after the step.

k_1 = func_in(t_in, y_in, param);
k_2 = func_in(t_in + h/2, y_in + h * k_1/2, param);
k_3 = func_in(t_in + h/2, y_in + h * k_2/2, param);
k_4 = func_in(t_in + h,   y_in + h * k_3, param);

y_out = y_in + (h / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
t_out = t_in + h;
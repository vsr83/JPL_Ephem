function [y_out, t_out] = runge8(func_in, t_in, y_in, h, param)
% RUNGE8 - Perform numerical integration step with 8th order Runge-Kutta.
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

t_2 = t_in + h * (4.0/27.0);
t_3 = t_in + h * (2.0/9.0);
t_4 = t_in + h * (1.0/3.0);
t_5 = t_in + h * (1.0/2.0);
t_6 = t_in + h * (2.0/3.0);
t_7 = t_in + h * (1.0/6.0);
t_8 = t_in + h;
t_9 = t_in + h * (5.0/6.0);
t_10= t_in + h;

k_1 = func_in(t_in, y_in, param);
y_2 = y_in + h*(4.0/27.0) * k_1;
k_2 = func_in(t_2, y_2, param);
y_3 = y_in + h*(1.0/18.0) * (k_1 + 3 * k_2);
k_3 = func_in(t_3, y_3, param);
y_4 = y_in + h*(1.0/12.0) * (k_1 + 3 * k_3);
k_4 = func_in(t_4, y_4, param);
y_5 = y_in + h*(1.0/8.0) * (k_1 + 3 * k_4);
k_5 = func_in(t_5, y_5, param);
y_6 = y_in + h*(1.0/54.0) * (13*k_1 - 27*k_3 + 42*k_4 + 8*k_5);
k_6 = func_in(t_6, y_6, param);
y_7 = y_in + h*(1.0/4320.0) * (389*k_1 - 54*k_3 + 966*k_4 - 824*k_5 + 243*k_6);
k_7 = func_in(t_7, y_7, param);
y_8 = y_in + h*(1.0/20.0) * (-234*k_1 + 81*k_3 - 1164*k_4 + 656*k_5 - 122*k_6 + 800*k_7);
k_8 = func_in(t_8, y_8, param);
y_9 = y_in + h*(1.0/288.0) * (-127*k_1 + 18*k_3 - 678*k_4 + 456*k_5 - 9*k_6 + 576*k_7 + 4*k_8);
k_9 = func_in(t_9, y_9, param);
y_10 = y_in + h*(1.0/820.0) * (1481*k_1 - 81*k_3 + 7104*k_4 - 3376*k_5 + 72*k_6 - 5040*k_7 - 60*k_8 + 720*k_9);
k_10 = func_in(t_10, y_10, param);

y_out = y_in + (h/840) * (41*k_1 + 27*k_4 + 272*k_5 + 27*k_6 + 216*k_7 + 216*k_9 + 41*k_10);
t_out = t_in + h;
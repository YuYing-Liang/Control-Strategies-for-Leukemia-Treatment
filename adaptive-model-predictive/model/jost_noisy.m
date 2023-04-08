function [f] = jost_noisy(t,x,d,theta,num_days_for_cycle,on_drug_days,var_w,t_s)
%JOST state ODE for 6-MP and 6-TGN leukopoiesis with system noise
%   8 state model described in Jost et al. 2020
%   d = dose
%   theta = vector of patient parameters (4 by 1)
%   var_w = noise variance

if get_cycle_day(t, num_days_for_cycle) <= on_drug_days
    u = d;
else
    u = 0;
end

bsa = 1.71; % change to match patient!
theta = [
    31.2;
    12.72;
    0.019;
    9.9216;
    0.219*(bsa^1.16);
    theta(1);
    theta(2);
    theta(3);
    theta(4);
    2.3765;
];

% define matrices and nonlinear components
A = [-theta(1)                0 0 0 0 0 0 0; ...
    theta(1) -theta(2)          0 0 0 0 0 0; ...
    0 theta(3)*theta(4) -theta(5) 0 0 0 0 0; ...
    0 0 0 -theta(7)                 0 0 0 0; ...
    0 0 0 theta(7) -theta(7)          0 0 0; ...
    0 0 0 0 theta(7) -theta(7)          0 0; ...
    0 0 0 0 0 theta(7) -theta(7)          0; ...
    0 0 0 0 0 0 theta(7) -theta(10)];

B = [0.22;0;0;0;0;0;0;0];

f_hat = [...
    0;...
    0;...
    0;...
    ((theta(6)/(x(8)+0.00001))^theta(9))*(theta(7)*x(4)-theta(7)*theta(8)*x(3)*x(4));...
    0;...
    0;...
    0;...
    0];

f = A*x + B*u + f_hat;
f = x + t_s * f + normrnd(0,var_w,[8,1]);
end
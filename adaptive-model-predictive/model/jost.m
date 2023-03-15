function [f] = jost(t,x,d,theta,num_days_for_cycle,on_drug_days)
%JOST state ODE for 6-MP and 6-TGN leukopoiesis
%   8 state model described in Jost et al. 2020
%   d = dose
%   theta = vector of parameters

if mod(t,num_days_for_cycle) <= on_drug_days
    u = d;
else
    u = 0;
end


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

f_hat = [0;...
    0;...
    0;...
    theta(7)*((theta(6)/(x(8)+1e-6))^theta(9))*x(4)-theta(7)*theta(8)*((theta(6)/(x(8)+1e-6))^theta(9))*x(3)*x(4);...
    0;...
    0;...
    0;...
    0];

f = A*x + B*u + f_hat;

end
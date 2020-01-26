% Description: Numerical Approximation of RLC Circuit
% Author: John Nicholas
% This script performs a numerical approximation of current
% in an RLC circuit with given initial conditions using a
% step size of dt over N steps.
% The second-order differential equation defining the system:
%     i'' + (R/L)i' + (1/LC)i = 0
% Is simplified into a system of equations:
%     i' = j
%     j' = -(R/L)j - (1/LC)i
% Using Euler's method, we approximate j(t) and i(t)
%     j(t+dt) ~ j(t) + j'(t)*dt
%     i(t+dt) ~ i(t) + i'(t)*dt
% Substituting the above expressions for i' and j':
%     j(t+dt) ~ j(t) + [-(R/L)j - (1/LC)i]*dt
%     i(t+dt) ~ i(t) + i(t)*dt

clear;
% Set Up Initial values
R = 5; % Ohms
L = 10; % Henrys
C = 100; % Farads
N = 10000; % Number of steps
i(N) = 0; % Initialize array
j(N) = 0; % Initialize array
i(1) = 5; % Amperes
j(1) = 0; % Amperes per second
dt = 0.001; % s

for t = 1:N
   j(t+1) = j(t) + (-(R/L)*j(t) - (1/L*C)*i(t))*dt; % Euler's method for j(t)
   i(t+1) = i(t) + j(t)*dt; % Euler's method for i(t)
end

plot(0:dt:N*dt,i,'k-.','LineWidth',1);
title('Current Through RLC Circuit');
xlabel('Time (s)');
ylabel('Current (A)');

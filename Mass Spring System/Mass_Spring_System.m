% Description: Numerical Approximation of Mass-Spring system
% Author: John Nicholas
% This script performs a numerical approximation of the position
% and velocity of an idea mass-spring system with given initial
% conditions using a step size of 'dt' over N steps.
% The second-order differential equation defining the system:
%     x'' = -kx/m
% Is simplified into a system of equations:
%     x' = v
%     v' = -kx/m
% Using Euler's method, we approximate v(t) and x(t)
%     v(t+dt) ~ v(t) + v'(t)*dt
%     x(t+dt) ~ x(t) + x'(t)*dt
% Substituting the above expressions for x' and v':
%     v(t+dt) ~ v(t) + -kx(t)*dt/m
%     x(t+dt) ~ x(t) + v(t)*dt

clear;
% Set Up Initial values
k = 40; % N/m
m = 4; % kg
N = 5000; % Number of steps
x(N) = 0; % Initialize array
v(N) = 0; % Initialize array
x(1) = -0.1; % m ; Initial position
v(1) = 0; % m/s; Initial velocity
dt = 0.001; % s

for t = 1:N
   v(t+1) = v(t) + dt*(-k/m*x(t)); % Euler's method for v(t)
   x(t+1) = x(t) + v(t)*dt; % Euler's method for x(t)
end

figure;
subplot(2,1,1);
plot(0:dt:N*dt,x,'k-.','LineWidth',1);
xlabel('Time (s)');
ylabel('Position (m)');
title('Position of Spring-Mass System');
subplot(2,1,2);
plot(0:dt:N*dt,v,'r-.','LineWidth',1);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity of Spring-Mass System');

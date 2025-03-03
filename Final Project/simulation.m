% HW4 Ruigang Chen & Ben Sarfati
clear all; close all; clc

%% globals and parameters

% global sgn_slip
global mu
% sgn_slip = 1;
mu = 100;
rez = 0.01;

%ODE parameters
tspan = [0 10]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         
op_slip = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_slip);         

%% testing 
    
X0 = [0 0 pi/16 pi/4 0 0 0 0]';  %[x; y; th1; th2; dx; dy; dth1; dth2];
[t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 

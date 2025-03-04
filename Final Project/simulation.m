% Final Project Ruigang Chen & Ben Sarfati
clear all; close all; clc

%% globals and parameters

[~, ~, ~, l, ~, ~, ~] = model_params();
rez = 0.01;

%ODE parameters
tspan = [0 10]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         

%% testing 

%[x; y; th1; th2; dx; dy; dth1; dth2];
% X0 = [0 0 pi/16 pi/4 0 0 0 0]';
% X0 = [0 0 pi/16 deg2rad(89) 0 0 0 0]';
% X0 = [0 0 pi/8 -pi/3 0 0 0 0]';

%Initial guess from figure 9b of Gamus and Or 2015m fig. 9b
% X0 = [0 0 0.75 -0.5 0 0 0 0];

[t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 

%% collecting

x = X(:,1);
y = X(:,2);
th1 = X(:,3);
th2 = X(:,4);
x_d = X(:,5);
y_d = X(:,6);
th1_d = X(:,7);
th2_d = X(:,8);
xtilde = x+2*l*sin(th1)+2*l*sin(th2);
ytilde = y+2*l*cos(th1)-2*l*cos(th2);
xtilde_d = x_d+2*l*th1_d.*cos(th1)+2*l*th2_d.*cos(th2);
ytilde_d = y_d-2*l*th1_d.*sin(th1)+2*l*th2_d.*sin(th2);

%% plotting everything

close all;

figure;
h1 = plot(t,x,'LineWidth',2); hold on
h1 = plot(t,xtilde,'LineWidth',2); hold on
set(gcf,'color','w');
title('x vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('x', 'Interpreter', 'latex', 'fontsize', 20);
legend('$x$','$\tilde x$','Interpreter','latex','fontsize',20,'location','ne')
grid on;

figure;
h1 = plot(t,y,'LineWidth',2); hold on
h1 = plot(t,ytilde,'LineWidth',2); hold on
set(gcf,'color','w');
title('y vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('y', 'Interpreter', 'latex', 'fontsize', 20);
legend('$y$','$\tilde y$','Interpreter','latex','fontsize',20,'location','ne')
grid on;

figure;
h1 = plot(t,th1*180/pi,'LineWidth',2); hold on
h2 = plot(t,th2*180/pi,'LineWidth',2);
set(gcf,'color','w');
title('Angles vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Angle[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('$\theta_1$','$\theta_2$','Interpreter','latex','fontsize',20,'location','nw')
grid on;

% figure;
% h1 = plot(t,y_d,'LineWidth',2); hold on
% h1 = plot(t,ytilde_d,'LineWidth',2); hold on
% set(gcf,'color','w');
% title('$\dot y$ vs. Time','fontsize',20,'Interpreter','latex')
% xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
% ylabel('$\dot y$', 'Interpreter', 'latex', 'fontsize', 20);
% legend('$\dot y$','$\dot{\tilde y}$','Interpreter','latex','fontsize',20,'location','ne')
% grid on;


%% Testing with impact 

finalTime = 10;
finalTimes = [];
finalX = [];
impactTimes = [];
% finalLambda = [];
% finalq_dd = []; %for cross-checking

%Solve ODE up to terminal event
[t,X,te,ye,ie] = ode45(@sys_stick, [0,finalTime], X0, op_stick);  

%Record progress
finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))]; %transition state will be considered to belong to the next state
finalTimes = [finalTimes;t(1:end-1);NaN];
% finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))]; %NaNs added for plotting disconinuities
% finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];

while t(end) < finalTime
    switch ie(end)
        case 3
            warning('simulation interrupted due to falling');
            break;
        case 4
            %Record impact 
            impactTimes(end+1) = t(end);

            %Account for impact
            Xplus = impact_law(ye(end,:));
            
            %Solve ODE up to next terminal event
            [t,X,te,ye,ie] = ode45(@sys_stick, [te, finalTime], Xplus, op_stick);
            
            %Record progress
            finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))]; %transition state will be considered to belong to the next state
            finalTimes = [finalTimes;t(1:end-1);NaN];
            % finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))]; %NaNs added for plotting disconinuities
            % finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];

    end
end

%Replace final NaN
finalX(end,:) = X(end,:);
finalTimes(end,:) = t(end,:);
% finalLambda(end,:) = Lambda(end,:);
% finalq_dd(end,:) = q_dd(end,:);

X = finalX;
t = finalTimes;

%% collecting

x = X(:,1);
y = X(:,2);
th1 = X(:,3);
th2 = X(:,4);
x_d = X(:,5);
y_d = X(:,6);
th1_d = X(:,7);
th2_d = X(:,8);
xtilde = x+2*l*sin(th1)+2*l*sin(th2);
ytilde = y+2*l*cos(th1)-2*l*cos(th2);
xtilde_d = x_d+2*l*th1_d.*cos(th1)+2*l*th2_d.*cos(th2);
ytilde_d = y_d-2*l*th1_d.*sin(th1)+2*l*th2_d.*sin(th2);

%% plotting everything

close all;

figure;
h1 = plot(t,x,'LineWidth',2); hold on
h1 = plot(t,xtilde,'LineWidth',2); hold on
set(gcf,'color','w');
title('x vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('x', 'Interpreter', 'latex', 'fontsize', 20);
legend('$x$','$\tilde x$','Interpreter','latex','fontsize',20,'location','ne')
grid on;

figure;
h1 = plot(t,y,'LineWidth',2); hold on
h1 = plot(t,ytilde,'LineWidth',2); hold on
set(gcf,'color','w');
title('y vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('y', 'Interpreter', 'latex', 'fontsize', 20);
legend('$y$','$\tilde y$','Interpreter','latex','fontsize',20,'location','ne')
grid on;

figure;
h1 = plot(t,th1*180/pi,'LineWidth',2); hold on
h2 = plot(t,th2*180/pi,'LineWidth',2);
set(gcf,'color','w');
title('Angles vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Angle[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('$\theta_1$','$\theta_2$','Interpreter','latex','fontsize',20,'location','nw')
grid on;

% figure;
% h1 = plot(t,y_d,'LineWidth',2); hold on
% h1 = plot(t,ytilde_d,'LineWidth',2); hold on
% set(gcf,'color','w');
% title('$\dot y$ vs. Time','fontsize',20,'Interpreter','latex')
% xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
% ylabel('$\dot y$', 'Interpreter', 'latex', 'fontsize', 20);
% legend('$\dot y$','$\dot{\tilde y}$','Interpreter','latex','fontsize',20,'location','ne')
% grid on;
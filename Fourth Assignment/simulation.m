% HW2 Q3 Ruigang Chen & Ben Sarfati
clear all; close all; clc

%% globals and parameters

global sgn_slip
global mu
global R

sgn_slip = 1;
mu = 0.05;
R = 0.6;

rez = 0.01;

mu = 0.05;

%ODE parameters
tspan = [0 10]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         
op_slip = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_slip);         


%% Q5

%run
te = [];
wSlip = 0;
while isempty(te)
    wSlip = wSlip+rez;
    X0 = [0 R 0 0 0 0 0 wSlip];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
    [t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 
end
w0 = 0.9*wSlip;
X0 = [0; R; 0; 0; 0; 0; 0; w0];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
[t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 
Lambda = zeros(length(t),2);
for i = 1:length(t)
    [~,Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
end

%label
th = X(:,3);
ph = X(:,4);
lambdan = Lambda(:,2);
lambdat = Lambda(:,1);

%plot (a)
figure;
h1 = plot(t,th*180/pi,'LineWidth',2); hold on
h2 = plot(t,ph*180/pi,'LineWidth',2);
% plot(te,the*180/pi,'.','Color',h1.Color,'MarkerSize',25)
% plot(te,phe*180/pi,'.','Color',h2.Color,'MarkerSize',25)
set(gcf,'color','w');
title('Angles vs. Time; $0.9\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 30);
legend('$\theta $','$\phi $','Interpreter','latex','fontsize',30,'location','ne')
xlim(tspan)
grid on;
% saveas(gcf, 'q5astick.png');

%plot (c)
figure;
plot(t,lambdan,'LineWidth',2);
set(gcf,'color','w');
title('Normal Force vs. Time; $0.9\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_n$ [N]', 'Interpreter', 'latex', 'fontsize', 30);
xlim(tspan)
grid on;
% saveas(gcf, 'q5cstick.png');

%plot (d)
figure;
plot(t,lambdat./lambdan,'LineWidth',2);
set(gcf,'color','w');
title('Force ratio vs. Time; $0.9\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_t}{\lambda_n}$', 'Interpreter', 'latex', 'fontsize', 30);
xlim(tspan)
grid on;
yline(mu)
yline(-mu)
legend('ratio','ratio limits','fontsize',30,'location','ne')
% saveas(gcf, 'q5dstick.png');

w0 = wSlip;
X0 = [0; R; 0; 0; 0; 0; 0; w0];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
[t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 
Lambda = zeros(length(t),2);
for i = 1:length(t)
    [~,Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
end

%label
th = X(:,3);
the = ye(3);
ph = X(:,4);
phe = ye(4);
lambdan = Lambda(:,2);
lambdat = Lambda(:,1);

%plot (a)
figure;
h1 = plot(t,th*180/pi,'LineWidth',2); hold on
h2 = plot(t,ph*180/pi,'LineWidth',2);
plot(te,the*180/pi,'.','Color',h1.Color,'MarkerSize',25)
plot(te,phe*180/pi,'.','Color',h2.Color,'MarkerSize',25)
set(gcf,'color','w');
title('Angles vs. Time; $\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 30);
legend('$\theta $','$\phi $','$\theta_{slip} $','$\phi_{slip} $','Interpreter','latex','fontsize',30,'location','ne')
xlim(tspan)
grid on;
% saveas(gcf, 'q5aslip.png');

%plot (c)
figure;
plot(t,lambdan,'LineWidth',2);
set(gcf,'color','w');
title('Normal Force vs. Time; $\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_n$ [N]', 'Interpreter', 'latex', 'fontsize', 30);
xlim(tspan)
grid on;
% saveas(gcf, 'q5cslip.png');

%plot (d)
figure;
plot(t,lambdat./lambdan,'LineWidth',2);
set(gcf,'color','w');
title('Force ratio vs. Time; $\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_t}{\lambda_n}$', 'Interpreter', 'latex', 'fontsize', 30);
xlim(tspan)
grid on;
yline(mu)
yline(-mu)
legend('ratio','ratio limits','fontsize',30,'location','ne')
% saveas(gcf, 'q5dslip.png');

% %plot (b)
% figure;
% h1 = plot(t,vt,'LineWidth',2);
% set(gcf,'color','w');
% title('Slip Velocity vs. Time; $0.9\omega_{slip}$','fontsize',20,'Interpreter','latex')
% xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
% ylabel('[$V_t$]', 'Interpreter', 'latex', 'fontsize', 30);
% xlim(tspan)
% grid on;
% % saveas(gcf, 'q5bstick.png');



%% Q6

% wBreak = wSlip;
wBreak = 1.32;
reachedEndTime = 1;
while reachedEndTime
    wBreak = wBreak+rez;
    [t,X,Lambda,q_dd,stickInds,reachedEndTime] = solveHybridDynamics(wBreak);
    wbreak
end

%label
th = X(:,3);
ph = X(:,4);
lambdan = Lambda(:,2);
lambdat = Lambda(:,1);
vt = X(:,5)+R*X(:,7);

%plot (a)
figure;
h1 = plot(t(stickInds),th(stickInds)*180/pi,'LineWidth',2); hold on
h2 = plot(t(stickInds),ph(stickInds)*180/pi,'LineWidth',2);
plot(t(~stickInds),th(~stickInds)*180/pi,'--','Color',h1.Color)
plot(t(~stickInds),ph(~stickInds)*180/pi,'--','Color',h2.Color)
set(gcf,'color','w');
title('Angles vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 30);
legend('$\theta_{stick}$','$\phi_{stick}$','$\theta_{slip}$','$\phi_{slip}$','Interpreter','latex','fontsize',30,'location','ne')
grid on;
% saveas(gcf, 'q6a.png');

%plot (b)
figure;
h1 = plot(t(stickInds),vt(stickInds),'LineWidth',2); hold on
plot(t(~stickInds),vt(~stickInds),'--','Color',h1.Color)
set(gcf,'color','w');
title('Slip Velocity vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('[$V_t$]', 'Interpreter', 'latex', 'fontsize', 30);
legend('$v_{t,stick}$','$v_{t,slip}$','Interpreter','latex','fontsize',30,'location','ne')
grid on;
% saveas(gcf, 'q6b.png');

%plot (c)
figure;
h1 = plot(t(stickInds),lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Normal Force vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_{n}$ [N]','Interpreter', 'latex', 'fontsize', 30);
legend('stick','slip','Interpreter','latex','fontsize',30,'location','ne')
grid on;
% saveas(gcf, 'q6c.png');

%plot (d)
figure;
h1 = plot(t(stickInds),lambdat(stickInds)./lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdat(~stickInds)./lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Force ratio vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_{t}}{\lambda_{n}}$', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
yline(mu)
yline(-mu)
legend('stick','slip','ratio limits','fontsize',30,'location','ne')
% saveas(gcf, 'q6d.png');

%% Q7

w0 = (3*wSlip+wBreak)/4;



%label
t = finalTimes;
th = finalX(:,3);
ph = finalX(:,4);
lambdan = finalLambda(:,2);
lambdat = finalLambda(:,1);
vt = finalX(:,5)+R*finalX(:,7);

%plot (a)
figure;
h1 = plot(t(stickInds),th(stickInds)*180/pi,'LineWidth',2); hold on
h2 = plot(t(stickInds),ph(stickInds)*180/pi,'LineWidth',2);
plot(t(~stickInds),th(~stickInds)*180/pi,'--','Color',h1.Color)
plot(t(~stickInds),ph(~stickInds)*180/pi,'--','Color',h2.Color)
set(gcf,'color','w');
title('Angles vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 30);
legend('$\theta_{stick}$','$\phi_{stick}$','$\theta_{slip}$','$\phi_{slip}$','Interpreter','latex','fontsize',30,'location','ne')
grid on;
% saveas(gcf, 'q6a.png');

%plot (b)
figure;
h1 = plot(t(stickInds),vt(stickInds),'LineWidth',2); hold on
plot(t(~stickInds),vt(~stickInds),'--','Color',h1.Color)
set(gcf,'color','w');
title('Slip Velocity vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('[$V_t$]', 'Interpreter', 'latex', 'fontsize', 30);
legend('$v_{t,stick}$','$v_{t,slip}$','Interpreter','latex','fontsize',30,'location','ne')
grid on;
% saveas(gcf, 'q6b.png');

%plot (c)
figure;
h1 = plot(t(stickInds),lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Normal Force vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_{n}$ [N]','Interpreter', 'latex', 'fontsize', 30);
legend('stick','slip','Interpreter','latex','fontsize',30,'location','ne')
grid on;
% saveas(gcf, 'q6c.png');

%plot (d)
figure;
h1 = plot(t(stickInds),lambdat(stickInds)./lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdat(~stickInds)./lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Force ratio vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_{t}}{\lambda_{n}}$', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
yline(mu)
yline(-mu)
legend('stick','slip','ratio limits','fontsize',30,'location','ne')
% saveas(gcf, 'q6d.png');

%plot (d)
figure;
h1 = plot(t(stickInds),lambdat(stickInds)./lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdat(~stickInds)./lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Force ratio vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_{t}}{\lambda_{n}}$', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
yline(mu)
yline(-mu)
legend('stick','slip','ratio limits','fontsize',30,'location','ne')
% saveas(gcf, 'q6d.png');
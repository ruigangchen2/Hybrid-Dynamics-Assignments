% HW4 Ruigang Chen & Ben Sarfati
clear all; close all; clc

%% globals and parameters

global sgn_slip
global mu
sgn_slip = 1;
mu = 0.05;
rez = 0.01;
R = 0.6;

%ODE parameters
tspan = [0 10]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         
op_slip = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_slip);         


te = [];
wSlip = 0;
while isempty(te)
    wSlip = wSlip+rez;
    X0 = [0 R 0 0 0 0 0 wSlip];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
    [t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 
end

%% Q5 (ω0 = 0.9ωslip)
w0 = 0.9*wSlip;
X0 = [0; R; 0; 0; 0; 0; 0; w0];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
[t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 
Lambda = zeros(length(t),2);
for i = 1:length(t)
    [~,Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
end

th = X(:,3);
ph = X(:,4);
lambdat = Lambda(:,1);
lambdan = Lambda(:,2);

%% Q5 (ω0 = 0.9ωslip) --- plot a)
figure;
h1 = plot(t,th*180/pi,'LineWidth',2); hold on
h2 = plot(t,ph*180/pi,'LineWidth',2);
set(gcf,'color','w');
title('Angles vs. Time; $0.9\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Angle[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('$\theta(t)$','$\phi(t)$','Interpreter','latex','fontsize',20,'location','ne')
xlim(tspan)
grid on;
saveas(gcf, 'q5a(ω0 = 0.9ωslip).png');

%% Q5 (ω0 = 0.9ωslip) --- plot c)
figure;
plot(t,lambdan,'LineWidth',2);
set(gcf,'color','w');
title('Normal Force vs. Time; $0.9\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_n$ [N]', 'Interpreter', 'latex', 'fontsize', 20);
xlim(tspan)
grid on;
saveas(gcf, 'q5c(ω0 = 0.9ωslip).png');

%% Q5 (ω0 = 0.9ωslip) --- plot d)
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
ylim([-0.1,0.1])
legend('$\frac{\lambda_t}{\lambda_n}$','$\mu$ limits','Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'q5d(ω0 = 0.9ωslip).png');

%% Q5 (ω0 = ωslip)
w0 = wSlip-0.01 % The maximum value at the current resolution
X0 = [0; R; 0; 0; 0; 0; 0; w0];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
[t,X,te,ye,ie] = ode45(@sys_stick, t_eval, X0, op_stick); 
Lambda = zeros(length(t),2);
for i = 1:length(t)
    [~,Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
end

th = X(:,3);
ph = X(:,4);
lambdan = Lambda(:,2);
lambdat = Lambda(:,1);

%% Q5 (ω0 = ωslip) --- plot a)
figure;
h1 = plot(t,th*180/pi,'LineWidth',2); hold on
h2 = plot(t,ph*180/pi,'LineWidth',2);
set(gcf,'color','w');
title('Angles vs. Time; $\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Angle[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('$\theta(t)$','$\phi(t)$','Interpreter','latex','fontsize',20,'location','ne')
xlim(tspan)
grid on;
saveas(gcf, 'q5a(ω0 = ωslip).png');

%% Q5 (ω0 = ωslip) --- plot c)
figure;
plot(t,lambdan,'LineWidth',2);
set(gcf,'color','w');
title('Normal Force vs. Time; $\omega_{slip}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_n$ [N]', 'Interpreter', 'latex', 'fontsize', 20);
xlim(tspan)
grid on;
saveas(gcf, 'q5c(ω0 = ωslip).png');

%% Q5 (ω0 = ωslip) --- plot d)
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
ylim([-0.1,0.1])
legend('$\frac{\lambda_t}{\lambda_n}$','$\mu$ limits','Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'q5d(ω0 = ωslip).png');

%% Q6

% wBreak = wSlip;
wBreak = 19.7;
lookingForWbreak = 1;
finalTime = 10;
while lookingForWbreak
    wBreak = wBreak+rez
    te = 0;
    ye = [0 R 0 0 0 0 0 wBreak];
    finalX = [];
    finalLambda = [];
    finalTimes = [];
    finalq_dd = []; %for cross-checking
    stickInds = logical([]);
    
    while ~isempty(te)
        if events_stick(te,ye)<0
            %Solve corresponding ODE
            [t,X,te,ye,ie] = ode45(@sys_stick, [te, finalTime], ye, op_stick);
            Lambda = zeros(length(t),2);
            q_dd = zeros(length(t),4);
            for i = 1:length(t)
                [q_dd(i,:),Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
            end
            
            %Record results
            finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))]; %transition state will be considered to belong to the next state
            finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))]; %NaNs added for plotting disconinuities
            finalTimes = [finalTimes;t(1:end-1);NaN];
            finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];

            %Record and update state
            stickInds = [stickInds; true(length(t),1)];
            sgn_slip = sign(1.5-ie);
        else
            %Solve corresponding ODE
            [t,X,te,ye,ie] = ode45(@sys_slip, [te, finalTime], ye, op_slip);
            Lambda = zeros(length(t),2);
            q_dd = zeros(length(t),4);
            for i = 1:length(t)
                [q_dd(i,:),Lambda(i,2)] = dyn_sol_slip(X(i,1:4)',X(i,5:8)',t(i));
                Lambda(i,1) = -sgn_slip*mu*Lambda(i,2);
            end

            %Record results
            finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))];
            finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))];
            finalTimes = [finalTimes;t(1:end-1);NaN];
            finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];
           
            %Record and update state
            stickInds = [stickInds; false(length(t),1)];
            if ie-1
                lookingForWbreak = 0;
                break;
            end
            sgn_slip = -sgn_slip;
        end
    end
end
wBreak
finalX(end,:) = X(end,:);
finalLambda(end,:) = Lambda(end,:);
finalTimes(end,:) = t(end,:);
finalq_dd(end,:) = q_dd(end,:);


t = finalTimes;
th = finalX(:,3);
ph = finalX(:,4);
lambdan = finalLambda(:,2);
lambdat = finalLambda(:,1);
vt = finalX(:,5)+R*finalX(:,7);

%% Q6 plot a)
figure;
h1 = plot(t(stickInds),th(stickInds)*180/pi,'LineWidth',2); hold on
h2 = plot(t(stickInds),ph(stickInds)*180/pi,'LineWidth',2);
plot(t(~stickInds),th(~stickInds)*180/pi,'--','Color',h1.Color)
plot(t(~stickInds),ph(~stickInds)*180/pi,'--','Color',h2.Color)
set(gcf,'color','w');
title('Angles vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Angle[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('$\theta_{stick}$','$\phi_{stick}$','$\theta_{slip}$','$\phi_{slip}$','Interpreter','latex','fontsize',20,'location','ne')
grid on;
saveas(gcf, 'q6a.png');

%% Q6 plot b)
figure;
h1 = plot(t(stickInds),vt(stickInds),'LineWidth',2); hold on
plot(t(~stickInds),vt(~stickInds),'--','Color',h1.Color)
set(gcf,'color','w');
title('Slip Velocity vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$v_t [\frac{m}{s}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('stick','slip','Interpreter','latex','fontsize',20,'location','ne')
grid on;
saveas(gcf, 'q6b.png');

%% Q6 plot c)
figure;
h1 = plot(t(stickInds),lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Normal Force vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_{n}$ [N]','Interpreter', 'latex', 'fontsize', 20);
legend('stick','slip','Interpreter','latex','fontsize',20,'location','ne')
grid on;
saveas(gcf, 'q6c.png');

%% Q6 plot d)
figure;
h1 = plot(t(stickInds),lambdat(stickInds)./lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdat(~stickInds)./lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Force ratio vs. Time; $\omega_{break}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_{t}}{\lambda_{n}}$', 'Interpreter', 'latex', 'fontsize', 20);
grid on;
yline(mu)
yline(-mu)
ylim([-0.1,0.1])
legend('stick','slip','$\mu$ limits','Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'q6d.png');

%% Q7

w0 = (3*wSlip+wBreak)/4;
te = 0;
ye = [0 R 0 0 0 0 0 w0];
finalX = [];
finalLambda = [];
finalTimes = [];
finalq_dd = []; %for cross-checking
stickInds = logical([]);

while ~isempty(te)
    if events_stick(te,ye)<0
        %Solve corresponding ODE
        
        [t,X,te,ye,ie] = ode45(@sys_stick, [te finalTime], ye, op_stick);
        Lambda = zeros(length(t),2);
        q_dd = zeros(length(t),4);
        for i = 1:length(t)
            [q_dd(i,:),Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
        end
        
        %Record results
        finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))]; %transition state will be considered to belong to the next state
        finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))]; %NaNs added for plotting disconinuities
        finalTimes = [finalTimes;t(1:end-1);NaN];
        finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];

        %Record and update state
        stickInds = [stickInds; true(length(t),1)];
        sgn_slip = sign(1.5-ie);
    else
        %Solve corresponding ODE

        [t,X,te,ye,ie] = ode45(@sys_slip, [te finalTime], ye, op_slip);
        Lambda = zeros(length(t),2);
        q_dd = zeros(length(t),4);
        for i = 1:length(t)
            [q_dd(i,:),Lambda(i,2)] = dyn_sol_slip(X(i,1:4)',X(i,5:8)',t(i));
            Lambda(i,1) = -sgn_slip*mu*Lambda(i,2);
        end

        %Record results
        finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))];
        finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))];
        finalTimes = [finalTimes;t(1:end-1);NaN];
        finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];
        
        %Record and update state
        stickInds = [stickInds; false(length(t),1)];
        if ie-1
            lookingForWbreak = 0;
            break;
        end
        sgn_slip = -sgn_slip;
    end
end
finalX(end,:) = X(end,:);
finalLambda(end,:) = Lambda(end,:);
finalTimes(end,:) = t(end,:);
finalq_dd(end,:) = q_dd(end,:);

t = finalTimes;
x = finalX(:,1);
th = finalX(:,3);
ph = finalX(:,4);
lambdan = finalLambda(:,2);
lambdat = finalLambda(:,1);
vt = finalX(:,5)+R*finalX(:,7);

%% Q7 plot a)
figure;
h1 = plot(t(stickInds),th(stickInds)*180/pi,'LineWidth',2); hold on
h2 = plot(t(stickInds),ph(stickInds)*180/pi,'LineWidth',2);
plot(t(~stickInds),th(~stickInds)*180/pi,'--','Color',h1.Color)
plot(t(~stickInds),ph(~stickInds)*180/pi,'--','Color',h2.Color)
set(gcf,'color','w');
title('Angles vs. Time; $\frac{3\omega_{slip}+\omega_{break}}{4}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Angle[$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('$\theta_{stick}$','$\phi_{stick}$','$\theta_{slip}$','$\phi_{slip}$','Interpreter','latex','fontsize',20,'location','ne')
grid on;
saveas(gcf, 'q7a.png');

%% Q7 plot b)
figure;
h1 = plot(t(stickInds),vt(stickInds),'LineWidth',2); hold on
plot(t(~stickInds),vt(~stickInds),'--','Color',h1.Color)
set(gcf,'color','w');
title('Slip Velocity vs. Time; $\frac{3\omega_{slip}+\omega_{break}}{4}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$v_t [\frac{m}{s}$]', 'Interpreter', 'latex', 'fontsize', 20);
legend('stick','slip','Interpreter','latex','fontsize',20,'location','ne')
grid on;
saveas(gcf, 'q7b.png');

%% Q7 plot c)
figure;
h1 = plot(t(stickInds),lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Normal Force vs. Time; $\frac{3\omega_{slip}+\omega_{break}}{4}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_{n}$ [N]','Interpreter', 'latex', 'fontsize', 20);
legend('stick','slip','Interpreter','latex','fontsize',20,'location','ne')
grid on;
saveas(gcf, 'q7c.png');

%% Q7 plot d)
figure;
h1 = plot(t(stickInds),lambdat(stickInds)./lambdan(stickInds),'LineWidth',2); hold on;
plot(t(~stickInds),lambdat(~stickInds)./lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);
set(gcf,'color','w');
title('Force ratio vs. Time; $\frac{3\omega_{slip}+\omega_{break}}{4}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_{t}}{\lambda_{n}}$', 'Interpreter', 'latex', 'fontsize', 20);
grid on;
yline(mu)
yline(-mu)
ylim([-0.1,0.1])
legend('stick','slip','$\mu$ limits','Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'q7d.png');
lambdan(stickInds)

%% Q7 plot e)
figure;
h1 = plot(t(stickInds),x(stickInds),'LineWidth',2); hold on
plot(t(~stickInds),x(~stickInds),'--','Color',h1.Color)
set(gcf,'color','w');
title('Horizontal Displacement vs. Time; $\frac{3\omega_{slip}+\omega_{break}}{4}$','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$x$ [m]', 'Interpreter', 'latex', 'fontsize', 20);
legend('stick','slip','Interpreter','latex','fontsize',20,'location','se')
grid on;
saveas(gcf, 'q7e.png');

%% Q8

muCrit = 0.01;
rez = 0.01;
rez2 = (1e-1)*pi/6;
alpha1 = 1; alpha2 = 1;
muCritNotReached = 1;
while muCritNotReached
    muCrit = muCrit + rez;
    theta = -pi/6:rez2:pi/6;
    phi = -pi:rez2:pi;
    [T,P] = meshgrid(theta,phi);
    alpha1 = zeros(size(T));
    alpha2 = zeros(size(T));
    for i = 1:numel(T)
            %We only need M, and M depends only on theta and phi
            [M,~,~,~,wn,wt] = dynamics_mat([NaN NaN T(i) P(i)],NaN(4,1));
            alpha1(i) = (wn/M)*(wn-muCrit*wt)';
            alpha2(2) = (wn/M)*(wn+muCrit*wt)';

            (wn/M)*(wn-muCrit*wt)'
            (wn/M)*(wn+muCrit*wt)'

            if alpha1(i)<0 || alpha2(i)<0
                muCritNotReached = 0;
            end
    end
end

% %Show smoothness of surfaces
% figure;
% surf(T,P,alpha1)
% figure;
% surf(T,P,alpha2)

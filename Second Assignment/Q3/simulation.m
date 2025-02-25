% HW2 Q3 Ruigang Chen &  Ben Sarfati
clear all; close all; clc

%% solve ode45
X = [0; 0; 0; 0; 0; 0; 0; 0];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
tspan = [0 100]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);         
[t, X] = ode45(@state_eq, t_eval, X, options); 

Lambda = zeros(length(t_eval),2);
slippageVelocity = Lambda;
X_d = zeros(size(X(:,1:4)));
for i = 1:tspan(2)/dt+1
    [X_d(i,1:4), Lambda(i,1:2)] = dyn_sol(X(i,1:4)',X(i,5:8)',t(i));
    
    %Slippage velocities
    [~,~,W,~] = dynamics_mat(X(i,1:4)',X(i,5:8)');
    slippageVelocity(i,1:2) = (W*(X(i,5:8)'))';
end

%% Label data

m = 30;  
d = 0.25;

x = X(:,1);
y = X(:,2);
th = X(:,3);
ph = X(:,4);
x_d = X(:,5);
y_d = X(:,6);
th_d = X(:,7);
ph_d = X(:,8);
x_dd = X_d(:,1);
y_dd = X_d(:,2);
th_dd = X_d(:,3);
ph_dd = X_d(:,4);
lambda1 = Lambda(:,1);
lambda2 = Lambda(:,2);

rP_d = [x_d y_d];
rC_dd = [x_dd-th_d.^2*d.*cos(th)-th_dd*d.*sin(th) y_dd-th_d.^2*d.*sin(th)+th_dd*d.*cos(th)];
e1tag = [cos(th) sin(th)];
e1tagtag = [cos(th+ph) sin(th+ph)];

%% plots Q3 a)
 
close all; 

figure;
plot(t,ph*180/pi,'LineWidth',2);
set(gcf,'color','w');
title('Steering Angle vs. Time','fontsize',20)
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\phi $ [$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
saveas(gcf, 'q3a.png');

%% plots Q3 b)
figure;
plot(t,th*180/pi,'LineWidth',2);
set(gcf,'color','w');
title('Body Orientation Angle vs. Time','fontsize',20)
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\theta $ [$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
saveas(gcf, 'q3b.png');

%% plots Q3 c)
figure;
plot(t,dot(rP_d,e1tag,2),'LineWidth',2);
set(gcf,'color','w');
title('Velocity of P parallel to back link vs. Time','fontsize',20)
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity [$\frac{m}{s}$]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
saveas(gcf, 'q3c.png');

%% plots Q3 d)
figure;
plot(x,y,'LineWidth',2);
set(gcf,'color','w');
title('Trajectory of P','fontsize',20)
xlabel('$\mathbf{r_P}\cdot\mathbf{e_1}$ [m]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\mathbf{r_P}\cdot\mathbf{e_2}$ [m]', 'Interpreter', 'latex', 'fontsize', 20);
grid on;
axis equal;
saveas(gcf, 'q3d.png');

%% plots Q3 e)
figure;
plot(t,lambda1,'LineWidth',2); hold on;
plot(t,lambda2,'LineWidth',2);
set(gcf,'color','w');
title('Constraint Forces vs Time','fontsize',20)
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Constraint Force [N]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
lgd = legend('$\lambda_1$','$\lambda_2$','Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
saveas(gcf, 'q3e.png');

%% plots Q3 f)
figure;
plot(t(1:100:100001),m*rC_dd(1:100:100001,1),'o','LineWidth',2); hold on;
plot(t,-lambda1.*e1tag(:,2)-lambda2.*e1tagtag(:,2),'LineWidth',1);
set(gcf,'color','w');
title('Ground Reaction Force in $\mathbf{e_1}$ vs Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$GRF_x$ [N]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
lgd = legend('$m\mathbf{\ddot r_C}\cdot\mathbf{e_1}$',...
    '$(\lambda_1\mathbf{e^,_2}+\lambda_2\mathbf{e^{,,}_2}+F_d)\cdot\mathbf{e_1}$','Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
saveas(gcf, 'q3f.png');

%% plots Q3 g)
figure;
plot(t,slippageVelocity(:,1),'LineWidth',2); hold on;
plot(t,slippageVelocity(:,2),'LineWidth',2);
set(gcf,'color','w');
title('Slippage velocity of back and front wheels vs Time','fontsize',20)
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity [$\frac{m}{s}$]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
lgd = legend('Back wheels','Front wheel','Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
saveas(gcf, 'q3g.png');
% HW2 Ruigang Chen Ben Sarfati
clear all; close all; clc

%% solve ode45
X = [0; 0; 0; 0; 0; 0; 0; 0];  %  [x; y; theta; phi; dx; dy; dtheta; dphi];
tspan = [0 100]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);         
[t, X] = ode45(@state_eq, t_eval, X, options); 


%% phi data
phi_data = [];
phi_d_data = [];
phi_dd_data = [];
for i = 1:length(t)
    [phi, phi_d, phi_dd] = angles_input(t(i)); 
    phi_data = [phi_data, phi];
    phi_d_data = [phi_d_data, phi_d];
    phi_dd_data = [phi_dd_data, phi_dd];
end
phi1 = phi_data(1,:)';
phi2 = phi_data(2,:)';
phi1_d = phi_d_data(1,:)';
phi2_d = phi_d_data(2,:)';
phi2_dd = phi_dd_data(2,:)';
%% data
x = X(:, 1);
y = X(:, 2);
theta = X(:, 3);
x_d = X(:, 4);
y_d = X(:, 5);
theta_d = X(:, 6);

q = [x y theta phi1 phi2];
q_d = [x_d y_d theta_d phi1_d phi2_d];

qb = [x y theta];
qb_d = [x_d y_d theta_d];


%% COM
r_com_data = [];
r_d_com_data = [];
for i = 1:length(t)
    [r_com, r_d_com] = center_of_mass(q(i,:), q_d(i,:));
    r_com_data = [r_com_data, r_com];
    r_d_com_data = [r_d_com_data, r_d_com];
end


%% tau and qb_dd

tau_data = [];
qb_dd_data = [];
for i = 1:length(t)
    [qb_dd, tau] = dyn_sol(q(i,:)',q_d(i,:)',t(i));
    tau_data = [tau_data, tau];
    qb_dd_data = [qb_dd_data, qb_dd];
end
theta_dd = qb_dd_data(3,:)';


m = 1;
l = 0.1;
w = 1;
g = 10;
%% Plot the angle of middle link 

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, theta * 180 / pi, 'LineWidth', 2);
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\theta$ [$^{\circ}$]', 'Interpreter', 'latex', 'fontsize', 20);
ylim([-230 50])
grid on;
title('Angle of the middle link', 'fontsize', 20);
saveas(gcf, 'theta.png');

%%  x(t)/l and xc(t)/l 

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, x/l, 'LineWidth', 2, 'DisplayName', '$\frac{x(t)}{l}$');
plot(t, r_com_data(1,:)/l, '--', 'LineWidth', 2, 'DisplayName', '$\frac{x_c(t)}{l}$');
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{x}{l}$', 'Interpreter', 'latex', 'fontsize', 30);
ylim([-0.8 0.3])
lgd = legend('Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
title('Normalized horizontal position between the middle link and COM', 'fontsize', 20);
saveas(gcf, 'xt_xct.png');


%% x_dot(t)/lw and xc_dot(t)/lw 

f0 = [];
f1 = [];
f2 = [];
Hc = [];
for i = 1:length(t)
    f0 = [f0, (m*l^2)/3*(13+2*cos(phi1(i)-phi2(i))+6*(cos(phi1(i))+cos(phi2(i))))];
    f1 = [f1, (m*l^2)/3*(3+cos(phi1(i)-phi2(i))+3*cos(phi1(i)))];
    f2 = [f2, (m*l^2)/3*(3+cos(phi1(i)-phi2(i))+3*cos(phi2(i)))];
    Hc = [Hc, f0(i)*theta_d(i)+f1(i)*phi1_d(i)+f2(i)*phi2_d(i)];
end

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, x_d/(l*w), '-', 'LineWidth', 2, 'DisplayName', '$\frac{\dot{x}(t)}{l \omega}$');
plot(t, r_d_com_data(1,:)'/(l*w), '-', 'LineWidth', 4, 'DisplayName', '$\frac{\dot{x}_c(t)}{l \omega}$');
plot(t, Hc(1,:)/(m*l^2*w), '-.', 'LineWidth', 2, 'DisplayName', '$\frac{H_{c}(t)}{ml^2 \omega}$');
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\dot x}{l \cdot \omega}$', 'Interpreter', 'latex', 'fontsize', 30);
ylim([-0.8 1.2])
lgd = legend('Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
title('Normalized horizontal velocity between the middle link and COM', 'fontsize', 20);
saveas(gcf, 'xd_xcd.png');


%% y_dot(t)/lw and yc_dot(t)/lw

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, y_d/(l*w), 'LineWidth', 2, 'DisplayName', '$\frac{\dot{y}(t)}{l \omega}$');
plot(t, r_d_com_data(2,:)'/(l*w), '--', 'LineWidth', 2, 'DisplayName', '$\frac{\dot{y}_c(t)}{l \omega}$');
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\dot y}{l \cdot \omega}$', 'Interpreter', 'latex', 'fontsize', 30);
lgd = legend('Location','NorthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
grid on;
title('Normalized vertical velocity between the middle link and COM', 'fontsize', 20);
saveas(gcf, 'yd_ycd.png');

%%  thtea/w
cal_dtheta = [];
for i = 1:length(t)
    cal_dtheta = [cal_dtheta, (0 - f1(i)*phi1_d(i) - f2(i)*phi2_d(i))/f0(i)];
end

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, theta_d/(w), 'LineWidth', 2, 'DisplayName', '$\frac{\dot{\theta}(t)}{\omega^2}$')
plot(t, cal_dtheta/(w), '-.', 'LineWidth', 2, 'DisplayName', '$\frac{\dot \theta_{d.}(t)}{\omega^2}$')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\dot \theta}{\omega}$', 'Interpreter', 'latex', 'fontsize', 30);
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20; 
ylim([-0.8 0.8])
grid on;
title('Normalized angular velocity', 'fontsize', 20);
saveas(gcf, 'theta_w.png');


%% tau

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, tau_data(1,:), 'LineWidth', 2, 'DisplayName', '$\tau_1$')
plot(t, tau_data(2,:), 'LineWidth', 2, 'DisplayName', '$\tau_2$')
xlabel('Time [s]', 'Interpreter', 'latex','fontsize', 20);
ylabel('$\tau$ [N$\cdot$m]', 'Interpreter', 'latex','fontsize', 20);
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20; 
ylim([-0.03 0.03])
grid on;
title('Actuation torques at the joints \tau_1 and  \tau_2','fontsize', 20);
saveas(gcf, 'tau.png');

%%  theta_dd/w
I = 1/3 * m * l^2;
x_dd = qb_dd_data(1,:)';
y_dd = qb_dd_data(2,:)';
tau2 = tau_data(2,:);
cal_ddtheta = [];
for i = 1:length(t)
    cal_ddtheta = [cal_ddtheta, (tau2(i) - m*g*l*cos(theta(i) + phi2(i)) - (4*l^2*m/3)*phi2_dd(i) - l^2*m*sin(phi2(i))*theta_d(i)^2 - l*m*cos(phi2(i) + theta(i))*y_dd(i) + l*m*sin(phi2(i) + theta(i))*x_dd(i))/((4*l^2*m/3) + l^2*m*cos(phi2(i)))];
end

figure('Position', [250, 250, 800, 400]);
hold on;
plot(t, theta_dd/w^2, 'LineWidth', 2, 'DisplayName', '$\frac{\ddot{\theta}(t)}{\omega^2}$')
plot(t, cal_ddtheta/w^2, '-.','LineWidth', 2, 'DisplayName', '$\frac{\ddot\theta_{d.}(t)}{\omega^2}$')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\ddot \theta}{\omega^2}$', 'Interpreter', 'latex', 'fontsize', 30);
ylim([-1.5 1.5])
lgd = legend;  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20; 
grid on;
title('Normalized angular acceleration', 'fontsize', 20);
saveas(gcf, 'theta_dd_w.png');


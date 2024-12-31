% HW2 Q6 Ruigang Chen & Ben Sarfati
clear all; close all; clc

%% solve ode45
X0 = [0; 0; 0; 0; 0; 0;];  %  [x; y; theta; dx; dy; dtheta;];
tspan = [0 60]; 
dt = 0.001; 
t_eval = tspan(1):dt:tspan(2);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

c_term = 0; % change c here
[t, X] = ode45(@(t, X) state_eq(t, X, c_term), t_eval, X0, options); 
c_term = 1; % change c here
[t_c, X_c] = ode45(@(t, X) state_eq(t, X, c_term), t_eval, X0, options);

%% data
phi_data = [];
phi_d_data = [];
for i = 1:length(t)
    [phi, phi_d, phi_dd] = angles_input(t(i)); 
    phi_data = [phi_data, phi];
    phi_d_data = [phi_d_data, phi_d];
end
phi = phi_data(1,:)';
phi_d = phi_d_data(1,:)';

x = X(:, 1);
y = X(:, 2);
theta = X(:, 3);
x_d = X(:, 4);
y_d = X(:, 5);
theta_d = X(:, 6);
q = [x y theta phi];
q_d = [x_d y_d theta_d phi_d];

x_c = X_c(:, 1);
y_c = X_c(:, 2);
theta_c = X_c(:, 3);
x_d_c = X_c(:, 4);
y_d_c = X_c(:, 5);
theta_d_c = X_c(:, 6);
q_c = [x_c y_c theta_c phi];
q_d_c = [x_d_c y_d_c theta_d_c phi_d];

tau_data = [];
qp_dd_data = [];
lambda_data = [];
for i = 1:length(t)
    [qp_dd, tau, lambda] = dyn_sol(q(i,:)',q_d(i,:)',t(i),0);
    tau_data = [tau_data, tau];
    qp_dd_data = [qp_dd_data, qp_dd];
    lambda_data = [lambda_data, lambda];
end
x_dd = qp_dd_data(1,:)';
y_dd = qp_dd_data(2,:)';
theta_dd = qp_dd_data(3,:)';
lambda1 = lambda_data(1,:)';
lambda2 = lambda_data(2,:)';


tau_data_c = [];
qp_dd_data_c = [];
lambda_data_c = [];
for i = 1:length(t_c)
    [qp_dd_c, tau_c, lambda_c] = dyn_sol(q_c(i,:)',q_d_c(i,:)',t_c(i),1);
    tau_data_c = [tau_data_c, tau_c];
    qp_dd_data_c = [qp_dd_data_c, qp_dd_c];
    lambda_data_c = [lambda_data_c, lambda_c];
end
x_dd_c = qp_dd_data_c(1,:)';
y_dd_c = qp_dd_data_c(2,:)';
theta_dd_c = qp_dd_data_c(3,:)';
lambda1_c = lambda_data_c(1,:)';
lambda2_c = lambda_data_c(2,:)';



d = 0.15;
l = 0.6;
b = 0.2;
w = 0.25;

rP_d = [x_d y_d];
rC_dd = [x_dd-theta_d.^2*d.*cos(theta)-theta_dd*d.*sin(theta) y_dd-theta_d.^2*d.*sin(theta)+theta_dd*d.*cos(theta)];
e1tag = [cos(theta) sin(theta)];
e1tagtag = [cos(theta+phi) sin(theta+phi)];

rP_d_c = [x_d_c y_d_c];
rC_dd_c = [x_dd_c - theta_d_c.^2 .* d .* cos(theta_c) - theta_dd_c .* d .* sin(theta_c) 
           y_dd_c - theta_d_c.^2 .* d .* sin(theta_c) + theta_dd_c .* d .* cos(theta_c)];
e1tag_c = [cos(theta_c) sin(theta_c)];
e1tagtag_c = [cos(theta_c+phi) sin(theta_c+phi)];


rF_d_c = [x_d_c - theta_d_c .* l .* sin(theta_c) - b .* sin(theta_c + phi) .* (theta_d_c + phi_d) y_d_c + theta_d_c .* l .* cos(theta_c) + b .* cos(theta_c + phi) .* (theta_d_c + phi_d)];
rR_d_c = [x_d_c + 0.5 .* w .* theta_d_c .* cos(theta_c) y_d_c + 0.5 .* w .* theta_d_c .* sin(theta_c)];
rL_d_c = [x_d_c - 0.5 .* w .* theta_d_c .* cos(theta_c) y_d_c - 0.5 .* w .* theta_d_c .* sin(theta_c)];

F_d_c = (rF_d_c+rR_d_c+rL_d_c).*e1tag_c;

%% plots Q6 c)
figure;
plot(t,dot(rP_d,e1tag,2),'b','LineWidth',2);
hold on;
plot(t_c,dot(rP_d_c,e1tag_c,2),'r','LineWidth',2);
set(gcf,'color','w');
title('Velocity of P parallel to back link vs. Time','fontsize',20)
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Velocity [$\frac{m}{s}$]', 'Interpreter', 'latex', 'fontsize', 30);
grid on;
lgd = legend('c=0','c=1','Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20; 
saveas(gcf, 'q6c.png');
%% plots Q6 d)
figure;
plot(x,y,'b','LineWidth',2);
hold on;
plot(x_c,y_c,'r','LineWidth',2);
set(gcf,'color','w');
title('Trajectory of P','fontsize',20)
xlabel('$\mathbf{r_P}\cdot\mathbf{e_1}$ [m]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\mathbf{r_P}\cdot\mathbf{e_2}$ [m]', 'Interpreter', 'latex', 'fontsize', 20);
grid on;
lgd = legend('c=0','c=1','Location','SouthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20; 
axis equal;
saveas(gcf, 'q6d.png');

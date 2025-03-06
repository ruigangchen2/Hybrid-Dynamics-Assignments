clear; clc;

%% find periodic solution

Z0 = [-0.149, 0.733, -0.501].';

for i = 1:100
    [Z0_periodic, ~, ~, ~, ~] = fsolve(@(Z)(Poincare_map(Z) - Z), Z0);
    Z0 = Z0_periodic;
end

disp("th1:")
disp(Z0(1))
disp("th2:")
disp(Z0(1))
disp("th1_d:")
disp(Z0(2))
disp("th2_d:")
disp(Z0(3))

X0 = [0 0 Z0(1) Z0(1) 0 0 Z0(2) Z0(3)]

%% ODE
[m, mh, Ic, l, g, a, mu] = model_params;

impact_status = 0; % 0: stick  1: impact

t_start = 0;
t_stop = 10;
t_current = 0;
options = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t,X)events_stick(t,X));

X_data = [];
t_data = [];

while abs(t_stop-t_current) > 0.0001

    if impact_status == 0
        [t,X,~,~,ie] = ode45(@(t,X)sys_stick(t,X),[t_start,t_stop],X0,options);
        X_data = [X_data; X];
        t_data = [t_data; t];
        t_start = t(end);

        X0 = [X(end,1),X(end,2),X(end,3),X(end,4),X(end,5),X(end,6),X(end,7),X(end,8)]';
        
        if ~isempty(ie)
            if ie(end) == 3 % hip  collision
                break;
            elseif ie(end) == 4 % foot's collision
                impact_status = 1; 
                Xold = X0;
            end
        else
            break; 
        end

    end

    if impact_status == 1
        X0 = impact_law(Xold);
    end

    impact_status = 0;
    t_current = t_data(end);

end
    
%%
x = X_data(:, 1);
y = X_data(:, 2);
th1 = X_data(:, 3);
th2 = X_data(:, 4);
x_d = X_data(:, 5);
y_d = X_data(:, 6);
th1_d = X_data(:, 7);
th2_d = X_data(:, 8);

%% Plot a)

figure;
plot(t_data(1), th1(1), 'r-', 'LineWidth', 2)
hold on
plot(t_data(1), th2(1), 'b-', 'LineWidth', 2)

plot(t_data(1), th1(1), 'r--', 'LineWidth', 2) 
plot(t_data(1), th2(1), 'b--', 'LineWidth', 2)

% Scuffing
indx_scuffing = find(abs(th1-th2) < 0.001 & th1 < 0.1 & th1 > -0.1);
t_scuff = t_data(indx_scuffing);
th_scuff = th1(indx_scuffing);
plot(t_scuff, th_scuff, 'k*', 'LineWidth', 2)

s = 1;
th1_prev = th1(1);
th2_prev = th2(1);
for i = 2:length(t_data) 
    th1_next = th1(i);
    th2_next = th2(i);
    if (abs(th1_next - th1_prev) > 0.01) % choose to escape the data
        plot(t_data(s:i-1), th1(s:i-1), 'r-', 'LineWidth', 2)        
        plot(t_data(i-1:i), th1(i-1:i), 'r', 'LineWidth', 2)
        plot(t_data(s:i-1), th2(s:i-1), 'b-', 'LineWidth', 2)        
        plot(t_data(i-1:i), th2(i-1:i), 'b--', 'LineWidth', 2)
        s = i;
    end
    th1_prev = th1(i);
    th2_prev = th2(i);
end

plot(t_data(s:end), th1(s:end), 'r-', 'LineWidth', 2)
plot(t_data(s:end), th2(s:end), 'b-', 'LineWidth', 2)
plot(t_scuff, th_scuff, 'k*', 'LineWidth', 2)

title('Angles vs. Time inlcuding scuffing','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\theta}$ [rad]', 'Interpreter', 'latex', 'fontsize', 20);
ylim([-0.3 0.7]);
legend("${\theta_1(t)}$", "${\theta_2(t)}$", "${\theta_1(t)}$ Relabling","${\theta_2(t)}$ Relabling","Scuffing",'Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'a.png');

%% Plot b)
figure;
plot(th1(1), th1_d(1), 'r-', 'LineWidth', 2)
hold on
plot(th2(1), th2_d(1), 'b-', 'LineWidth', 2)
plot(th1(1), th1_d(1), 'r-.', 'LineWidth', 2)
plot(th2(1), th2_d(1), 'b-.', 'LineWidth', 2)
% Scuffing:
plot(th1(indx_scuffing), th1_d(indx_scuffing), 'ko', 'LineWidth', 3)

% Collision:
indx_collision = find(abs(th1-th2) < 0.001 & (th1 > 0.1 | th1 < -0.1));
plot(th1(indx_collision), th1_d(indx_collision), 'kx', 'LineWidth', 3)


s = 1;
th1_prev = th1(1);
th2_prev = th2(1);
for i = 2:length(t_data) 
    th1_next = th1(i);
    th2_next = th2(i);
    if (abs(th1_next - th1_prev) > 0.01)
        plot(th1(s:i-1), th1_d(s:i-1), 'r-', 'LineWidth', 2)        
        plot(th1(i-1:i), th1_d(i-1:i), 'r-.', 'LineWidth', 2)
        plot(th2(s:i-1), th2_d(s:i-1), 'b-', 'LineWidth', 2)        
        plot(th2(i-1:i), th2_d(i-1:i), 'b-.', 'LineWidth', 2)
        s = i;
    end
    th1_prev = th1(i);
    th2_prev = th2(i);
end
plot(th1(s:end), th1_d(s:end), 'r-', 'LineWidth', 2)
plot(th2(s:end), th2_d(s:end), 'b-', 'LineWidth', 2)

% Scuffing:
plot(th1(indx_scuffing), th1_d(indx_scuffing), 'ko', 'LineWidth', 3)
plot(th2(indx_scuffing), th2_d(indx_scuffing), 'ko', 'LineWidth', 3)

% Collision:
plot(th1(indx_collision), th1_d(indx_collision), 'kx', 'LineWidth', 3)
plot(th2(indx_collision), th2_d(indx_collision), 'kx', 'LineWidth', 3)


title('Phase plane trajectories','fontsize',20,'Interpreter','latex')
xlabel('${\theta}$ [rad]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\dot{\theta}}$ [rad/s]', 'Interpreter', 'latex', 'fontsize', 20);
ylim([-2 2]);
xlim([-0.25 0.55]);
legend("${\theta_1(t)}$", "${\theta_2(t)}$", "Relabeling", "Relabeling", "Scuffing",  "Collision",'Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'b.png');


%% Plot c)
lambda_t_data = [];
lambda_n_data = [];
for i = 1:length(t_data)
    [~, lambda] = dyn_sol_stick(t_data(i), X_data(i, :));
    lambda_t_data = [lambda_t_data, lambda(1)];
    lambda_n_data = [lambda_n_data, lambda(2)];
end

figure;
plot(t_data, lambda_n_data, "r-", 'LineWidth', 2); 
hold on;
plot([t_data(1) t_data(end)], [0 0], 'k-.', 'LineWidth', 2);
title('$\lambda_n$ vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\lambda_n$ [N]', 'Interpreter', 'latex', 'fontsize', 20);
legend("$\lambda_n$",'Interpreter','latex','fontsize',20,'location','ne')
ylim([-10 140]);
saveas(gcf, 'c.png');

%% Plot d)

figure;
plot(t_data, lambda_t_data./lambda_n_data, '-', 'LineWidth', 2)
hold on

LAMBDA_t = [];
LAMBDA_n = [];


for i = 1:length(indx_collision)
    
    [M, ~, ~, ~, ~, Wtilde] = dynamics_mat(X_data(indx_collision(i),:));
    A = (Wtilde/M)*Wtilde';
    vp_minus = Wtilde*X_data(indx_collision(i),5:end).';
    LAMBDA = -A\vp_minus;
    LAMBDA_t = [LAMBDA_t, LAMBDA(1)];
    LAMBDA_n = [LAMBDA_n, LAMBDA(2)];
end

LAMBDA_impact_index = find(LAMBDA_t./LAMBDA_n < 0);

plot(t_data(indx_collision(LAMBDA_impact_index)), LAMBDA_t(LAMBDA_impact_index)./LAMBDA_n(LAMBDA_impact_index), 'kx', 'LineWidth', 3)
plot([t_data(1) t_data(end)], mu*[1 1], 'k-.', 'LineWidth', 2)
plot([t_data(1) t_data(end)], -mu*[1 1], 'k-.', 'LineWidth', 2)


title('Force ratio vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Force ratio', 'Interpreter', 'latex', 'fontsize', 20);
ylim([-1.2 1.2]);
legend("${\lambda_t/\lambda_n}$", "${\Lambda_t/\Lambda_n}$", "$\pm \mu$", 'Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'd.png');

%% Q6

min_lam_ratio = max(abs(lambda_t_data ./ lambda_n_data));
min_Lam_ratio = max(abs(LAMBDA_t./LAMBDA_n));
max_abs_ratio = max(min_lam_ratio, min_Lam_ratio);
fprintf('min_mu：%.4f\n', max_abs_ratio);


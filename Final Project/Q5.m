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
t_stop = 5;
t_current = 0;
options = odeset('reltol',1e-8,'abstol',1e-8,'Events',@(t,X)events_stick(t,X));

X_data = [];
t_data = [];

while abs(t_stop-t_current) > 10^-3

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
x_dot = X_data(:, 5);
y_dot = X_data(:, 6);
th1_dot = X_data(:, 7);
th2_dot = X_data(:, 8);

r_p = [X_data(:,1)+2*l*sin(X_data(:,3))+2*l*sin(X_data(:,4)) X_data(:,2)+2*l*cos(X_data(:,3))-2*l*cos(X_data(:,4))];

%% Plot a)

figure;
plot(t_data(1), th1(1), 'r-', 'LineWidth', 2)
hold on
plot(t_data(1), th2(1), 'b-', 'LineWidth', 2)

% Scuffing
indx_rp_zero = find(abs(r_p(:,2)) < 0.0001);
i = 0:floor(length(indx_rp_zero)/3);
indx_scuff = indx_rp_zero(2+3*i);

t_scuff = t_data(indx_scuff);
th_scuff = th1(indx_scuff);
plot(t_scuff, th_scuff, 'k*', 'LineWidth', 2)

% Plotting
s = 1;
x_prev = th1(1);
for i = 1:length(t_data) 
    x_next = th1(i);
    if (abs(x_next - x_prev) > 0.1)
        plot(t_data(s:i-1), th1(s:i-1), 'r-', 'LineWidth', 2)        
        plot(t_data(i-1:i), th1(i-1:i), 'r-.', 'LineWidth', 1)
        s = i;
    end
    x_prev = th1(i);
end

plot(t_data(s:end), th1(s:end), 'r-', 'LineWidth', 2)

s = 1;
x_prev = th2(1);
for i = 2:length(t_data) 
    x_next = th2(i);
    if (abs(x_next - x_prev) > 0.1)
        plot(t_data(s:i-1), th2(s:i-1), 'b-', 'Color', '#0070C0', 'LineWidth', 2)        
        plot(t_data(i-1:i), th2(i-1:i), 'b-.', 'Color', '#0070C0', 'LineWidth', 2)
        s = i;
    end
    x_prev = th2(i);
end
plot(t_data(s:end), th2(s:end), 'b-', 'Color', '#0070C0', 'LineWidth', 2)


title('Angles vs. Time inlcuding scuffing','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\theta}$ [rad]', 'Interpreter', 'latex', 'fontsize', 20);
ylim([-0.3 0.5]);
legend("${\theta_1(t)}$", "${\theta_2(t)}$", "Scuffing",'Interpreter','latex','fontsize',20,'location','ne')
saveas(gcf, 'a.png');
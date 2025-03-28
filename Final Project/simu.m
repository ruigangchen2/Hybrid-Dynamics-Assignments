%Final project --- Ruigang Chen, Ben Sarfati
clear; clc;

global mu ;
mu = 1;

%% find periodic solution

Z0 = [-0.149, 0.733, -0.501].';

%Converges after 1 iteration
numIters = 1;
Z0_periodic = zeros(3,numIters);
for i = 1:numIters
    [Z0_periodic(:,i), ~, ~, ~, ~] = fsolve(@(Z)(Poincare_map(Z) - Z), Z0);
    Z0 = Z0_periodic(:,i);
end

% figure; 
% hold on;
% plot(Z0_periodic,'markerFaceColor','k');

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
[m, mh, Ic, l, g, a] = model_params;

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

title('Angles vs. Time including scuffing','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\theta}$ [rad]', 'Interpreter', 'latex', 'fontsize', 20);
ylim([-0.3 0.7]);
legend("${\theta_1(t)}$", "${\theta_2(t)}$", "${\theta_1(t)}$ Relabeling","${\theta_2(t)}$ Relabeling","Scuffing",'Interpreter','latex','fontsize',20,'location','ne')
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
% min_Lam_ratio = max(abs(LAMBDA_t./LAMBDA_n));
% max_abs_ratio = max(min_lam_ratio, min_Lam_ratio);
fprintf('min_mu：%.4f\n', min_lam_ratio);

%% Q11

%% find periodic solution

mu = 1.5*min_lam_ratio;
Z0slip = [-0.149, 0.733, -0.501, 0].';

numIters = 1;
Z0slip_periodic = zeros(4,numIters);
for i = 1:numIters
    [Z0slip_periodic(:,i), ~, ~, ~, ~] = fsolve(@(Z)(Poincare_map2(Z) - Z), Z0slip);
    Z0slip = Z0slip_periodic(:,i);
end

% figure; 
% hold on;
% plot(Z0_periodic,'markerFaceColor','k');

disp("th1:")
disp(Z0slip(1))
disp("th2:")
disp(Z0slip(1))
disp("th1_d:")
disp(Z0slip(2))
disp("th2_d:")
disp(Z0slip(3))
disp("x_d:")
disp(Z0slip(4))

X0 = [0 0 Z0slip(1) Z0slip(1) Z0slip(4) 0 Z0slip(2) Z0slip(3)]

%% ODE

global sgn_slip
sgn_slip = sign(X0(4));

finalX = [];
finalTimes = [];
finalLambda = [];
stickInds = logical([]);
relabelInds = logical([]);

Xe = X0;
te = 0;
finalTime = 10;
op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         
op_slip = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_slip);         
iter = 1;

while te<finalTime
    %Check slipping or sticking
    [value, ~, direction] = events_stick(te,Xe);
    if value(1:3).*direction(1:3) < 0
        [t,X,te,Xe,ie] = ode45(@sys_stick,[te, finalTime], Xe, op_stick);

        Lambda = zeros(length(t),2);
        for i = 1:length(t)
            [~,Lambda(i,1:2)] = dyn_sol_stick(t(i),X(i,:)');
        end            

        %Record results 
        finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))];
        finalTimes = [finalTimes;t(1:end-1);NaN];
        finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))];
        stickInds = [stickInds; true(length(t),1)];
        relabelInds = [relabelInds; false(length(t),1)];

        %Check for end
        if ~(t(end)<finalTime)
            break;
        end
        te = te(end);
        Xe = Xe(end,:);
        ie = ie(end);
        
        %Update state
        switch ie
            case 1
                sgn_slip = 1;
            case 2
                sgn_slip = -1;
            case 3
                error('failure; falling')
            case 4
                Xe = impact_law(Xe');
                relabelInds(end-1) = true;
        end
    else
        [t,X,te,Xe,ie] = ode45(@sys_slip, [te, finalTime], Xe, op_slip);

        Lambda = zeros(length(t),2);
        for i = 1:length(t)
            [~,Lambda(i,2)] = dyn_sol_slip(t(i),X(i,:)');
            Lambda(i,1) = -sgn_slip*mu*Lambda(i,2);
        end
        
        %Record results 
        finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))];
        finalTimes = [finalTimes;t(1:end-1);NaN];
        finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))];
        stickInds = [stickInds; false(length(t),1)];
        relabelInds = [relabelInds; false(length(t),1)];

        %Check for end
        if ~(t(end)<finalTime)
            break;
        end
        te = te(end);
        Xe = Xe(end,:);
        ie = ie(end);
        
        %Update state
        switch ie(end)
            case 1
                sgn_slip = -sgn_slip;
            case 2
                error('failure; stance foot separation')
            case 3
                error('failure; falling')
            case 4
                Xe = impact_law(Xe');
                relabelInds(end-1) = true;
        end
    end
    
    iter = iter +1;
end

finalX(end,:) = X(end,:);
finalLambda(end,:) = Lambda(end,:);
finalTimes(end,:) = t(end,:);


t = finalTimes;
th1 = finalX(:,3);
th2 = finalX(:,4);
th1_d = finalX(:,7);
th2_d = finalX(:,8);
lambdan = finalLambda(:,2);
lambdat = finalLambda(:,1);

indx_scuffing = find(abs(th1-th2) < 0.001 & th1 < 0.1 & th1 > -0.1);
t_scuff = t(indx_scuffing);
th_scuff = th1(indx_scuffing);
indx_collision = find(abs(th1-th2) < 0.001 & (th1 > 0.1 | th1 < -0.1));

%% plot (a) 

figure;
h1 = plot(t(stickInds),th1(stickInds),'LineWidth',2,'color','r'); hold on
h2 = plot(t(stickInds),th2(stickInds),'LineWidth',2,'color','b');
plot(t(~stickInds),th1(~stickInds),'--','Color',h1.Color)
plot(t(~stickInds),th2(~stickInds),'--','Color',h2.Color)

set(gcf,'color','w');
title('Angles vs. Time including scuffing','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\theta}$ [rad]', 'Interpreter', 'latex', 'fontsize', 20);

plot(t_scuff, th_scuff, 'k*', 'LineWidth', 2)

for i = find(relabelInds)'
    plot(t(i:2:i+2),th1(i:2:i+2),':','LineWidth',2,'Color',h1.Color)
    plot(t(i:2:i+2),th2(i:2:i+2),':','LineWidth',2,'Color',h2.Color)
end
legend("${\theta_1(t)}$ stick", "${\theta_2(t)}$ stick", "${\theta_1(t)}$ slip","${\theta_2(t)}$ slip","Scuffing","Relabel","Relabel",'Interpreter','latex','fontsize',20,'location','ne')

%% plot (b)

figure;
h1 = plot(th1(stickInds),th1_d(stickInds),'LineWidth',2,'color','r'); hold on
h2 = plot(th2(stickInds),th2_d(stickInds),'LineWidth',2,'color','b');
plot(th1(~stickInds),th1_d(~stickInds),'--','Color',h1.Color)
plot(th2(~stickInds),th2_d(~stickInds),'--','Color',h2.Color)

for i = find(relabelInds,1)'
    plot(th1(i:2:i+2),th1_d(i:2:i+2),':','LineWidth',2,'Color',h1.Color)
    plot(th2(i:2:i+2),th2_d(i:2:i+2),':','LineWidth',2,'Color',h2.Color)
end

% Scuffing:
plot(th1(indx_scuffing), th1_d(indx_scuffing), 'ko', 'LineWidth', 3)
plot(th2(indx_scuffing), th2_d(indx_scuffing), 'ko', 'LineWidth', 3,'HandleVisibility', 'off')

% Collision:
plot(th1(indx_collision), th1_d(indx_collision), 'kx', 'LineWidth', 3)
plot(th2(indx_collision), th2_d(indx_collision), 'kx', 'LineWidth', 3,'HandleVisibility', 'off')

legend("${\theta_1(t)}$ stick", "${\theta_2(t)}$ stick", "${\theta_1(t)}$ slip","${\theta_2(t)}$ slip","Scuffing","collision",'Interpreter','latex','fontsize',20,'location','se')
set(gcf,'color','w');
title('Phase Plane Trajectories','fontsize',20,'Interpreter','latex')
xlabel('${\theta}$ [rad/s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\dot\theta}$ [rad]', 'Interpreter', 'latex', 'fontsize', 20);

%% plot (c)

figure;
h1 = plot(t(stickInds),lambdan(stickInds),'LineWidth',2,'color','r'); hold on
plot(t(~stickInds),lambdan(~stickInds),'--','Color',h1.Color)

set(gcf,'color','w');
title('${\lambda_n}$ vs. time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('${\lambda}$ [N]', 'Interpreter', 'latex', 'fontsize', 20);

for i = find(relabelInds)'
    plot(t(i:2:i+2),lambdan(i:2:i+2),':','LineWidth',2,'Color',h1.Color)
end
plot([t_data(1) t_data(end)], [0 0], 'k-.', 'LineWidth', 2);

legend("$\lambda_n$",'relabel','Interpreter','latex','fontsize',20,'location','ne')

%% plot (d)

figure;
h1 = plot(t(stickInds),lambdat(stickInds)./lambdan(stickInds),'LineWidth',2,'color','r'); hold on
plot(t(~stickInds),lambdat(~stickInds)./lambdan(~stickInds),'--','LineWidth',2,'Color',h1.Color);

set(gcf,'color','w');
title('Force ratio vs. Time','fontsize',20,'Interpreter','latex')
xlabel('Time [s]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('$\frac{\lambda_{t}}{\lambda_{n}}$', 'Interpreter', 'latex', 'fontsize', 20);

for i = find(relabelInds)'
    plot(t(i:2:i+2),lambdat(i:2:i+2)./lambdan(i:2:i+2),':','LineWidth',2,'Color',h1.Color)
end
yline(mu)
yline(-mu)
legend("$\frac{\lambda_t}{\lambda_n}$",'relabel','Interpreter','latex','fontsize',20,'location','ne')
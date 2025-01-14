%HW3 Ben Sarfati Ruigang Chen
clear all; clc; close all;

%% Task 1

r1 = [-1; 1]; 
r2 = [2;2];
alpha = [pi/4; 3*pi/4];
n1 = [cos(alpha(1)) sin(alpha(1))];
t1 = [-sin(alpha(1)) cos(alpha(1))];
n2 = [cos(alpha(2)) sin(alpha(2))];
t2 = [-sin(alpha(2)) cos(alpha(2))];

mu = 0.4;
gamma = atan(mu);
[x_min,x_max] = fric_eq([r1 r2],alpha,mu*ones(2,1));
xc_min = x_min(end);
xc_max = x_max(end);

n1Plot = [linspace(r1(1),r1(1)+4*n1(1),10);linspace(r1(2),r1(2)+4*n1(2),10)];
n2Plot = [linspace(r2(1),r2(1)+4*n2(1),10);linspace(r2(2),r2(2)+4*n2(2),10)];
t1Plot = [linspace(r1(1)-0.5*t1(1),r1(1)+0.5*t1(1),10);linspace(r1(2)-0.5*t1(2),r1(2)+0.5*t1(2),10)];
t2Plot = [linspace(r2(1)-0.5*t2(1),r2(1)+0.5*t2(1),10);linspace(r2(2)-0.5*t2(2),r2(2)+0.5*t2(2),10)];
%cone plots
cone1 = [linspace(r2(1),r2(1)+2*(cos(alpha(2)-gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)-gamma)),10)];
cone1_ = [linspace(r2(1),r2(1)+2*(cos(alpha(2)+gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)+gamma)),10)];
cone2 = [linspace(r1(1),r1(1)+2*cos(alpha(1)-gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)-gamma),10)];
cone2_ = [linspace(r1(1),r1(1)+2*cos(alpha(1)+gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)+gamma),10)];

figure;
hold on;
plot(r1(1),r1(2),'k.','markersize',20); 
plot(n1Plot(1,:),n1Plot(2,:),'b:','linewidth',2);
plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
yl = ylim;
plot([xc_min,xc_min],[yl(1)-1,yl(2)+1],'g','linewidth',2);
patch([xc_min,xc_max,xc_max,xc_min],[yl(1)-1 yl(1)-1 yl(2)+1 yl(2)+1], 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.1);
plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
plot(r2(1),r2(2),'k.','markersize',20); 
plot(n2Plot(1,:),n2Plot(2,:),'b:','linewidth',2);
plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
plot(t2Plot(1,:),t2Plot(2,:),'k','linewidth',2);
plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
plot(cone1_(1,:),cone1_(2,:),'r','linewidth',2);
plot(cone2(1,:),cone2(2,:),'r','linewidth',2);
plot(cone2_(1,:),cone2_(2,:),'r','linewidth',2);

plot([xc_min,xc_min],[yl(1)-1,yl(2)+1],'g','linewidth',2);
plot([xc_max,xc_max],[yl(1)-1,yl(2)+1],'g','linewidth',2);
ylim([0,5])
axis equal;

xlabel('X [m]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Y [m]', 'Interpreter', 'latex', 'fontsize', 20);
set(gcf,'color','w');
title('Admissible COM Locations and Contact Visualization','fontsize',20);
grid on;
lgd = legend('contact point','contact normal','friction cone','boundary line','admissible COM location','Location','NorthWest');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
saveas(gcf, 'Task I (mu=0.4).png');

% change the mu
mu = 0.8;
gamma = atan(mu);
[x_min,x_max] = fric_eq([r1 r2],alpha,mu*ones(2,1));
xc_min = x_min(end);
xc_max = x_max(end);

n1Plot = [linspace(r1(1),r1(1)+4*n1(1),10);linspace(r1(2),r1(2)+4*n1(2),10)];
n2Plot = [linspace(r2(1),r2(1)+4*n2(1),10);linspace(r2(2),r2(2)+4*n2(2),10)];
t1Plot = [linspace(r1(1)-0.5*t1(1),r1(1)+0.5*t1(1),10);linspace(r1(2)-0.5*t1(2),r1(2)+0.5*t1(2),10)];
t2Plot = [linspace(r2(1)-0.5*t2(1),r2(1)+0.5*t2(1),10);linspace(r2(2)-0.5*t2(2),r2(2)+0.5*t2(2),10)];
%cone plots
cone1 = [linspace(r2(1),r2(1)+2*(cos(alpha(2)-gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)-gamma)),10)];
cone1_ = [linspace(r2(1),r2(1)+2*(cos(alpha(2)+gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)+gamma)),10)];
cone2 = [linspace(r1(1),r1(1)+2*cos(alpha(1)-gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)-gamma),10)];
cone2_ = [linspace(r1(1),r1(1)+2*cos(alpha(1)+gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)+gamma),10)];

figure;
hold on;
plot(r1(1),r1(2),'k.','markersize',20); 
plot(n1Plot(1,:),n1Plot(2,:),'b:','linewidth',2);
plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
yl = ylim;
plot([xc_min,xc_min],[yl(1),yl(2)],'g','linewidth',2);
patch([xc_min,xc_max,xc_max,xc_min],[yl(1)-1 yl(1)-1 yl(2)+1 yl(2)+1], 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.1);
plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
plot(r2(1),r2(2),'k.','markersize',20); 
plot(n2Plot(1,:),n2Plot(2,:),'b:','linewidth',2);
plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
plot(t2Plot(1,:),t2Plot(2,:),'k','linewidth',2);
plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
plot(cone1_(1,:),cone1_(2,:),'r','linewidth',2);
plot(cone2(1,:),cone2(2,:),'r','linewidth',2);
plot(cone2_(1,:),cone2_(2,:),'r','linewidth',2);
plot([xc_min,xc_min],[yl(1)-1,yl(2)+1],'g','linewidth',2);
plot([xc_max,xc_max],[yl(1)-1,yl(2)+1],'g','linewidth',2);
ylim([0,6])
axis equal;

xlabel('X [m]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Y [m]', 'Interpreter', 'latex', 'fontsize', 20);
set(gcf,'color','w');
title('Admissible COM Locations and Contact Visualization','fontsize',20);
grid on;
lgd = legend('contact point','contact normal','friction cone','boundary line','admissible COM location','Location','NorthEast');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 15;  
saveas(gcf, 'Task I (mu=0.8).png');

%% Task 2
% finished by handwriting.

%% Task 3
for mu = 0:0.1:3
    mu
    gamma = atan(mu);
    [x_min,x_max] = fric_eq([r1 r2],alpha,mu*ones(2,1));
    xc_min = x_min(end);
    xc_max = x_max(end);
    n1Plot = [linspace(r1(1),r1(1)+4*n1(1),10);linspace(r1(2),r1(2)+4*n1(2),10)];
    n2Plot = [linspace(r2(1),r2(1)+4*n2(1),10);linspace(r2(2),r2(2)+4*n2(2),10)];
    t1Plot = [linspace(r1(1)-0.5*t1(1),r1(1)+0.5*t1(1),10);linspace(r1(2)-0.5*t1(2),r1(2)+0.5*t1(2),10)];
    t2Plot = [linspace(r2(1)-0.5*t2(1),r2(1)+0.5*t2(1),10);linspace(r2(2)-0.5*t2(2),r2(2)+0.5*t2(2),10)];
    %cone plots
    cone1 = [linspace(r2(1),r2(1)+2*(cos(alpha(2)-gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)-gamma)),10)];
    cone1_ = [linspace(r2(1),r2(1)+2*(cos(alpha(2)+gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)+gamma)),10)];
    cone2 = [linspace(r1(1),r1(1)+2*cos(alpha(1)-gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)-gamma),10)];
    cone2_ = [linspace(r1(1),r1(1)+2*cos(alpha(1)+gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)+gamma),10)];
    
    figure;
    hold on;

    if isnan(x_max)
        if isnan(x_min)
            plot(r1(1),r1(2),'k.','markersize',20); 
            plot(n1Plot(1,:),n1Plot(2,:),'b:','linewidth',2);
            plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
            yl = ylim;
            patch([-5,5,5,-5],[yl(1)-9 yl(1)-9 yl(2)+9 yl(2)+9], 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.1);
            plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
            plot(r2(1),r2(2),'k.','markersize',20); 
            plot(n2Plot(1,:),n2Plot(2,:),'b:','linewidth',2);
            plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
            plot(t2Plot(1,:),t2Plot(2,:),'k','linewidth',2);
            plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
            plot(cone1_(1,:),cone1_(2,:),'r','linewidth',2);
            plot(cone2(1,:),cone2(2,:),'r','linewidth',2);
            plot(cone2_(1,:),cone2_(2,:),'r','linewidth',2);
            xticks([-4 -3 -2 -1 0 1 2 3 4])
            xticklabels({'-\infty','...','-2','-1','0','1','2','...','+\infty'})
            ylim([-0.5,5])
            xlim([-4.8,4.8])
            set(gcf,'color','w');
            title('Admissible COM Locations and Contact Visualization','fontsize',20);
            grid on;
            lgd = legend('contact point','contact normal','friction cone','admissible COM location','Location','NorthWest');  
            lgd.Interpreter = 'latex';  
            lgd.FontSize = 15;  
        else
            plot(r1(1),r1(2),'k.','markersize',20); 
            plot(n1Plot(1,:),n1Plot(2,:),'b:','linewidth',2);
            plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
            yl = ylim;
            plot([xc_min,xc_min],[yl(1)-10,yl(2)+10],'g','linewidth',2);
            patch([xc_min,9,9,xc_min],[yl(1)-10 yl(1)-10 yl(2)+10 yl(2)+10], 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.1);
            plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
            plot(r2(1),r2(2),'k.','markersize',20); 
            plot(n2Plot(1,:),n2Plot(2,:),'b:','linewidth',2);
            plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
            plot(t2Plot(1,:),t2Plot(2,:),'k','linewidth',2);
            plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
            plot(cone1_(1,:),cone1_(2,:),'r','linewidth',2);
            plot(cone2(1,:),cone2(2,:),'r','linewidth',2);
            plot(cone2_(1,:),cone2_(2,:),'r','linewidth',2);

            xticks([-1 0 1 2 3 4 5 6 7])
            xticklabels({'-1','0','1','2','3','4','5','...','\infty'})
            
            ylim([-1,7])
            xlim([-2,7])
            
            set(gcf,'color','w');
            title('Admissible COM Locations and Contact Visualization','fontsize',20);
            grid on;
            lgd = legend('contact point','contact normal','friction cone','boundary line','admissible COM location','Location','NorthWest');  
            lgd.Interpreter = 'latex';  
            lgd.FontSize = 15;  
        end
    else
        plot(r1(1),r1(2),'k.','markersize',20); 
        plot(n1Plot(1,:),n1Plot(2,:),'b:','linewidth',2);
        plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
        yl = ylim;
        plot([xc_min,xc_min],[yl(1)-10,yl(2)+10],'g','linewidth',2);
        patch([xc_min,xc_max,xc_max,xc_min],[yl(1)-10 yl(1)-10 yl(2)+10 yl(2)+10], 'black', 'FaceColor', 'cyan', 'FaceAlpha', 0.1);     
        plot([xc_max,xc_max],[yl(1)-10,yl(2)+10],'g','linewidth',2);
        
        plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
        plot(r2(1),r2(2),'k.','markersize',20); 
        plot(n2Plot(1,:),n2Plot(2,:),'b:','linewidth',2);
        plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
        plot(t2Plot(1,:),t2Plot(2,:),'k','linewidth',2);
        plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
        plot(cone1_(1,:),cone1_(2,:),'r','linewidth',2);
        plot(cone2(1,:),cone2(2,:),'r','linewidth',2);
        plot(cone2_(1,:),cone2_(2,:),'r','linewidth',2);
        
        set(gcf,'color','w');
        title('Admissible COM Locations and Contact Visualization','fontsize',20);
        grid on;
        if mu == 0
            lgd = legend('contact point','contact normal','friction cone','boundary line','Location','NorthWest');  
        else
            lgd = legend('contact point','contact normal','friction cone','boundary line','admissible COM location','Location','NorthEast');  
        end
        ylim([-1,6])
        xlim([-2,11])
        
        lgd.Interpreter = 'latex';  
        lgd.FontSize = 15;  
    end

    axis equal;
    xlabel('X [m]', 'Interpreter', 'latex', 'fontsize', 20);
    ylabel('Y [m]', 'Interpreter', 'latex', 'fontsize', 20);
    file_name = ['Task III (mu=' num2str(mu) ')' '.png'];%文件名称
    saveas(gcf, file_name);
end

%% Task 4
rc = [1;5];
mg = 1;
mu = 0.4;
gamma = atan(mu);
[f_min,f_max] = range_fx_eq([r1 r2 rc],alpha,mu*ones(2,1));
fx_min = f_min(end)
fx_max = f_max(end)

total_force = [linspace(rc(1),rc(1)+6*(sin(atan(fx_min/mg))),10);linspace(rc(2),rc(2)-6*(cos(atan(fx_min/mg))),10)];
total_force_ = [linspace(rc(1),rc(1)+6*(sin(atan(fx_max/mg))),10);linspace(rc(2),rc(2)-6*(cos(atan(fx_max/mg))),10)];

n1Plot = [linspace(r1(1),r1(1)+4*n1(1),10);linspace(r1(2),r1(2)+4*n1(2),10)];
n2Plot = [linspace(r2(1),r2(1)+4*n2(1),10);linspace(r2(2),r2(2)+4*n2(2),10)];
t1Plot = [linspace(r1(1)-0.5*t1(1),r1(1)+0.5*t1(1),10);linspace(r1(2)-0.5*t1(2),r1(2)+0.5*t1(2),10)];
t2Plot = [linspace(r2(1)-0.5*t2(1),r2(1)+0.5*t2(1),10);linspace(r2(2)-0.5*t2(2),r2(2)+0.5*t2(2),10)];
%cone plots
cone1 = [linspace(r2(1),r2(1)+2*(cos(alpha(2)-gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)-gamma)),10)];
cone1_ = [linspace(r2(1),r2(1)+2*(cos(alpha(2)+gamma)),10);linspace(r2(2),r2(2)+2*(sin(alpha(2)+gamma)),10)];
cone2 = [linspace(r1(1),r1(1)+2*cos(alpha(1)-gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)-gamma),10)];
cone2_ = [linspace(r1(1),r1(1)+2*cos(alpha(1)+gamma),10);linspace(r1(2),r1(2)+2*sin(alpha(1)+gamma),10)];

fig = figure;
set(fig, 'Position', [100, 100, 800, 600]);
hold on;
plot(rc(1),rc(2),'m.','markersize',20); 
plot(r1(1),r1(2),'k.','markersize',20); 
plot(n1Plot(1,:),n1Plot(2,:),'b:','linewidth',2);
plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
plot(total_force(1,:),total_force(2,:),'g','linewidth',2);
h = fill([-2,1,3,-1],[-1,5,-1,-1],'cyan');
plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
plot(r2(1),r2(2),'k.','markersize',20); 
plot(n2Plot(1,:),n2Plot(2,:),'b:','linewidth',2);
plot(t1Plot(1,:),t1Plot(2,:),'k','linewidth',2);
plot(t2Plot(1,:),t2Plot(2,:),'k','linewidth',2);
plot(cone1(1,:),cone1(2,:),'r','linewidth',2);
plot(cone1_(1,:),cone1_(2,:),'r','linewidth',2);
plot(cone2(1,:),cone2(2,:),'r','linewidth',2);
plot(cone2_(1,:),cone2_(2,:),'r','linewidth',2);
plot(total_force_(1,:),total_force_(2,:),'g','linewidth',2);
ylim([0,6])

set(h,'edgealpha',0,'facealpha',0.2) 
axis equal;
xlabel('X [m]', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Y [m]', 'Interpreter', 'latex', 'fontsize', 20);
set(gcf,'color','w');
title('Admissible f_x and External Force Visualization','fontsize',20);
grid on;
lgd = legend('COM','contact point','contact normal','friction cone','boundary line',' area of external force','Location','NorthWest');  
lgd.Interpreter = 'latex';  
lgd.FontSize = 20;  
saveas(gcf, 'Task IV.png');

%% Task 5

h = 5;
rc = [1;h];
mg = 1;
mu = 0.4;
gamma = atan(mu);

th_rez = 0.01*pi;
th = th_rez:th_rez:2*pi;
fd_max = range_fd_eq([r1 r2 rc],alpha,mu*ones(2,1),th);
fd_allowable = min(fd_max);

figure; 
plot(th,fd_max,'*-')
xticks(0:pi/2:3*pi); 
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
set(gcf,'color','w');
title('Maximum Disturbance Force Magnitude vs. Angle','fontsize',20);
ylabel('|F_d|','fontsize',16)
xlabel('[rad]','fontsize',16)
grid on;

h = 0:0.1:10;
fd_allowable = zeros(size(h));
for i = 1:length(h)
    rc = [1;h(i)];
    fd_max = range_fd_eq([r1 r2 rc],alpha,mu*ones(2,1),th);
    fd_allowable(i) = min(fd_max);
end

figure; 
plot(h,fd_allowable,'*-')
set(gcf,'color','w');
title('Maximum Disturbance Force Magnitude Across All Angles vs. COM height','fontsize',16);
ylabel('|F_d|','fontsize',16)
xlabel('h','fontsize',16)
grid on;

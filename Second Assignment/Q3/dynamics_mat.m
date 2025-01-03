function [M,B,W,W_d] = dynamics_mat(q,q_d)

    %Q3
    c = 0;

    m = 30;  
    l = 0.6;
    b = 0.2;
    d = 0.25;
    w = 0.25;
    Ic = 0.15;

    x = q(1);
    y = q(2);
    th = q(3);
    ph = q(4);
    x_d = q_d(1);
    y_d = q_d(2);
    th_d = q_d(3);
    ph_d = q_d(4);

    M = zeros(4,4);
    M(1,1) = m;
    M(2,2) = m;
    M(1,3) = -m*d*sin(th);
    M(3,1) = -m*d*sin(th);
    M(2,3) = m*d*cos(th);
    M(3,2) = m*d*cos(th);
    M(3,3) = Ic+m*d^2;
    
    B = zeros(4,1);
    B(1) = c*(-b*sin(ph+th)*(ph_d+th_d)+3*x_d-l*sin(th)*th_d)-m*d*th_d^2*cos(th);
    B(2) = c*(b*cos(ph+th)*(ph_d+th_d)+3*y_d+l*cos(th)*th_d)-m*d*th_d^2*sin(th);
    B(3) = c*(b^2*(ph_d+th_d)+l^2*th_d+0.5*w^2*th_d+b*cos(ph+th)*y_d-b*sin(ph+th)*x_d+l*cos(th)*y_d-l*sin(th)*x_d+b*l*cos(ph)*ph_d+2*b*l*cos(ph)*th_d);
    B(4) = b*c*(-x_d*sin(ph+th)+y_d*cos(ph+th)+b*(ph_d+th_d)+l*th_d*cos(ph));
    
    W = zeros(2,4);
    W(1,1) = -sin(th);
    W(1,2) = cos(th);
    W(2,1) = -sin(ph+th);
    W(2,2) = cos(ph+th);
    W(2,3) = b+l*cos(ph);
    W(2,4) = b;
    
    W_d = zeros(2,4);
    W_d(1,1) = -th_d*cos(th);
    W_d(1,2) = -th_d*sin(th);
    W_d(2,1) = -(ph_d+th_d)*cos(ph+th);
    W_d(2,2) = -(ph_d+th_d)*sin(ph+th);
    W_d(2,3) = -l*ph_d*sin(ph);
    
    
end
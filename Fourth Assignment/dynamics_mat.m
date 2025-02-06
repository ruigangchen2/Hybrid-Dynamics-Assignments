function [M,B,G,W,wn,wt] = dynamics_mat(q,q_d)

    R = 0.6;
    m1 = 5;
    m2 = 15;
    g = 10;
    h = 0.4;
    l = 0.2;
    J1 = 0.1*m1*R^2;
    J2 = 0.5*m2*l^2;

    % x = q(1);
    % y = q(2);
    th = q(3);
    ph = q(4);
    % x_d = q_d(1);
    % y_d = q_d(2);
    th_d = q_d(3);
    ph_d = q_d(4);

    M = zeros(4,4);
    M(1,1) = m1+m2;
    M(2,2) = m1+m2;
    M(3,1) = h*m1*cos(th)+l*m2*cos(th+ph);
    M(1,3) = M(3,1);
    M(1,4) = l*m2*cos(th+ph);
    M(4,1) = M(1,4);
    M(2,4) = l*m2*sin(th+ph);
    M(4,2) = M(2,4);
    M(4,4) = m2*l^2+J2;
    M(3,4) = M(4,4);
    M(4,3) = M(4,4);
    M(3,3) = M(4,4)+m1*h^2+J1;

    B = zeros(4,1);
    B(1) = -(h*m1*sin(th)*th_d+l*m2*sin(th+ph)*(th_d+ph_d))*th_d-l*m2*sin(th+ph)*(th_d+ph_d)*ph_d;
    B(2) = (h*m1*cos(th)*th_d+l*m2*cos(th+ph)*(th_d+ph_d))*th_d+l*m2*cos(th+ph)*(th_d+ph_d)*ph_d;
    
    G = zeros(4,1);
    G(2) = g*(m1+m2);
    G(3) = g*(m2*l*sin(th+ph)+m1*h*sin(th));
    G(4) = g*m2*l*sin(th+ph);

    wt = [1 0 R 0];
    wn = [0 1 0 0];
    W = [wt;wn];
    
end
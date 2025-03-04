function [M,B,G,W,Wtilde] = dynamics_mat(q,q_d)

    [m, mh, Ic, l, g, a] = model_params();

    th1 = q(3);
    th2 = q(4);
    th1_d = q_d(3);
    th2_d = q_d(4);

    M = zeros(4,4);
    M(1,1) = 2*m+mh;
    M(2,2) = M(1,1);
    M(1,3) = (3*m+2*mh)*l*cos(th1);
    M(3,1) = M(1,3);
    M(2,3) = -(3*m+2*mh)*l*sin(th1);
    M(3,2) = M(2,3);
    M(1,4) = l*m*cos(th2);
    M(4,1) = M(1,4);
    M(2,4) = l*m*sin(th2);
    M(4,2) = M(2,4);
    M(3,3) = Ic+5*l^2*m+4*l^2*mh;
    M(3,4) = 2*m*l^2*cos(th1+th2);
    M(4,3) = M(3,4);
    M(4,4) = m*l^2+Ic;

    B = zeros(4,1);

    B(1) = -3*l*m*sin(th1)*th1_d^2-l*m*sin(th2)*th2_d^2-2*l*mh*sin(th1)*th1_d^2;
    B(2) = -3*l*m*cos(th1)*th1_d^2+l*m*cos(th2)*th2_d^2-2*l*mh*cos(th1)*th1_d^2;
    B(3) = -2*m*l^2*(th2_d)^2*sin(th1+th2);
    B(4) = -2*m*l^2*(th1_d)^2*sin(th1+th2);
    
    G = zeros(4,1);
    G(1) = -g*sin(a)*(2*m+mh);
    G(2) = g*cos(a)*(2*m+mh);
    G(3) = -g*l*sin(a+th1)*(3*m+2*mh);
    G(4) = -m*g*l*sin(a-th2);

    wt = [1 0 0 0];
    wn = [0 1 0 0];
    W = [wt;wn];
    
    wtildet = [1 0 2*l*cos(th1) 2*l*cos(th2)];
    wtilden = [0 1 -2*l*sin(th1) 2*l*sin(th2)];
    Wtilde = [wtildet;wtilden];
    
end
function [M,B,W,W_d,Fd] = dynamics_mat(q,q_d)

    %Q3
    c = 0;

    m = 30;  
    l = 0.6;
    b = 0.2;
    d = 0.25;
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
    B(1) = -m*d*th_d^2*cos(th);
    B(2) = -m*d*th_d^2*sin(th);
    
    Fd = c*[x_d-b*(th_d+ph_d)*sin(th+ph)-l*th_d*sin(th)
            y_d+b*(th_d+ph_d)*cos(th+ph)+l*th_d*cos(th)
            -x_d*l*sin(th)+y_d*l*cos(th)+l^2*th_d+b^2*(th_d+ph_d)+b*l*th_d*cos(ph)
            -x_d*b*sin(th+ph)+y_d*b*cos(th+ph)+b^2*(th_d+ph_d)+b*l*th_d*cos(ph)];
    
    W = [-sin(th) cos(th) 0 0 
        -sin(th+ph) cos(th+ph) b+l*cos(ph) b];
    
    W_d = [-th_d*cos(th) -th_d*sin(th) 0 0 
        -(th_d+ph_d)*cos(th+ph) -(th_d+ph_d)*sin(th+ph) -l*ph_d*sin(ph) 0];
end
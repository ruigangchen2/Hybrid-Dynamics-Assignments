function [M,B,W,W_d] = dynamics_mat(q,q_d,c_term)

    c = c_term;

    m = 30;  
    l = 0.6;
    b = 0.2;
    d = 0.15;
    w = 0.25;
    Ic = 0.15;

    x = q(1);
    y = q(2);
    theta = q(3);
    phi = q(4);
  
    x_d = q_d(1);
    y_d = q_d(2);
    theta_d = q_d(3);
    phi_d = q_d(4);

    M=zeros(4);
    M(1,1)=m;
    M(1,3)=-d*m*sin(theta);
    M(2,2)=m;
    M(2,3)=d*m*cos(theta);
    M(3,1)=-d*m*sin(theta);
    M(3,2)=d*m*cos(theta);
    M(3,3)=Ic+d^2*m;

    B = zeros(4,1);
    B(1)=-c*(l*theta_d*sin(theta) - 3*x_d + b*phi_d*sin(phi + theta) + b*theta_d*sin(phi + theta)) - d*m*theta_d^2*cos(theta);
    B(2)=c*(3*y_d + l*theta_d*cos(theta) + b*phi_d*cos(phi + theta) + b*theta_d*cos(phi + theta)) - d*m*theta_d^2*sin(theta);
    B(3)=(c*(2*b^2*phi_d + 2*b^2*theta_d + 2*l^2*theta_d + theta_d*w^2 + 2*l*y_d*cos(theta) - 2*l*x_d*sin(theta) + 2*b*y_d*cos(phi + theta) - 2*b*x_d*sin(phi + theta) + 2*b*l*phi_d*cos(phi) + 4*b*l*theta_d*cos(phi)))/2;
    B(4)=b*c*(b*phi_d + b*theta_d + y_d*cos(phi + theta) - x_d*sin(phi + theta) + l*theta_d*cos(phi));

    W = zeros(2,4);
    W(1,1)=-sin(theta);
    W(1,2)=cos(theta);
    W(2,1)=-sin(phi + theta);
    W(2,2)=cos(phi + theta);
    W(2,3)=b + l*cos(phi);
    W(2,4)=b;

    W_d = zeros(2,4);
    W_d(1,1)=-theta_d*cos(theta);
    W_d(1,2)=-theta_d*sin(theta);
    W_d(2,1)=- phi_d*cos(phi + theta) - theta_d*cos(phi + theta);
    W_d(2,2)=- phi_d*sin(phi + theta) - theta_d*sin(phi + theta);
    W_d(2,3)=-l*phi_d*sin(phi);
end
function [M_matrix,B_matrix,G_matrix] = dynamics_mat(q,q_dot)

    m = 1; 
    g = 10; 
    l = 0.1; 

    theta = q(3,1);
    phi1 = q(4,1);
    phi2 = q(5,1);

    dtheta = q_dot(3,1);
    dphi1 = q_dot(4,1);
    dphi2 = q_dot(5,1);

    M_matrix = zeros(5,5);
    M_matrix(1,1) = 3*m;
    M_matrix(1,2) = 0;
    M_matrix(1,3) = l*m*(sin(phi1+theta)-sin(phi2+theta));
    M_matrix(1,4) = l*m*sin(phi1+theta);
    M_matrix(1,5) = -l*m*sin(phi2+theta);

    M_matrix(2,1) = 0;
    M_matrix(2,2) = 3*m;
    M_matrix(2,3) = -l*m*(cos(phi1+theta)-cos(phi2+theta));
    M_matrix(2,4) = -l*m*cos(phi1+theta);
    M_matrix(2,5) = l*m*cos(phi2+theta);

    M_matrix(3,1) = l*m*(sin(phi1+theta)-sin(phi2+theta));
    M_matrix(3,2) = -l*m*(cos(phi1+theta)-cos(phi2+theta));
    M_matrix(3,3) = 5*l^2*m + 2*l^2*m*cos(phi1) + 2*l^2*m*cos(phi2);
    M_matrix(3,4) = (4*l^2*m)/3+l^2*m*cos(phi1);
    M_matrix(3,5) = (4*l^2*m)/3+l^2*m*cos(phi2);

    M_matrix(4,1) = l*m*sin(phi1+theta);
    M_matrix(4,2) = -l*m*cos(phi1+theta);
    M_matrix(4,3) = (4*l^2*m)/3+l^2*m*cos(phi1);
    M_matrix(4,4) = 4*l^2*m/3;
    M_matrix(4,5) = 0;    

    M_matrix(5,1) = -l*m*sin(phi2+theta);
    M_matrix(5,2) = l*m*cos(phi2+theta);
    M_matrix(5,3) = (4*l^2*m)/3+l^2*m*cos(phi2);
    M_matrix(5,4) = 0;
    M_matrix(5,5) = 4*l^2*m/3;


    B_matrix = zeros(5,1);
    B_matrix(1,1) = l*m*cos(phi1 + theta)*dphi1^2 - l*m*cos(phi2 + theta)*dphi2^2 + ...
                    l*m*cos(phi1 + theta)*dtheta^2 - l*m*cos(phi2 + theta)*dtheta^2 + ...
                    2*l*m*cos(phi1 + theta)*dphi1*dtheta - 2*l*m*cos(phi2 + theta)*dphi2*dtheta;
    B_matrix(2,1) = l*m*sin(phi1 + theta)*dphi1^2 - l*m*sin(phi2 + theta)*dphi2^2 + ...
                    l*m*sin(phi1 + theta)*dtheta^2 - l*m*sin(phi2 + theta)*dtheta^2 + ...
                    2*l*m*sin(phi1 + theta)*dphi1*dtheta - 2*l*m*sin(phi2 + theta)*dphi2*dtheta;
    B_matrix(3,1) = -l^2*m*sin(phi1)*dphi1^2 - l^2*m*sin(phi2)*dphi2^2 - 2*l^2*m*sin(phi1)*dphi1*dtheta - 2*l^2*m*sin(phi2)*dphi2*dtheta;
    B_matrix(4,1) = l^2*m*sin(phi1)*dtheta^2;
    B_matrix(5,1) = l^2*m*sin(phi2)*dtheta^2;

    G_matrix = zeros(5,1);
    G_matrix(1,1) = 0;
    G_matrix(2,1) = 3*m*g; 
    G_matrix(3,1) = g*l*m*(cos(phi2+theta)-cos(phi1+theta));
    G_matrix(4,1) = -g*l*m*cos(phi1+theta);
    G_matrix(5,1) = g*l*m*cos(phi2+theta);

end
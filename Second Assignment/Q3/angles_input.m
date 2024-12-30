function [phi_matrix, phi_dot_matrix, phi_ddot_matrix] = angles_input(t)

    phi_matrix=zeros(1,1);
    phi_dot_matrix=zeros(1,1);
    phi_ddot_matrix=zeros(1,1);

    phi_matrix(1) = pi/4*cos(t);
    phi_dot_matrix(1) = -pi/4*sin(t);
    phi_ddot_matrix(1) = -pi/4*cos(t);
    
end
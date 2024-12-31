function [phi_matrix, phi_dot_matrix, phi_ddot_matrix] = angles_input(t)

    phi_matrix = (pi/4)*cos(t);
    phi_dot_matrix = -(pi/4)*sin(t);
    phi_ddot_matrix = -(pi/4)*cos(t);
    
end
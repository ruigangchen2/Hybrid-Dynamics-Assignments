function [phi_matrix, phi_dot_matrix, phi_ddot_matrix] = angles_input(t)

    phi_matrix = (pi/2)*cos(t);
    phi_dot_matrix = -(pi/2)*sin(t);
    phi_ddot_matrix = -(pi/2)*cos(t);
    
end
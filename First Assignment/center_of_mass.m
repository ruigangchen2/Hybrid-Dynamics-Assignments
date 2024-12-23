function [r_com, r_d_com]=center_of_mass(q, q_d)

    l = 0.1; 

    x = q(:,1); 
    y = q(:,2);
    theta = q(:,3); 
    phi1 = q(:,4); 
    phi2 = q(:,5);

    x_d = q_d(:,1); 
    y_d = q_d(:,2); 
    theta_d = q_d(:,3); 
    phi1_d = q_d(:,4); 
    phi2_d = q_d(:,5);

    r_com = zeros(2,1);
    r_com(1) = (x - 1/3*l*cos(theta+phi1)+1/3*l*cos(theta+phi2));  
    r_com(2) = (y - 1/3*l*sin(theta+phi1)+1/3*l*sin(theta+phi2));


    r_d_com = zeros(2,1);
    r_d_com(1) = (x_d + 1/3*l*(theta_d+phi1_d)*sin(theta+phi1) - 1/3*l*(theta_d+phi2_d)*sin(theta+phi2));
    r_d_com(2) = (y_d - 1/3*l*(theta_d+phi1_d)*cos(theta+phi1) + 1/3*l*(theta_d+phi2_d)*cos(theta+phi2));
end 

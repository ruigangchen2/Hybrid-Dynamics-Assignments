function [phi_matrix, phi_dot_matrix, phi_ddot_matrix] = angles_input(t)

    omega = 1;
    alpha = pi / 4; 
    beta = pi / 2; 
    psi = pi / 4; 
    t_f = 2 * pi / omega;

    phi_matrix=zeros(2,1);
    phi_dot_matrix=zeros(2,1);
    phi_ddot_matrix=zeros(2,1);


    if t <= t_f
        phi_matrix(1) = -(t^3*(alpha + beta*sin(psi - omega*t))*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5;
        phi_matrix(2) = (t^3*(alpha + beta*sin(psi + omega*t))*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5;

        phi_dot_matrix(1) = (beta*omega*t^3*cos(psi - omega*t)*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 - (3*t^2*(alpha + beta*sin(psi - omega*t))*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 - (t^3*(12*t - 15*t_f)*(alpha + beta*sin(psi - omega*t)))/t_f^5;
        phi_dot_matrix(2) = (t^3*(12*t - 15*t_f)*(alpha + beta*sin(psi + omega*t)))/t_f^5 + (3*t^2*(alpha + beta*sin(psi + omega*t))*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 + (beta*omega*t^3*cos(psi + omega*t)*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5;

        phi_ddot_matrix(1) = (2*beta*omega*t^3*cos(psi - omega*t)*(12*t - 15*t_f))/t_f^5 - (6*t^2*(12*t - 15*t_f)*(alpha + beta*sin(psi - omega*t)))/t_f^5 - (6*t*(alpha + beta*sin(psi - omega*t))*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 - (12*t^3*(alpha + beta*sin(psi - omega*t)))/t_f^5 + (6*beta*omega*t^2*cos(psi - omega*t)*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 + (beta*omega^2*t^3*sin(psi - omega*t)*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5;
        phi_ddot_matrix(2) = (12*t^3*(alpha + beta*sin(psi + omega*t)))/t_f^5 + (6*t^2*(12*t - 15*t_f)*(alpha + beta*sin(psi + omega*t)))/t_f^5 + (6*t*(alpha + beta*sin(psi + omega*t))*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 + (2*beta*omega*t^3*cos(psi + omega*t)*(12*t - 15*t_f))/t_f^5 + (6*beta*omega*t^2*cos(psi + omega*t)*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5 - (beta*omega^2*t^3*sin(psi + omega*t)*(6*t^2 - 15*t*t_f + 10*t_f^2))/t_f^5;

    else
        phi_matrix(1) = - alpha - beta*sin(psi - omega*t);
        phi_matrix(2) = alpha + beta*sin(psi + omega*t);

        phi_dot_matrix(1) = beta*omega*cos(psi - omega*t);
        phi_dot_matrix(2) = beta*omega*cos(psi + omega*t);

        phi_ddot_matrix(1) =  beta*omega^2*sin(psi - omega*t);
        phi_ddot_matrix(2) = -beta*omega^2*sin(psi + omega*t);

    end
end
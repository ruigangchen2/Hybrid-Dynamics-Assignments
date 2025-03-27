function [value, isterminal, direction] = events_stick(t,X)
    [~, ~, ~, l, ~, ~, mu] = model_params;
    [~, lambda]=dyn_sol_stick(t,X);
    
    y = X(2); 
    th1 = X(3); 
    th2 = X(4);
    
    lambda_t = lambda(1);
    lambda_n = lambda(2);
    
    %Slip onset of stance foot dir 1
    value(1) = -lambda(1)-mu*lambda(2);
    isterminal(1) = 1;
    direction(1) = 1;
    
    %Slip onset of stance foot dir 2
    value(2) = lambda(1)-mu*lambda(2);
    isterminal(2) = 1;
    direction(2) = 1;
    
    % Falling
    value(3) = y+2*l*cos(th1);
    isterminal(3) = 1;
    direction(3) = -1;

    %Swing foot impact 
    value(4) = th1-th2;
    isterminal(4) = 1;
    direction(4) = 1;
    

end
function dXdt = sys_stick(t, X)

    q = X(1:4);  
    dq = X(5:8);
    [ddq, ~] = dyn_sol_stick(t, X);
    dXdt = [dq;ddq];

end

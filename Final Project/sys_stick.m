function dXdt = sys_stick(t, X)

    dq = X(5:8);
    [ddq, ~] = dyn_sol_stick(t, X);
    dXdt = [dq;ddq];

end

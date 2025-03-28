function dXdt = sys_slip(t,X)
% equivalent to state_eq.m

q_d = X(5:8);

[q_dd,~] = dyn_sol_slip(t,X);

dXdt = [q_d;q_dd];
function dXdt = sys_slip(t,X)
% equivalent to state_eq.m

q = X(1:4);
q_d = X(5:8);

[q_dd,~] = dyn_sol_slip(q,q_d,t);

dXdt = [q_d;q_dd];
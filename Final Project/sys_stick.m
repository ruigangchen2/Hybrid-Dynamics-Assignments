function dXdt = sys_stick(t,X)
% equivalent to state_eq.m

q = X(1:4);
q_d = X(5:8);

[q_dd,~] = dyn_sol_stick(q,q_d,t);

dXdt = [q_d;q_dd];
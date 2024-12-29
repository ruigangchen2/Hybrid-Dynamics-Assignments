function X_d=state_eq(t,X)

q = X(1:4);
q_d = X(5:8);
[q_dd,~] = dyn_sol(q,q_d,t);
X_d = [q_d;q_dd];
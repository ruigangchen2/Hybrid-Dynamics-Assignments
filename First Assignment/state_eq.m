function X_d=state_eq(t,X)

qb = X(1:3);
qb_d = X(4:6);
[qs,qs_d,~]=angles_input(t); %here is where the redundancy comes in
q = [qb;qs];
q_d = [qb_d;qs_d];
[qb_dd,~]=dyn_sol(q,q_d,t);

X_d = [qb_d;qb_dd];
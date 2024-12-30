function X_d=state_eq(t,X,c_term)

qp = X(1:3);
qp_d = X(4:6);
[qa,qa_d,~]=angles_input(t);
q = [qp;qa];
q_d = [qp_d;qa_d];
[qp_dd,~,~]=dyn_sol(q,q_d,t,c_term);
X_d = [qp_d;qp_dd];
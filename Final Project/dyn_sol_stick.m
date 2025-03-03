function [q_dd, lambda]=dyn_sol_stick(q,q_d,t)

[M,B,G,W,~] = dynamics_mat(q,q_d);
lambda = (W*(M\(W')))\(W*(M\(B+G)));
q_dd = M\(W'*lambda-B-G);
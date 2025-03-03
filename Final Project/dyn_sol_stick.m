function [q_dd, Lambda]=dyn_sol_stick(q,q_d,t)

[M,B,G,W,~,~] = dynamics_mat(q,q_d);
Lambda = (W*(M\(W')))\(W*(M\(B+G)));
q_dd = M\(W'*Lambda-B-G);
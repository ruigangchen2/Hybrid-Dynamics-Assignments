function [q_dd, lambdan]=dyn_sol_slip(q,q_d,t)

global sgn_slip
global mu;

[M,B,G,~,wt,wn] = dynamics_mat(q,q_d);
sig = sgn_slip;

alpha = wn*(M\((wn-sig*mu*wt)'));
beta = wn*(M\(B+G));
lambdan = alpha/beta;
q_dd = M\(-B-G+(wn-sig*mu*wt)'*lambdan);

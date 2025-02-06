function [q_dd, lambdan]=dyn_sol_slip(q,q_d,t)

global sgn_slip
global mu;

[M,B,G,~,wn,wt] = dynamics_mat(q,q_d);
% alpha = wn*(M\((wn-sgn_slip*mu*wt)'));
% beta = wn*(M\(B+G));
alpha = (wn/M)*(wn-sgn_slip*mu*wt)';
beta = (wn/M)*(B+G);
% alpha = (wn)*inv(M)*(wn-sgn_slip*mu*wt)';
% beta = (wn)*inv(M)*(B+G);

lambdan = beta/alpha;
q_dd = M\(-B-G+(wn-sgn_slip*mu*wt)'*lambdan);
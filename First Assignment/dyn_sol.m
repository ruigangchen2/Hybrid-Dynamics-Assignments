function [qb_dd,tau]=dyn_sol(q,q_d,t)

b = 1:3;
s = 4:5;

% qs = q(s); redundancy in angles_input:
% qs_d = q_d(s); 1: we don't need qs and qs_d; 2: we can get them from this
[~,~,qs_dd]=angles_input(t);
[M,B,G]=dynamics_mat(q,q_d);
Mbb = M(b,b);
Mbs = M(b,s);
Mss = M(s,s);
Bb = B(b);
Bs = B(s);
Gb = G(b);
Gs = G(s);

qb_dd = Mbb\(-Mbs*qs_dd-Bb-Gb);
tau = Mbs'*qb_dd+Mss*qs_dd+Bs+Gs;
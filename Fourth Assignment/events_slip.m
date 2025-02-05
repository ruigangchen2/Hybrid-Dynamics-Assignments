function [value,isterminal,direction] = events_slip(t,X)

global sgn_slip
global mu

q = X(1:4);
q_d = X(5:8);
[M,B,G,~,wn,wt] = dynamics_mat(q,q_d);
alpha = wn*(M\((wn-sgn_slip*mu*wt)'));
beta = wn*(M\(B+G));
lambdan = alpha/beta;

vt = wt*q_d;

value(1) = lambdan; %consistent slip
isterminal(1) = 1;
direction(1) = -1;

% value(2) = -beta; %consistent separation
% isterminal(2) = 1;
% direction(2) = -1;

value(2) = vt; %change in slip dir
isterminal(2) = 1;
direction(2) = 0;


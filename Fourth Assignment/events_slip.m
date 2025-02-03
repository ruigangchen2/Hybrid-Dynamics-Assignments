function [value,isterminal,direction] = events_slip(t,X)
% event 1 is left slip (?) AKA when -lambda_t equals or exceeds mu*lambda_n
% event 2 is right slip AKA when lambda_t equals or exceeds mu*lambda_n

global sgn_slip
global mu

q = X(1:4);
q_d = X(5:8);
[M,B,G,~,wn,wt] = dynamics_mat(q,q_d);
alpha = wn*(M\((wn-sgn_slip*mu*wt)'));
beta = wn*(M\(B+G));
lambdan = alpha/beta;

value(1) = lambdan; %consistent slip
isterminal(1) = 1;
direction(1) = -1;

value(2) = -beta; %consistent separation
isterminal(2) = 1;
direction(2) = -1;


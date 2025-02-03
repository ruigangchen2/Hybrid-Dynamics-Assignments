function [value,isterminal,direction] = events_stick(t,X)
% event 1 is slipping right AKA when -lambda_t equals or exceeds mu*lambda_n
% event 2 is slipping left AKA when lambda_t equals or exceeds mu*lambda_n

q = X(1:4);
q_d = X(5:8);
[M,B,G,W,~,~] = dynamics_mat(q,q_d);
global mu
lambda = (W*(M\W'))\(W*(M\(B+G)));

value(1) = -lambda(1)-mu*lambda(2);
isterminal(1) = 1;
direction(1) = 1;

value(2) = lambda(1)-mu*lambda(2);
isterminal(2) = 1;
direction(2) = 1;


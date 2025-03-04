function [value,isterminal,direction] = events_stick(t,X)
% event 1 is slipping right AKA when -lambda_t equals or exceeds mu*lambda_n
% event 2 is slipping left AKA when lambda_t equals or exceeds mu*lambda_n

q = X(1:4);
q_d = X(5:8);
[M,B,G,W,~] = dynamics_mat(q,q_d);
[~, ~, ~, ~, ~, ~, mu] = model_params();
lambda = (W*(M\W'))\(W*(M\(B+G)));

y = q(2); th1 = q(3); th2 = q(4); l = 0.8;
yh = y+2*l*cos(th1);
ytilde = y+2*l*cos(th1)-2*l*cos(th2);

value(1) = -lambda(1)-mu*lambda(2);
isterminal(1) = 1;
direction(1) = 1;

value(2) = lambda(1)-mu*lambda(2);
isterminal(2) = 1;
direction(2) = 1;

%Falling
value(3) = yh;
isterminal(3) = 1;
direction(3) = -1;

%Impact
value(4) = ytilde;
isterminal(3) = 1;
direction(3) = -1;
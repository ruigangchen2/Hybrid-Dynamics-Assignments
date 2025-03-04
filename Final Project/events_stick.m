function [value,isterminal,direction] = events_stick(t,X)
% event 1 is slipping right AKA when -lambda_t equals or exceeds mu*lambda_n
% event 2 is slipping left AKA when lambda_t equals or exceeds mu*lambda_n

q = X(1:4);
q_d = X(5:8);
[M,B,G,W,~] = dynamics_mat(q,q_d);
% [m, mh, Ic, l, g, alpha, mu] = model_params
[~, ~, ~, l, ~, ~, mu] = model_params();
lambda = (W*(M\W'))\(W*(M\(B+G)));

y = q(2); th1 = q(3); th2 = q(4);
yh = y+2*l*cos(th1);
% ytilde = y+2*l*cos(th1)-2*l*cos(th2);

value(1) = -lambda(1)-mu*lambda(2);
isterminal(1) = 0; %Turned OFF for no-slip
direction(1) = 1;

value(2) = lambda(1)-mu*lambda(2);
isterminal(2) = 0; %Turned OFF for no-slip
direction(2) = 1;

%Falling
value(3) = yh;
isterminal(3) = 1;
direction(3) = -1;

%Impact of swing foot in front of stance foot - impact_law
value(4) = th2-th1;
isterminal(4) = 1;
direction(4) = -1;

%Impact of swing foot *behind* stance foot - scuffing onset
value(5) = th2-th1;
isterminal(5) = 0;
direction(5) = 1;

%Impact of swing foot "on top of" stance foot - end scuffing or about to
%fall
value(6) = th2+th1;
isterminal(6) = 0;
direction(6) = 1;
function [value,isterminal,direction] = events_slip(t,X)

global sgn_slip
global mu


[~, ~, ~, l, ~, ~] = model_params;    
y = X(2); 
th1 = X(3); 
th2 = X(4);
q_d = X(5:8);
[M, B, G, W, ~, ~] = dynamics_mat(X);
wt = W(1,:);
wn = W(2,:);
alpha = wn*(M\((wn-sgn_slip*mu*wt)'));
beta = wn*(M\(B+G));
lambdan = beta/alpha;

vt = wt*q_d;

%Change in slip dir
value(1) = vt*sgn_slip; 
isterminal(1) = 1;
direction(1) = -1;

%Separation
value(2) = lambdan; 
isterminal(2) = 1;
direction(2) = -1;

%Falling
value(3) = y+2*l*cos(th1);
isterminal(3) = 1;
direction(3) = -1;

% %Swing foot impact 
value(4) = th1-th2;
isterminal(4) = 0;
direction(4) = 1;

%Scuffing 
value(5) = th2-th1;
isterminal(5) = 0;
direction(5) = 1;
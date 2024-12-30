function [qp_dd, F_qa, lambda]=dyn_sol(q,q_d,t,c_term)

p = 1:3;
a = 4;

qp_d = q_d(p);
qa_d = q_d(a);

[~,~,qa_dd]=angles_input(t);
[M,B,W,W_d]=dynamics_mat(q,q_d,c_term);

Mpp = M(p,p);
Mpa = M(p,a);
Maa = M(a,a);

Bp = B(p);
Ba = B(a);

Wp = W(1:2,p);
Wa = W(1:2,a);

Wp_d = W_d(1:2,p);
Wa_d = W_d(1:2,a);

A = [Mpp, zeros(3,1), -Wp'; 
     Mpa', -eye(1,1), -Wa'; 
     Wp, zeros(2,1), zeros(2,2)];

B = -[Mpa*qa_dd + Bp; 
     Maa*qa_dd + Ba; 
     Wa*qa_dd + Wp_d*qp_d + Wa_d*qa_d];

dX = A\B;

qp_dd = dX(1:3);   
F_qa = dX(4);
lambda = dX(5:6);

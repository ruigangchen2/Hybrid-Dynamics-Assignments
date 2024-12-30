function [q_dd, Lambda]=dyn_sol(q,q_d,t)

[M,B,W,W_d]=dynamics_mat(q,q_d);
Fq = [0;0;0;0.5*cos(t)];

x = [M -W'; W zeros(2,2)]\[(Fq-B); -W_d*q_d];
q_dd = x(1:4);
Lambda = x(5:6);

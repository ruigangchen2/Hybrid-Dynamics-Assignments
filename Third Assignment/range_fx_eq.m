function [f_min,f_max]=range_fx_eq(r,alpha,mu)

mg = 1;

A = [-sin(alpha(1))-mu(1)*cos(alpha(1)) cos(alpha(1))-mu(1)*sin(alpha(1)) 0 0 0 
     sin(alpha(1))-mu(1)*cos(alpha(1)) -cos(alpha(1))-mu(1)*sin(alpha(1)) 0 0 0
     0 0 -sin(alpha(2))-mu(2)*cos(alpha(2)) cos(alpha(2))-mu(2)*sin(alpha(2)) 0
     0 0 sin(alpha(2))-mu(2)*cos(alpha(2)) -cos(alpha(2))-mu(2)*sin(alpha(2)) 0];

b = zeros(4,1);

Atilde = [1 0 1 0 1
          0 1 0 1 0 
          -r(2,1) r(1,1) -r(2,2) r(1,2) -r(2,3)];

btilde = [0; mg; mg*r(1,3)];
 
c = [0;0;0;0;1];

f_min = linprog(c,A,b,Atilde,btilde);
f_max = linprog(-c,A,b,Atilde,btilde);

end
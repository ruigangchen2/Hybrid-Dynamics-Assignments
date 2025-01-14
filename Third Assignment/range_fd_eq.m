function [fd_max]=range_fd_eq(r,alpha,mu,th)

mg = 1;

A = [-sin(alpha(1))-mu(1)*cos(alpha(1)) cos(alpha(1))-mu(1)*sin(alpha(1)) 0 0 0 
     sin(alpha(1))-mu(1)*cos(alpha(1)) -cos(alpha(1))-mu(1)*sin(alpha(1)) 0 0 0
     0 0 -sin(alpha(2))-mu(2)*cos(alpha(2)) cos(alpha(2))-mu(2)*sin(alpha(2)) 0
     0 0 sin(alpha(2))-mu(2)*cos(alpha(2)) -cos(alpha(2))-mu(2)*sin(alpha(2)) 0];

b = zeros(4,1);

fd_max = zeros(size(th));
for i = 1:length(th)
    Atilde = [1 0 1 0 cos(th(i))
              0 1 0 1 sin(th(i))
              -r(2,1) r(1,1) -r(2,2) r(1,2) -r(2,3)*cos(th(i))+r(1,3)*sin(th(i))];

    btilde = [0; mg; mg*r(1,3)];

    c = [0;0;0;0;1];

    f_max = linprog(-c,A,b,Atilde,btilde);
    
    if isempty(f_max)
        fd_max(i) = NaN;
    elseif f_max(end)>10
        fd_max(i) = 10;
    else
        fd_max(i) = f_max(end);
    end
end
end
function X_relabeled = relabel(X)

l = 0.8;
x = X(1);
y = X(2);
th1 = X(3);
th2 = X(4);
x_d = X(5);
y_d = X(6);
th1_d = X(7);
th2_d = X(8);

xtilde = x+2*l*sin(th1)+2*l*sin(th2);
ytilde = y+2*l*cos(th1)-2*l*cos(th2);
xtilde_d = x_d+2*l*th1_d*cos(th1)+2*l*th2_d*cos(th2);
ytilde_d = y_d-2*l*th1_d*sin(th1)+2*l*th2_d*sin(th2);

newx = xtilde;
newy = ytilde;
newth1 = th2;
newth2 = th1;
newx_d = xtilde_d;
newy_d = ytilde_d;
newth1_d = th2_d;
newth2_d = th1_d;

X_relabeled = [newx;newy;newth1;newth2;newx_d;newy_d;newth1_d;newth2_d];

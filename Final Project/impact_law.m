function Xnew= impact_law(Xold)
    [m, mh, Ic, l, ~, ~, ~] = model_params;
    [Mc, ~, ~, ~, ~, Wtilde] = dynamics_mat([Xold;NaN]);
    x = Xold(1);
    y = Xold(2);
    th1 = Xold(3);
    th2 = Xold(4);
    
    Xold_upd = (eye(4) - Mc\Wtilde'/(Wtilde/Mc*Wtilde')*Wtilde)*Xold(5:8);

    Xnew = [x + 2*l*sin(th1) + 2*l*sin(th2);
            y + 2*l*cos(th1) - 2*l*cos(th2);
            -th2;
            -th1;
            0;
            0;
            -Xold_upd(4);
            -Xold_upd(3);]; 
end



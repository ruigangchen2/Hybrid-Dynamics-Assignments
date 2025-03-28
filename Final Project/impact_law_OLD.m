function Xnew= impact_law_OLD(Xold)
    [m, mh, Ic, l, ~, ~, mu] = model_params;
    [Mc, ~, ~, ~, ~, Wtilde] = dynamics_mat([Xold;NaN]);
    wtilden = Wtilde(2,:);
    x = Xold(1);
    y = Xold(2);
    th1 = Xold(3);
    th2 = Xold(4);
    
    %Under plastic collision, LambdaHat simplifies to LambdaII (et=en=0)
    %plugging LambdaHat=LambdaII=Ac\(Wtilde*qdot) into impulse momentum
    %relation is equivalent to the below eqn:
    Ac = (Wtilde/Mc*Wtilde');
    betaStick = (eye(4) - Mc\Wtilde'/Ac*Wtilde);
    Xold_upd = betaStick*Xold(5:8);
     
    Xnew = [x + 2*l*sin(th1) + 2*l*sin(th2);
            y + 2*l*cos(th1) - 2*l*cos(th2);
            -th2;
            -th1;
            Wtilde*Xold_upd; %2 values; x_d and y_d of new stance foot, and 0 by design
            -Xold_upd(4);
            -Xold_upd(3);]; 
        
    [~, ~, ~, ~, ~, WtildeNew] = dynamics_mat([Xnew;NaN]);
    pdot = WtildeNew*Xnew(5:8);

    if pdot(2) <= 0
        error('failure; no rear foot separation');
    end
end



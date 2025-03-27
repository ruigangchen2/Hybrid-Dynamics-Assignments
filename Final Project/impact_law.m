function Xnew= impact_law(Xold)
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

%         LambdaII = -Ac\(Wtilde*Xold(5:8));
%     Xold_upd1 = Xold(5:8)+Mc\Wtilde'*LambdaII;
%     slipping = abs(lambda(1))-mu*lambda(2);
%     if slipping bla bla 
%     Gamma = [1;(sgn_slip? sgn_slip_tilde?)*mu]
%     betaSlip = (eye(4) - (Mc\Wtilde'*Gamma)/(wtilden*Xold(1:4)'/Mc*Wtilde'*Gamma)*wtilden*Xold(1:4)');
%     Xold_upd = betaSlip*Xold(5:8);
     
    Xnew = [x + 2*l*sin(th1) + 2*l*sin(th2);
            y + 2*l*cos(th1) - 2*l*cos(th2);
            -th2;
            -th1;
            Wtilde*Xold_upd; %2 values; x_d and y_d of new stance foot, and 0 by design
            -Xold_upd(4);
            -Xold_upd(3);]; 
        
%     if Xnew(6) <= 0
%         error('failure; no rear foot separation');
%     end
end



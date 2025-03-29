function Xnew= impact_law(Xold)
    global mu
    [m, mh, Ic, l, ~, ~] = model_params;
    [Mc, ~, ~, ~, ~, Wtilde] = dynamics_mat([Xold;NaN]);
    wtilden = Wtilde(2,:);
    x = Xold(1);
    y = Xold(2);
    th1 = Xold(3);
    th2 = Xold(4);
    pdot = Wtilde*Xold(5:8)

    Ac = (Wtilde/Mc*Wtilde');
    LambdaI = [0; -pdot(2)/Ac(2,2)];
    LambdaII = -Ac\pdot;
    LambdaHat = LambdaII; %Since CoR's are 0
    
    if abs(LambdaHat(1)) <= mu*LambdaHat(2)
        kappa = 1;
    else
        kappa = mu*LambdaI(2)/(abs(LambdaII(1))-mu*(LambdaII(2)-LambdaI(2)));
    end
    
    Lambda = LambdaI+kappa*(LambdaII-LambdaI);
%     %Check when kappa not = 1 that these are equal 
%     Lambda(1)
%     Lambda(2)*mu*sign(LambdaHat(1))
    Deltaq_d = (Mc\(Wtilde'))*Lambda;
    qplus_d = Deltaq_d+Xold(5:8);

    Xnew = [x + 2*l*sin(th1) + 2*l*sin(th2);
            y + 2*l*cos(th1) - 2*l*cos(th2);
            -th2;
            -th1;
            Wtilde*qplus_d; %2 values; x_d and y_d of new stance foot, and 0 by design (in stick)
            -qplus_d(4);
            -qplus_d(3);]; 

    [~, ~, ~, ~, ~, WtildeNew] = dynamics_mat(Xnew);
    pdotNew = WtildeNew*Xnew(5:8);
%     pdotNew = Wtilde*Xnew(5:8)?

    if pdotNew(2) <= 0
        warning('failure; double-foot impact');
        Xnew = NaN(size(Xnew));
    end
end



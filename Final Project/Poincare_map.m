function Znew = Poincare_map(Zold)

    status = 0; 
    t_start = 0;
    t_stop = 10;
    X0 = [0, 0, Zold(1), Zold(1), 0, 0, Zold(2), Zold(3)]; 
    options = odeset('reltol', 1e-8, 'abstol', 1e-8, 'Events', @(t, X) events_stick(t, X));
    while status < 1  
        [t, X, ~, ~, ie] = ode45(@(t, X) sys_stick(t, X), [t_start, t_stop], X0, options);
        X0 = X(end, :).';
        t_start = t(end);
        if ie(end) == 3   % hip collision
            Znew = ones(3, 1) * 100;
            return;
        elseif ie(end) == 4
            status = 1; 
            Xold = X0;
        end
    end
    if status == 1
        temp = impact_law(Xold);
        Znew = temp([3, 7, 8]);
    end

end
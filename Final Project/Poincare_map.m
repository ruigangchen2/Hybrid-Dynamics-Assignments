function Znew = Poincare_map(Zold)

    % (By design, post-impact new stance foot veloc. = 0)
    %                             |  |
    % (arbitrary; periodic solns  |  | don't require periodic x,y)
    %     |  |                    |  |
    %     V  V                    V  V
    X0 = [0, 0, Zold(1), Zold(1), 0, 0, Zold(2), Zold(3)]; 
    t_start = 0;
    t_stop = 10;
    options = odeset('reltol', 1e-8, 'abstol', 1e-8, 'Events', @(t, X) events_stick(t, X));
    [~, X, ~, ~, ie] = ode45(@(t, X) sys_stick(t, X), [t_start, t_stop], X0, options);

    if ie(end) == 3   % hip collision
        error('failure; falling')
    end
   
    temp = impact_law(X(end, :).');
    Znew = temp([3, 7, 8]);


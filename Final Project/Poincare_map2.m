function Znew = Poincare_map2(Zold)

global sgn_slip
sgn_slip = sign(Zold(4));

    % (By design, post-impact new stance foot normal veloc. = 0)
    %                                      |
    % (arbitrary; periodic solns don't     | require periodic x,y)
    %     |  |                             |
    %     V  V                             V
    Xe = [0, 0, Zold(1), Zold(1), Zold(4), 0, Zold(2), Zold(3)]; 
    te = 0;
    finalTime = 1000;
    op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         
    op_slip = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_slip);         
    firstIter = true;

    while te<finalTime
        %Check slipping or sticking
        [value, ~, direction] = events_stick(te,Xe);
        if value(1:3).*direction(1:3) < 0
            [t,X,te,Xe,ie] = ode45(@sys_stick,[te, finalTime], Xe, op_stick);
            te = te(end);
            Xe = Xe(end,:);
            ie = ie(end);
            
            %Update state
            switch ie
                case 1
                    sgn_slip = 1;
                case 2
                    sgn_slip = -1;
                case 3
                    warning('failure; falling')
                    Znew = NaN(1,4);
                    break;
                case 4
                    temp = impact_law(Xe');
%                     temp = impact_law(X(end, :).');
                    Znew = temp([3, 7, 8, 5]);
                    break;
            end
        else
            [t,X,te,Xe,ie] = ode45(@sys_slip, [te, finalTime], Xe, op_slip);
            te = te(end);
            Xe = Xe(end,:);
            ie = ie(end);
            
            %Update state
            switch ie(end)
                case 1
                    sgn_slip = -sgn_slip;
                case 2
                    warning('failure; stance foot separation')
                    Znew = NaN(1,4);
                case 3
                    warning('failure; falling')
                    Znew = NaN(1,4);
                case 4
                    temp = impact_law(Xe');
%                     temp = impact_law(X(end, :).');
                    Znew = temp([3, 7, 8, 5]);
                    break;
            end
        end
    end
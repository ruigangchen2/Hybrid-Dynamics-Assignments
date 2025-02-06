function [finalTimes,finalX,finalLambda,finalq_dd,stickInds,reachedEndTime] = solveHybridDynamics(w0)

global R
global mu

op_stick = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_stick);         
op_slip = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'Events',@events_slip); 
endTime = 10;
reachedEndTime = 1;
te = 0;
ye = [0 R 0 0 0 0 0 w0];
finalX = [];
finalLambda = [];
finalTimes = [];
finalq_dd = []; %for cross-checking
stickInds = logical([]);

while ~isempty(te)
    if events_stick(te,ye)<0
        %Solve corresponding ODE
        [t,X,te,ye,ie] = ode45(@sys_stick, [te endTime], ye, op_stick);
        Lambda = zeros(length(t),2);
        q_dd = zeros(length(t),4);
        for i = 1:length(t)
            [q_dd(i,:),Lambda(i,1:2)] = dyn_sol_stick(X(i,1:4)',X(i,5:8)',t(i));
        end
        
        %Record results
        finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))]; %transition state will be considered to belong to the next state
        finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))]; %NaNs added for plotting disconinuities
        finalTimes = [finalTimes;t(1:end-1);NaN];
        finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];

        %Record and update state
        stickInds = [stickInds; true(length(t),1)];
        sgn_slip = sign(1.5-ie);
    else
        %Solve corresponding ODE
        [t,X,te,ye,ie] = ode45(@sys_slip, [te endTime], ye, op_slip);
        Lambda = zeros(length(t),2);
        q_dd = zeros(length(t),4);
        for i = 1:length(t)
            [q_dd(i,:),Lambda(i,2)] = dyn_sol_slip(X(i,1:4)',X(i,5:8)',t(i));
            Lambda(i,1) = -sgn_slip*mu*Lambda(i,2);
        end

        %Record results
        finalX = [finalX;X(1:end-1,:);NaN(1,size(X,2))];
        finalLambda = [finalLambda;Lambda(1:end-1,:);NaN(1,size(Lambda,2))];
        finalTimes = [finalTimes;t(1:end-1);NaN];
        finalq_dd = [finalq_dd;q_dd(1:end-1,:);NaN(1,size(q_dd,2))];
        
        %Record and update state
        stickInds = [stickInds; false(length(t),1)];
        if length(ie)>1
            [~,firstEvent] = min(te);
            ie = ie(firstEvent);
            te = te(firstEvent);
            ye = ye(firstEvent,:);
        end
        
        if ie-1
            reachedEndTime = 0;
            break;
        end
        sgn_slip = -sgn_slip;
    end
end

finalX(end,:) = X(end,:);
finalLambda(end,:) = Lambda(end,:);
finalTimes(end,:) = t(end);
finalq_dd(end,:) = q_dd(end,:);
function [q_dd, lambda]=dyn_sol_stick(t,X)

    [M, B, G, W, ~, ~] = dynamics_mat(X);    
    lambda = (W*(M\(W')))\(W*(M\(B+G)));
    q_dd = M\(W'*lambda-B-G);
    
end


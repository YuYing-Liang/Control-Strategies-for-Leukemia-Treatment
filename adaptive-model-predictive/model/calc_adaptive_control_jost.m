function [f] = calc_adaptive_control_jost(p_star, rho, x_k, u, q, r, Hp, Hu, k, n_cycle, n_drug)
    % In this case, Hp and Hu are the number of cycles NOT number of days
    % k is the current cycle (we are predicting u for cycle k+1)
    f=0;
    for i=1:Hp
        % this is the estimated model (from calc_p_star)
        % we pass in u to f_star which is the vector we want to minimize
        for j=1:n_cycle
            t = n_cycle*(k-1) + n_cycle*(i-1) + j;
            x_k = jost(t, x_k, u(i), p_star, n_cycle, n_drug);
            z = g_nonlin(x_k);
            if k+t <= length(rho(1,1,:))
                y_err = z - rho(:,:,k+t);
                f = f + q*norm(y_err);
            end
        end
    end
    
    for i=1:Hu
        f = f + r*u(i)^2;
    end
end
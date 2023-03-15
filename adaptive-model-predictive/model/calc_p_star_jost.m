function [f] = calc_p_star_jost(p_guess, p_actual, x_init, u, N, n_cycle, n_drug)
    % N is the number of measurements up to this point
    f=0;
    x_k = x_init;
    x_k_obs = x_init;
    k = 1;
    
    for i=1:N 
        % prediction
        % estimated parameters go in here
        x_k = jost(i, x_k, u(k), p_guess, n_cycle, n_drug); % f_star
        y_k = g_nonlin(x_k);

        % true system
        % true parameters that the controller doesn't know
        x_k_obs = jost(i, x_k_obs, u(k), p_actual, n_cycle, n_drug); % f_actual
        y_k_obs = g_nonlin(x_k_obs);
        
        f = f + norm(y_k - y_k_obs);
        
        if mod(i, n_cycle) == 0
            k = k + 1;
        end
    end
    
    f = (1/N)*f;
end
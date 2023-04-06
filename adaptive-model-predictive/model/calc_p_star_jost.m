function [f] = calc_p_star_jost(p_guess, p_actual, x_init, u, N_days, n_cycle, n_drug, t_s)
    % N is the number of measurements up to this point
    f=0;
    x_k = x_init;
    x_k_obs = x_init;
    k = 1;

    if N_days == 1
        N_steps = 1;
    else
        N_steps = N_days/t_s;
    end

    for i=1:N_steps
        % prediction
        % estimated parameters go in here
        x_k = jost_discrete(i*t_s, x_k, u(k), p_guess, n_cycle, n_drug, t_s); % f_star
        y_k = g_nonlin(x_k);

        % true system
        % true parameters that the controller doesn't know
        x_k_obs = jost_discrete(i*t_s, x_k_obs, u(k), p_actual, n_cycle, n_drug, t_s); % f_actual
        y_k_obs = g_nonlin(x_k_obs);
        
        f = f + norm(y_k - y_k_obs);
        
        if get_cycle_day(i*t_s, n_cycle) == n_cycle
            k = k + 1;
        end
    end
    
    f = (1/N_steps)*f;
end
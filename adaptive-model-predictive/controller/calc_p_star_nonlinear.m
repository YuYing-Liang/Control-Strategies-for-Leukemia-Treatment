function [f] = calc_p_star_nonlinear(p_guess, p_actual, x_init, mn, u, k, t_s)
    f=0;
    w=[1,1,1,1];
    x_k = x_init;
    x_k_obs = x_init;
    
    for i=1:k
        % prediction
        % estimated parameters go in here
        x_k = f_nonlin(x_k, u(i), p_guess, t_s); % f_star
        y_k = g_nonlin(x_k);
        
        % true system
        % true parameters that the controller doesn't know
        x_k_obs = f_nonlin(x_k_obs, u(i), p_actual, t_s); % f_actual
        y_k_obs = g_nonlin(x_k_obs) + mn(:,:,i);
        
        for j=1:length(y_k)
            f = f + ((y_k(j,1) - y_k_obs(j,1))/w(j))^2;
        end
    end
    
    f = (1/k)*f;
end
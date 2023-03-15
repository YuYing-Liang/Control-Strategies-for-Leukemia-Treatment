function [f] = calc_p_star(A_star, A_actual, b, C, x_init, mn, u, k)
    f=0;
    w=[1,1,1];
    x_k = x_init;
    x_k_obs = x_init;
    
    for i=1:k
        % prediction
        % estimated parameters go in here
        x_k = A_star*x_k + b*u(i);
        y_k = C*x_k;
        
        % true system
        % true parameters that the controller doesn't know
        x_k_obs = A_actual*x_k_obs + b*u(i);
        y_k_obs = C*x_k_obs + mn(:,:,i);
        
        for j=1:length(y_k)
            f = f + ((y_k(j,1) - y_k_obs(j,1))/w(j))^2;
        end
    end
    
    f = (1/k)*f;
end
function [f] = calc_adaptive_control_nonlinear(p_star, rho, x_k, u, q, r, Hp, Hu, k, t_s)
    f=0;
    for i=1:Hp
        % this is the estimated model (from calc_p_star)
        % we pass in u to f_star which is the vector we want to minimize
        x_k = f_nonlin(x_k, u(i), p_star, t_s);
        z = g_nonlin(x_k);
        if k+i <= length(rho(1,1,:))
            y_err = z - rho(:,:,k+i);
            f = f + q*norm(y_err);
        end
    end
    
    for i=1:Hu
        f = f + r*u(i)^2;
    end
end
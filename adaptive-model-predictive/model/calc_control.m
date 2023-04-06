function [f] = calc_control(p_star, rho, x_k, curr_k, u, q, r, Hp, Hu, n_cycle, n_drug, t_s)
    % In this case, Hp and Hu are the number of days
    f=0;
    n_steps_per_day = 1/t_s;
    cycle = 1;
    for k=curr_k+1:curr_k+Hp
        day_in_cycle = get_cycle_day(k, n_cycle);
        u_applied = 0;
        if day_in_cycle < n_drug
            if cycle > length(u)
                break
            end
            u_applied = u(cycle);
        end
        for i=1:n_steps_per_day
            j = n_steps_per_day*(k-1) + i;
            x_k = jost_discrete(j*t_s, x_k, u_applied, p_star, n_cycle, n_drug, t_s);
            z = g_nonlin(x_k);
            if j <= length(rho)
                y_err = z - rho(j, :);
                f = f + q*norm(y_err);
            end
        end
        if get_cycle_day(k, n_cycle) == n_cycle
            cycle = cycle + 1;
        end
    end
    
%     for i=2:Hu
%         f = f + r*norm(u(i)-u(i-1));
%     end
end
% function [f] = calc_adaptive_control_jost(p_star, rho, x_k, u, u_calc, q, r, Hp, Hu, k, n_cycle, n_drug, t_s, nominal)
%     % In this case, Hp and Hu are the number of cycles NOT number of days
%     % k is the current cycle (we are predicting u for cycle k+1)
%     f=0;
%     for i=1:Hp % i is the cycle into the future
%         if k == 1 && i == 1
%             start_step = 1;
%             end_step = start_step + n_cycle/t_s - 2;
%         else
%             % n_cycle is # days/cycle, (k+i-2) is # of cycles
%             total_days = n_cycle*(k+i-2);
%             start_step = total_days/t_s;
%             end_step = start_step + n_cycle/t_s - 1;
%         end
%         for j=start_step:end_step
%             % this is the estimated model (from calc_p_star)
%             % we pass in u to f_star which is the vector we want to minimize
%             % x_k at step 100, is rho at step 102 --> look into this more
%             % decrease the number steps (1000 -> 10)
%             % disp(j);
%             x_k = jost_discrete((j+1)*t_s, x_k, u(i), p_star, n_cycle, n_drug, t_s);
%             z = g_nonlin(x_k);
%             if j+1 <= length(rho)
%                 y_err = z - rho(j+1, :);
% %                 if norm(y_err) < 1e-4
% %                     y_err = 0;
% %                     x_k(8) = rho(j+1);
% %                 end
% %                 if y_err ~= 0
% %                     disp(y_err)
% %                 end
%                 f = f + q*norm(y_err);
%             end
%         end
%     end
%     
%     for i=2:Hu
%         f = f + r*norm(u(i)-u(i-1));
%     end
%     
%     for i=k-1:-1:1
%         f = f + r*norm(u(1) - u_calc(i));
%     end
% 
%     if u(1) == nominal
%         fprintf("reached nominal! %d", nominal)
%         disp(f)
%         disp(u)
%     end
% end
function [f] = calc_adaptive_control_jost(p_star, rho, x_k, curr_k, u, q, r, Hp, Hu, n_cycle, n_drug, t_s)
    % In this case, Hp is the number of days while Hu is the number of cycles
    f=0;
    n_steps_per_day = 1/t_s;
    cycle = 1;
    for k=curr_k+1:curr_k+Hp - 1
        for i=1:n_steps_per_day
            j = n_steps_per_day*(k-1) + i;
           
            x_k = jost_discrete(j*t_s, x_k, u(cycle), p_star, n_cycle, n_drug, t_s);
            z = g_nonlin(x_k);
            
            if j <= length(rho) 
                y_err = z - rho(j,:);
                f = f + q*norm(y_err);
            end
        end
        if get_cycle_day(k, n_cycle) == n_cycle
            cycle = cycle + 1;
        end
    end
    
    for i=1:Hu
        f = f + r*norm(u(i));
    end
end
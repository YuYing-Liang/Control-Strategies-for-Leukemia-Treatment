function [f] = calc_control_comparison(p_star, rho, x_k, curr_k, curr_cycle, u, q, r, Hp, Hu, patient_cycles, t_s)
    % Here, Hp and Hu are in units of number of cycles
    % patient_cycles is a lookup table of a list of cycles in order
    f=0;
    n_steps_per_day = 1/t_s;
    cycle = 1;
    days = 1;
    days_to_predict = 0;
    n_cycle = patient_cycles(curr_cycle,3);
    n_drug = n_cycle;
    
    % calculating the total number of days in the predicted number of
    % cycles
    for k=curr_cycle:curr_cycle+Hp - 1
        if k <= length(patient_cycles)
            days_to_predict = days_to_predict + patient_cycles(k,3);
        end
    end
    
    for k=curr_k+1:curr_k+days_to_predict   % units of days
        for i=1:n_steps_per_day             % units of steps
            j = n_steps_per_day*(k-1) + i;
           
            x_k = jost_discrete(j*t_s, x_k, u(cycle), p_star, n_cycle, n_drug, t_s);
            z = g_nonlin(x_k);
            
            if j <= length(rho) 
                y_err = z - rho(j,:);
                f = f + q*norm(y_err);
            end
        end
        if days > n_cycle
            cycle = cycle + 1;
            days = 1;
            if curr_cycle+cycle-1 <= length(patient_cycles)
                n_cycle = patient_cycles(curr_cycle+cycle-1,3);
                n_drug = n_cycle;
            else
                break
            end
        else
            days = days + 1; 
        end
    end
    
    for i=1:Hu
        f = f + r*norm(u(i));
    end
    

    % In this case, Hp is the number of days while Hu is the number of cycles
%     f=0;
%     n_steps_per_day = 1/t_s;
%     cycle = 1;
%     day = 1;
%     for k=curr_k+1:curr_k+Hp - 1
%         for i=1:n_steps_per_day
%             j = n_steps_per_day*(k-1) + i;
%            
%             x_k = jost_discrete(j*t_s, x_k, u(cycle), p_star, Hp, Hp, t_s);
%             z = g_nonlin(x_k);
%             
%             if j <= length(rho) 
%                 y_err = z - rho(j,:);
%                 f = f + q*norm(y_err);
%             end
%         end
%         if day == Hp
%             cycle = cycle + 1;
%         end
%         day = day + 1;
%     end
%     
%     for i=1:Hu
%         f = f + r*norm(u(i));
%     end
end
function [f] = calc_control_comparison(p_star, rho, x_k, curr_k, curr_cycle, u, q, r, Hp, Hu, patient_cycles, t_s)
    % Argument Description
    % p_star    is the estimated parameter vector from calc_params_comparison
    % rho       is the set point (desired number of neutrophils) or nominal trajectory
    % x_k       is the last calculated state (using the controller and true parameter vector)
    % curr_k    is the current day in the simulation
    % curr_cycle is the current cycle in the simulation
    % u         is the variable to optimize
    % q, r      penalization factors
    % Hp        is the prediction horizon (in units of cycles)
    % Hu        <check noble for term> (in units of cycles) - not as important
    % patient_cycles    is a vector of cycle lengths in the entire treatment
    %       --> patient_cycles(:,1) is the start time of cycle
    %       --> patient_cycles(:,2) is the end time of cycle
    %       --> patient_cycles(:,3) is the cycle length in days
    %    note:  patient_cycles(1,:) is the first cycle, 
    %           patient_cycles(2,:) is the second cycle, etc.
    % ts        is the time step in days

    f=0;
    n_steps_per_day = 1/t_s;
    cycle = 1;
    days = 1;
    days_to_predict = 0;
    n_cycle = patient_cycles(curr_cycle,3);
    n_drug = n_cycle;
    
    % calculating the total number of days in the predicted number of cycles
    % see example in documentation for more details
    for k=curr_cycle:curr_cycle+Hp - 1
        if k <= length(patient_cycles)
            days_to_predict = days_to_predict + patient_cycles(k,3);
        end
    end
    
    for k=curr_k+1:curr_k+days_to_predict   % units of days
        for i=1:n_steps_per_day             % units of steps
            j = n_steps_per_day*(k-1) + i;
           
            % Note: we are using the best possible guess for the parameter
            % vector of the virtual patient (p_star), we are NOT using p_actual
            
            % This line is take the initial state and predicting into the
            % future (repetitively) to figure out the best control input
            % for the next treatment cycle.
            
            % IMPORTANT NOTE: In practice x_k should only rely on past 
            % measurements from the patient but we do not have access 
            % to the full state in these measurements. So an important step
            % in the future is to build an observer to do so.
            x_k = jost_discrete(j*t_s, x_k, u(cycle), p_star, n_cycle, n_drug, t_s);
            y = g_nonlin(x_k);
            
            if j <= length(rho) 
                y_err = y - rho(j,:);
                f = f + q*norm(y_err);
            end
        end
        
        % This is statement is used to keep track of how many cycles we've 
        % predicted ahead so far.
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
end
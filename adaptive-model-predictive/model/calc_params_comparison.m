function [f] = calc_params_comparison(past_mmts, p_guess, x_init, u, prev_step, n_cycle, n_drug, t_s)
    % past_mmts is a vector of previous "measurements" taken from the virtual patient
    %       --> size: 3 x 1 x number of "measurements"
    %       --> (1,:,:) is the "measurement"
    %       --> (2,:,:) is the index (aka the time that it was taken)
    %       --> (3,:,:) is the time in days (since the treatment started) that the measurement was taken
    % p_guess   is the variable that we're trying to optimize
    % x_init    is the estimated initial state for the virtual patient
    % u         is a vector of control inputs up to this point in the simulation
    % prev_step is the number of steps up to this point
    % n_cycle   is the number of days in this cycle
    % n_drug    is the number of "on drug" days in this cycle
    % t_s       is the timestep in units of days
   
    % x_k       is the prediction: estimated parameters go in here
    % y_k_obs   is the true system: has the true parameters that the controller doesn't know

    f=0;
    x_k = x_init;
    k = 1;
    mmt_i = 1;
    
    for i=1:prev_step
        % Note: Since i is the number of time steps and t_s is in units of days we need to
        % multiply i by t_s to get a unit in days
        % This is updating the state using the parameter guess
        x_k = jost_discrete((i+1)*t_s, x_k, u(k), p_guess, n_cycle, n_drug, t_s);
        y_k = g_nonlin(x_k);
        
        if i == past_mmts(2,:,mmt_i) % when the index matches measurement time
            % use measurement from past_mmts, this is the measurement from the true virtual patient
            y_k_obs = past_mmts(1,:,mmt_i);

            % update cost variable
            f = f + norm(y_k - y_k_obs);
            mmt_i = mmt_i + 1;
            
            % break out of loop after last measurement is evaluated
            if mmt_i > size(past_mmts,3)
                break
            end
        end
        
        if get_cycle_day((i+1)*t_s, n_cycle) == n_cycle
            k = k + 1;
        end
    end
    
    % normalize the cost variable by the number of measurements taken
    % so far from virtual patient
    f = (1/size(past_mmts,3))*f;
end
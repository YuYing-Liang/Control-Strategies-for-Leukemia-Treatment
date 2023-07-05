function [d] = get_cycle_day(day, N_days_in_cycle)
    d = mod(day, N_days_in_cycle);
    if d == 0
        d = N_days_in_cycle;
    end
end
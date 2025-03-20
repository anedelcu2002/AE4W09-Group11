function [r_total] = gbx_ratio_total(Omega_LSS_max, Omega_LSS_min, Omega_HSS_max, Omega_HSS_min, r_limit)
    % total Gearbox ratio
    r_max = Omega_HSS_max/Omega_LSS_max;
    r_min = Omega_HSS_min/Omega_LSS_min;

    if r_min > r_max
        disp('Cant find propper scaling ratio')
    end

    r_total = r_min + 0.5*(r_max - r_min);

    if r_total > r_limit
        disp('Total Gearbox ratio exceeds set limit')
    end
end
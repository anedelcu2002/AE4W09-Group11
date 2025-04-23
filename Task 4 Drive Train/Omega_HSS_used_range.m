function [Omega_min, Omega_max] = Omega_HSS_used_range(r_total, Omega_LSS_min, Omega_LSS_max)
    % function that calculates the Omega HSS range
    Omega_min = r_total*Omega_LSS_min;
    Omega_max = r_total*Omega_LSS_max; 
end
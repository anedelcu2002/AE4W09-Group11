function [Omega_min, Omega_max] = Omega_HSS_allowed_range(n, f_min, f_max)
    % function that calculates the Omega HSS range
    Omega_min = (4*pi*f_min)/(n);
    Omega_max = (4*pi*f_max)/(n); 
end
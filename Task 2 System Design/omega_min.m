function [omeg_min] = omega_min(U_ci, lambda, D)
% function calcultes the minimal rotational speed of the rotor in rad/s based on the tip speed ratio,
% the cut-in wind speed in m/s and the rotor diameter in m

    omeg_min = (lambda * U_ci) / (D/2);

end
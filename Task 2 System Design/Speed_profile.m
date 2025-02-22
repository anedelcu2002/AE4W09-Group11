function [U_array]=Speed_profile(U_0, z_0, alpha, h_0, h_hub)
    % calculates the speed at the new hub height h_hub in m, based on the original speed U_0 in m/s at the original height h_0 in m and z_0 in m

    if h_hub<60 %assuming that the measured speed is real and not potential!
        U_array=U_0*(log(h/z_0)/log(h_0/z_0));
    else
        U_meso=U_0*(log(60/z_0)/log(h_0/z_0));
        U_array=U_meso*(h_hub/60)^alpha
    end
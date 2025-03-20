function [V_gen] = gen_sizing(F_d, P_r, Omega_max, eff_gbx)
    % function that calculates the generator volume
        V_gen= (P_r*eff_gbx)/(2*Omega_max*F_d); 
end
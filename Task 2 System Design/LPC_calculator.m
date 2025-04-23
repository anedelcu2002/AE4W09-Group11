function [LPC] = LPC_calculator(E_y,D,D_orig,P_r,P_orig)
    % calculates the scaled diameter as a function of the rated power P_r in W, original power P_orig in W, and original diameter D_orig in m
    % then, calculates the LPC as a function of the calculated diameter D in m, the scaled diameter, and the yearly generated energy E_y in J

    D_scale = D_orig*((P_r/P_orig)^(1/2));         
    LPC = (0.7+(0.3*((D/D_scale)^2.6)))/E_y;      
end 
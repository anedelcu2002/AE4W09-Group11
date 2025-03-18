function [omeg_max] = tipspeed_max(U_r, lambda, D)
% function that calculates the maximal tip speed in m/s, using as input the rated wind speed and the tip speed ratio

<<<<<<< HEAD:Task 2 System Design/tipspeed_max.m
    %tip_v_max=lambda*U_r / (D/2);
    tip_v_max=lambda*U_r;
=======
    omeg_max= (lambda*U_r) / (D/2);

>>>>>>> Alex-module-breakdown:Task 2 System Design/omega_max.m
end
function [tip_v_max] = tipspeed_max(U_r, lambda, D)
% function that calculates the maximal tip speed in m/s, using as input the rated wind speed and the tip speed ratio

    %tip_v_max=lambda*U_r / (D/2);
    tip_v_max=lambda*U_r;
end
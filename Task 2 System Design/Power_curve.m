function [P_curve,U_r] = Power_curve(P_r,rho,U_ci,U_co,C_p,D,eff)
    % calculates the power curve and rated speed based on the rated power P_r in W, the air density rho in kg/m3, the cut-in and cut-out speeds
    % U_ci and U_co in m/s, the power coefficient C_p, the diameter D in m, the efficiency eff, and the speed time-series U_array in m/s

    R = 0.5*D;

    U_r = (P_r/(0.5*eff*rho*C_p*pi*R^2))^(1/3);

    P_curve = [];

    for i = 1:U_co
        if  i < U_ci
            P_curve(i) = 0;
        elseif i <= U_r
            P_curve(i) = 0.5*eff*rho*C_p*pi*(R^2)*(i^3);
        elseif i < U_co
            P_curve(i) = P_r;
        else
            P_curve(i) = 0;
        end
    end

end
function [P_curve,U_r,U_range] = Power_curve(P_r_LSS,rho,U_ci,U_co,C_p,D,eff)
    % calculates the power curve and rated speed based on the rated power P_r_LSS in W, the air density rho in kg/m3, the cut-in and cut-out speeds
    % U_ci and U_co in m/s from part 2, the power coefficient C_p from part 3, the diameter D in m from part 2, the efficiency eff_total

    R = 0.5*D;

    U_r = (P_r_LSS/(0.5*rho*C_p*pi*R^2))^(1/3);

    P_curve = [];
    U_range = 0:0.1:U_co+2;

    for i = 1:length(U_range)
        if  U_range(i) < U_ci
            P_curve(i) = 0;
        elseif U_range(i) <= U_r
            P_curve(i) = 0.5*eff*rho*C_p*pi*(R^2)*(U_range(i)^3);
        elseif U_range(i) < U_co
            P_curve(i) = P_r_LSS*eff;
        else
            P_curve(i) = 0;
        end
    end

end
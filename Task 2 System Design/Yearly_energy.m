function [E_y] = Yearly_energy(P_curve,U_ci,U_co,f_curve)
    %Anual Energy Yield, using as input the power curve, cut-in and cut-out speed in m/s and the Weibull wind distrubtion

    E_y=0;
    for i=U_ci:U_co
        E_y=3600*24*365*P_curve(i).*f_curve(i);
    end

end
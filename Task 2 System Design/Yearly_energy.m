function [Ey_kWh] = Yearly_energy(P_curve,U_ci,U_co,f_curve)
%Anual Energy Yield, using as input the power curve, cut-in and cut-out speed in m/s and the Weibull wind distrubtion

Ey_kWh = 24*365*integral(P_curve*f_curve, U_ci, U_co);

end
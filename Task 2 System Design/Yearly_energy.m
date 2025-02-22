function [Ey_kWh] = Yearly_energy(P_curve,U_r,f_curve)
%Anual Energy Yield
deltat = 365*24*3600;

Ey = deltat*f_curve*P_curve'; %J
Ey_kWh = (Ey/1000)/3600; %kWh 

end
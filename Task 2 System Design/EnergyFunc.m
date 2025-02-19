function [Ey_kWh] = EnergyFunc(R,Uci,Uco,k,a_scale,Pr,rho,eff,Cp)
%Windspeeds to look at
U_list = linspace(0, 200, 500);

%Weibelcurve
f_list = [];

for i = linspace(1,length(U_list),length(U_list))
    f_list(i) = (k/a_scale)*((U_list(i)/a_scale)^(k-1))*exp(-(U_list(i)/a_scale)^k);
end

%Urated
Ur = (Pr/(0.5*eff*rho*Cp*pi*R*R))^(1/3);

%Power curve
P_list = [];

for i = linspace(1,length(U_list),length(U_list))
    if  U_list(i) < Uci
        P_list(i) = 0;
    elseif U_list(i) < Ur
        P_list(i) = 0.5*eff*rho*Cp*pi*(R^2)*(U_list(i)^3);
    elseif U_list(i) < Uco
        P_list(i) = 0.5*eff*rho*Cp*pi*(R^2)*(Ur^3);
    else
        P_list(i) = 0;
    end
end

%Anual Energy Yield
deltat = 365*24*3600;

Ey = deltat*f_list*P_list'; %J
Ey_kWh = (Ey/1000)/3600; %kWh 
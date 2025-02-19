function [P_curve,Ur] = PowerFunc(Pr,rho,Uci,Uco,Cp,D,eff,U_range)
R = 0.5*D;

%Urated
Ur = (Pr/(0.5*eff*rho*Cp*pi*R*R))^(1/3);

%Power curve
P_curve = [];

for i = 1:length(U_range)
    if  U_range(i) < Uci
        P_curve(i) = 0;
    elseif U_range(i) < Ur
        P_curve(i) = 0.5*eff*rho*Cp*pi*(R^2)*(U_range(i)^3);
    elseif U_range(i) < Uco
        P_curve(i) = 0.5*eff*rho*Cp*pi*(R^2)*(Ur^3);
    else
        P_curve(i) = 0;
    end
end

end
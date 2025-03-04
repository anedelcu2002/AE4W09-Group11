function [Q] = torque(P_curve, U_co, D, lambda)
% Function to obtain the torque Q in Nm, for any wind speed U in m/s, also using as inputs the power curve,
% the rotor diameter D in m and the tip speed ratio lambda

    omega=[];
    for i=1:U_co
        omega(i)=lambda*i/(D/2);
    end
    Q = P_curve/omega;

end
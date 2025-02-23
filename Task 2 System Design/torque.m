function [Q] = torque (P_curve, U, D, lambda)
% Function to obtain the torque Q in Nm, for any wind speed U in m/s, also using as inputs the power curve,
% the rotor diameter D in m and the tip speed ratio lambda

    omega = lambda * U / (D/2);
    Q = P_curve(U)/omega;

end
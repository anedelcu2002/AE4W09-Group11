function BulgAir = steady_state_control_design(BulgAir,pitch,Cp_max,lambda)
% This is so we can run the steady-state operating points for all wind
% speeds.

% We need to define the following parameters:
% BulgAir.Control.WindSpeed
%     Cutin: cut-in wind speed in m/s
%     Cutout: cut-out wind speed in m/s
% BulgAir.Control.Pitch:
%     Fine: fine pitch in deg
% BulgAir.Control.Torque:
%     Demanded: Rated generator torque in Nm
%     Limit: Don't know
%     Slewrate: Don't know in Nm/s
%     OptGain: Kw^2 gain in Nm/(rad/s)^2
%     SpeedA: startup generator speed in rpm (for cut-in)
%     SpeedB: start of region 2 in generator rpm
%     SpeedB2: end of region 2 in generator rpm
%     SpeedC: rated speed in generator rpm
%     Min: Don't know in Nm
%     LowPassCutOffFreq: Don't know in rad/s

% TODO: Double-check these values.
BulgAir.Control.WindSpeed.Cutin = 3;
BulgAir.Control.WindSpeed.Cutout = 25;

BulgAir.Control.Pitch.Fine = pitch;


% Generator power goal.
Pg_W = 3.5e6;

% Optimal gain.
G = BulgAir.Drivetrain.Gearbox.Ratio;
rho_kgpm3 = BulgAir.AirDensity;
R_m = BulgAir.Blade.Radius(end);
drivetrain_efficiency = BulgAir.Drivetrain.Generator.Efficiency * BulgAir.Drivetrain.Gearbox.Efficiency;
K = rho_kgpm3 * pi * R_m^5 * Cp_max * drivetrain_efficiency / (2 * G^3 * lambda^3);

% generator_torque = K * generator_speed^2.
% generator_power =  K * generator_speed^3 * generator_efficiency.
omega_g_radps = (Pg_W / K / BulgAir.Drivetrain.Generator.Efficiency)^(1/3);
T_g_Nm = Pg_W / omega_g_radps / BulgAir.Drivetrain.Generator.Efficiency;
omega_g_rpm = omega_g_radps * 60 / (2*pi);
SpeedC = omega_g_rpm;

omega_r_cutin_radps = lambda * BulgAir.Control.WindSpeed.Cutin / R_m;
omega_r_cutin_rpm = omega_r_cutin_radps * 60 / (2*pi);
omega_g_cutin_rpm = omega_r_cutin_rpm * G;
SpeedB = omega_g_cutin_rpm;


% For now, I assume no torque limit, so we don't have a region 2.5.
Limit = T_g_Nm * 1.25;
SpeedB2 = SpeedC*0.95;

% Also assume no region 1.5.
SpeedA = SpeedB*0.75;
Min = 200;

% Assign them to the Control.Torque struct.
BulgAir.Control.Torque.Demanded = T_g_Nm;
BulgAir.Control.Torque.Limit = Limit;
BulgAir.Control.Torque.OptGain = K;
BulgAir.Control.Torque.SpeedA = SpeedA;
BulgAir.Control.Torque.SpeedB = SpeedB;
BulgAir.Control.Torque.SpeedB2 = SpeedB2;
BulgAir.Control.Torque.SpeedC = SpeedC;
BulgAir.Control.Torque.Min = Min;

end

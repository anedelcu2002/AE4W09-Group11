%% Reset and load data.
close all; clearvars; clc;

% Part 1 of tak 5 is to complete the turbine design in FASTTool. This
% includes the Torque control and pitch control (but not the detailed
% control design) and the tower design.

% Let's load our current design (just the rotor)
BulgAir = load("BulgAirChordTwist.mat");

% Now load BulgAirChordTwist in the FASTTool and run 'Steady operating
% curves' for different pitch angles but not for different wind speeds
% (because we're only interested in the below-rated response, so the
% Cp-lambda curve). Now load the results from that analysis.
results = load("BulgAir_CpLambdaPitch.mat");


%% Let's find the optimal operating condition.
% This is at max Cp.
[Cp_max, Cp_max_index] = max(results.Rotor_cP, [], "all");
[row, col] = ind2sub(size(results.Rotor_cP), Cp_max_index);
lambda = results.Rotor_Lamda(col);
pitch = results.Rotor_Pitch(row);


% Clean the workspace a bit.
clearvars row col Cp_max_index

%% Design input from previous tasks.
BulgAir.Drivetrain.Generator.Efficiency = 0.975;
BulgAir.Drivetrain.Gearbox.Ratio = 0.95;
BulgAir.Drivetrain.Gearbox.Efficiency = 0.95;



%% Now we need to define some control parameters
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
%     SpeedB: start of region 2 in rpm
%     SpeedB2: end of region 2 in rpm
%     SpeedC: rated speed in rpm
%     Min: Don't know in Nm
%     LowPassCutOffFreq: Don't know in rad/s

% TODO: Double-check these values.
BulgAir.Control.WindSpeed.Cutin = 3;
BulgAir.Control.WindSpeed.Cutout = 25;

BulgAir.Control.Pitch.Fine = pitch;

% Let's start over.

% Generator power goal.
Pg_W = 3.5e6;


% Optimal gain.
G = BulgAir.Drivetrain.Gearbox.Ratio;
rho_kgpm3 = BulgAir.AirDensity;
R_m = BulgAir.Blade.Radius(end);
drivetrain_efficiency = BulgAir.Drivetrain.Generator.Efficiency * BulgAir.Drivetrain.Gearbox.Efficiency;

K = rho_kgpm3 * pi * R_m^5 * Cp_max * drivetrain_efficiency / (2 * G^3 * lambda^3);
% generator_torque = K * generator_speed^2.
% generator_power =  K * generator_speed^3.
omega_g_radps = (Pg_W / K)^(1/3);
T_g_Nm = Pg_W / omega_g_radps;
omega_r_radps = omega_g_radps / G;
omega_r_rpm = omega_r_radps * 60 / (2*pi);
SpeedC = omega_r_rpm;

omega_r_cutin_radps = lambda * BulgAir.Control.WindSpeed.Cutin / R_m;
omega_r_cutin_rpm = omega_r_cutin_radps * 60 / (2*pi);
SpeedB = omega_r_cutin_rpm;


% %% OLD.
% % The torque control parameters take a bit of calculation.
% Pr_W = 3.5e6;
% drivetrain_efficiency = BulgAir.Drivetrain.Generator.Efficiency * BulgAir.Drivetrain.Gearbox.Efficiency;
% Pg_W = 3.5e6 / drivetrain_efficiency;
% rho_kgpm3 = BulgAir.AirDensity;
% R_m = BulgAir.Blade.Radius(end);
% A_m2 = pi * R_m^2;
% Ur_mps = (Pg_W / (0.5 * rho_kgpm3 * A_m2 * Cp_max))^(1/3);
% 
% % lambda = omega_r_radps * R / Ur_mps, so:
% omega_r_radps = lambda * Ur_mps / R_m;
% omega_r_rpm = omega_r_radps * 60 / (2*pi);
% SpeedC = omega_r_rpm;
% omega_r_cutin_radps = lambda * BulgAir.Control.WindSpeed.Cutin / R_m;
% omega_r_cutin_rpm = omega_r_cutin_radps * 60 / (2*pi);
% SpeedB = omega_r_cutin_rpm;
% 
% % Optimal gain.
% G = BulgAir.Drivetrain.Gearbox.Ratio;
% K = rho_kgpm3 * pi * R_m^5 * Cp_max * drivetrain_efficiency / (2 * G^3 * lambda^3);
% T_r_Nm = K * omega_r_radps^2;
% fprintf('Rated power goal %g and actual %g.\n', Pr_W, T_r_Nm * omega_r_radps * drivetrain_efficiency)
% 
% % power = torque * generator_speed * drivetrain_efficiency
% % power = K * generator_speed^3 * drivetrain_efficiency
% omega_g_radps = (Pr_W / K / drivetrain_efficiency)^(1/3);
% omega_r_radps = omega_g_radps / G;
% fprintf('Rated power goal %g and actual %g.\n', Pg_W, T_r_Nm * omega_r_radps)



% For now, I assume no torque limit, so we don't have a region 2.5.
Limit = Inf;
SpeedB2 = SpeedC*0.98;

% Also assume no region 1.5.
SpeedA = SpeedB*0.5;
Min = 0;

% Assign them to the Control.Torque struct.
BulgAir.Control.Torque.Demanded = T_g_Nm;
BulgAir.Control.Torque.Limit = Limit;
BulgAir.Control.Torque.OptGain = K;
BulgAir.Control.Torque.SpeedA = SpeedA;
BulgAir.Control.Torque.SpeedB = SpeedB;
BulgAir.Control.Torque.SpeedB2 = SpeedB2;
BulgAir.Control.Torque.SpeedC = SpeedC;
BulgAir.Control.Torque.Min = Min;


%% Add the tower design.
%TODO


%% Save the turbine.
save("BulgAirComplete.mat", "-struct", "BulgAir")

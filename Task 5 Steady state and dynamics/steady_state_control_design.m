function BulgAir = steady_state_control_design(BulgAir,NREL5MW,pitch,Cp_max,lambda)
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

BulgAir.Control.WindSpeed.Cutin = 2;
BulgAir.Control.WindSpeed.Cutout = 25;

BulgAir.Control.Pitch.Fine = pitch;

% Load some things into variables we'll need more often.
G = BulgAir.Drivetrain.Gearbox.Ratio;
n_generator = BulgAir.Drivetrain.Generator.Efficiency;
n_gearbox = BulgAir.Drivetrain.Gearbox.Efficiency;
n_drivetrain = n_generator * n_gearbox;
rho_kgpm3 = BulgAir.AirDensity;
R_m = BulgAir.Blade.Radius(end);

% SpeedC (rated generator speed) can be determined either by the maximum
% tip speed or just by intersection of the Kw^2 curve at rated power.

% Tip speed = omega * R
max_tip_speed_mps = 100;
max_rotor_speed_radps = max_tip_speed_mps / BulgAir.Blade.Radius(end);
speed_limit_radps = max_rotor_speed_radps * G;
SpeedC_tipspeed = radps2rpm(speed_limit_radps);

% Optimal gain.
K = rho_kgpm3 * pi * R_m^5 * Cp_max * n_drivetrain / (2 * G^3 * lambda^3);

% generator_power = Kw^2 * w * n_generator.
Pg_W = 3.5e6;
rated_generator_speed_radps = (Pg_W / K / n_generator)^(1/3);
SpeedC_power = radps2rpm(rated_generator_speed_radps);
SpeedC = min(SpeedC_tipspeed, SpeedC_power);


% The demanded torque is now determined by the generator power.
% generator_power = speed * torque * generator_efficiency.
rated_torque_Nm = Pg_W / rpm2radps(SpeedC) / n_generator;
torque_limit = NREL5MW.Control.Torque.Limit / NREL5MW.Control.Torque.Demanded * rated_torque_Nm;


% The torque limit and B2 point, use the same ratio as what was used for
% the NREL5MW turbine.
% Since our Turbine does not have a region 2.5, we just set B2 very close
% to C.
SpeedB2 = 0.99*SpeedC;


% SpeedA is the cut-in speed.
omega_r_cutin_radps = lambda * BulgAir.Control.WindSpeed.Cutin / R_m;
omega_g_cutin_rpm = radps2rpm(omega_r_cutin_radps) * G;
SpeedA = omega_g_cutin_rpm;

% SpeedB is then define with the ratio of A-B of the NREL5MW turbine.
SpeedB = SpeedA * NREL5MW.Control.Torque.SpeedB / NREL5MW.Control.Torque.SpeedA;

% I don't know what to do with the minimum torque. Let's scale it with the
% turbine size.
Min = NREL5MW.Control.Torque.Min * BulgAir.Blade.Radius(end) / NREL5MW.Blade.Radius(end);

% Assign them to the Control.Torque struct.
BulgAir.Control.Torque.Demanded = rated_torque_Nm;
BulgAir.Control.Torque.Limit = torque_limit;
BulgAir.Control.Torque.OptGain = K;
BulgAir.Control.Torque.SpeedA = SpeedA;
BulgAir.Control.Torque.SpeedB = SpeedB;
BulgAir.Control.Torque.SpeedB2 = SpeedB2;
BulgAir.Control.Torque.SpeedC = SpeedC;
BulgAir.Control.Torque.Min = Min;




% Let's also make a plot.
N = 1000;
generator_speed_rpm = linspace(0, SpeedC*1.2, N);
generator_speed_radps = rpm2radps(generator_speed_rpm);
optimal_torque_Nm = K * generator_speed_radps.^2;

temp = linspace(SpeedB, SpeedB2, N);
generator_speed_op_rpm = [0, SpeedA, temp, SpeedC];
generator_torque_op_rpm = [0, 0, K * rpm2radps(temp).^2, rated_torque_Nm];
speedpoints = [SpeedA, SpeedB, SpeedB2, SpeedC];
torquepoints = [0, K * rpm2radps(SpeedB).^2, K * rpm2radps(SpeedB2).^2, rated_torque_Nm];


figure; grid on; hold on;
plot(generator_speed_rpm, optimal_torque_Nm, '--');
plot(generator_speed_op_rpm, generator_torque_op_rpm)
xline(SpeedC_tipspeed)
plot(speedpoints, torquepoints, 'kx')

legend('Optimal', 'Variable-speed controller', 'tip speed limit')
xlabel('Generator speed (rpm)')
ylabel('Generator torque (Nm)')



end


function rpm = radps2rpm(radps)
rpm = radps * 60 / (2*pi);
end

function radps = rpm2radps(rpm)
radps = rpm * 2*pi / 60;
end

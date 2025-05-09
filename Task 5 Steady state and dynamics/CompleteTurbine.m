%% Reset and load data.
close all; clearvars; clc;

% Part 1 of tak 5 is to complete the turbine design in FASTTool. This
% includes the Torque control and pitch control (but not the detailed
% control design) and the tower design.

% Let's load our current design (just the rotor)
BulgAir = load("../Task 3 Rotor Design/BulgAirChordTwist.mat");
NREL5MW = load("../FastTool/NREL5MW.mat");

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
BulgAir.Drivetrain.Gearbox.Ratio = 100;
BulgAir.Drivetrain.Gearbox.Efficiency = 0.95;
BulgAir.Drivetrain.Generator.HSSInertia = 683;
BulgAir.Nacelle.Hub.Overhang = 5.675;

% TODO: Blade stiffness.
% Should be good to go since Alex added it in Task 3


%% Now we need to define some control parameters
BulgAir = steady_state_control_design(BulgAir,NREL5MW,pitch,Cp_max,lambda);
exportgraphics(gcf, 'generator_torque_curve.pdf')



%% Add the tower design.
BulgAir.Tower.Height = linspace(0, 120, 11);
BulgAir.Tower.Diameter = -0.02185*BulgAir.Tower.Height + 7.012;
BulgAir.Tower.WallThickness = -0.000104167*BulgAir.Tower.Height + 0.0409;
BulgAir.Tower.HubHeight = BulgAir.Tower.Height(end);


%% Save the turbine.
save("BulgAirComplete.mat", "-struct", "BulgAir")


%% Run the steady-state analysis on the FAST tool.
addpath('..\FASTTool\subfunctions');
% ModalAnalysis(BulgAir) TODO

%% Now make a plot
close all; clearvars; clc;
load('BulgAirComplete_steady_state_results.mat')

figure; hold on; grid on;
plot(WindSpeed, ElectricalPower./1e6)

ylabel('Electrical power (MW)')
xlabel('Wind speed (m/s)')

exportgraphics(gcf, 'power_curve.pdf')

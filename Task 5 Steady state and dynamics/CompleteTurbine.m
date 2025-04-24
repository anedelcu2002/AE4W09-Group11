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
BulgAir = steady_state_control_design(BulgAir,pitch,Cp_max,lambda);


%% Add the tower design.
%TODO


%% Save the turbine.
save("BulgAirComplete.mat", "-struct", "BulgAir")

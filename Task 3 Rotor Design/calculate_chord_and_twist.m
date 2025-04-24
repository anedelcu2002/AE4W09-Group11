%% Reset and load the Airfoil data
%close all; clearvars; clc;

% BulgAir includes the data of our new airfoils in the airfoil struct and
% defines the airfoil distribution in the blade in the blade struct.
BulgAir = load("BulgAir.mat").BulgAir;

% Also load the NREL5MW turbine for scaling.
NREL5MW = load("..\FASTTool\NREL5MW.mat");


%% Use the same indices.
% Before we start, make sure that all the r/R indices are the same and
% let's stick to the convention of the NREL5MW definition.
r_n = NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end);
r_n_bulgAir = BulgAir.Blade.Radius / BulgAir.Blade.Radius(end);
BulgAir.Blade.Radius = interp1(r_n_bulgAir, BulgAir.Blade.Radius, r_n);
BulgAir.Blade.NFoil = interp1(r_n_bulgAir, BulgAir.Blade.NFoil, r_n, 'nearest');
r_n_bulgAirNew = BulgAir.Blade.Radius / BulgAir.Blade.Radius(end);

% To see this effect:
figure; hold on;
plot(r_n, 'x')
plot(r_n_bulgAir, 'x')
plot(r_n_bulgAirNew, '+')
xlabel('index')
ylabel('r/R')
legend('NREL5MW', 'old', 'new')


%% Let's run our analytic and scaling functions.
lambda = 8;
lambda_old = 7;

R = 143/2;
R_old = 126/2;

analyticBlade = analytic_chord_and_twist_for_blade(BulgAir.Blade, BulgAir.Airfoil, lambda);
scaledBlade = scale_chord_and_twist_for_blade(NREL5MW.Blade, lambda_old, lambda, R);


%% Let's make some plots.
figure;
tiledlayout(2, 1)

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Twist);
plot(analyticBlade.Radius / analyticBlade.Radius(end), analyticBlade.Twist);
plot(scaledBlade.Radius / scaledBlade.Radius(end), scaledBlade.Twist);
grid;
ylabel('Twist (deg)')
xlabel('r/R (-)')

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Chord);
plot(analyticBlade.Radius / analyticBlade.Radius(end), analyticBlade.Chord);
plot(scaledBlade.Radius / scaledBlade.Radius(end), scaledBlade.Chord);
grid;
ylabel('Chord (m)')
xlabel('r/R (-)')


legend('NREL 5 MW', 'analytic 3.5 MW', 'scaled 3.5 MW');


%% Let's now make our final twist and chord distribution.
r_n = BulgAir.Blade.Radius / BulgAir.Blade.Radius(end);  % r/R

% For the twist, we want a smooth twist distribution that fits the optimal
% aerodynamic solution as well as possible. Near the root the optimal twist
% is undefined so we ignore it.
tempTwist = analyticBlade.Twist;
noTwist = isnan(tempTwist);
tempTwist = tempTwist(~noTwist);

% Fit a polynomial to this data.
twist_p = polyfit(r_n(~noTwist), tempTwist, 5);
twist = polyval(twist_p, r_n);

% Set the noTwist region so that it smoothly switches to the twisted
% region.
twist(noTwist) = twist(find(~noTwist, 1));


% For the chord, we want our scaled scaled solution near the root because
% structural effects dominate. Near the tip, we want our chord to be closer
% to the analytic solution for optimal aerodynamic performance. 
w = zeros(size(r_n));
w(r_n > 0.6) = 1;
w(r_n > 0.95) = 0;

tempChord = analyticBlade.Chord;
tempChord(isinf(tempChord)) = 0;

chord_weighted = w .* tempChord + (1 - w) .* scaledBlade.Chord;

chord_p = polyfit(r_n, chord_weighted, 5);
chord = polyval(chord_p, r_n);



%% Let's plot that.
figure;
tiledlayout(2, 1)

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Twist);
plot(analyticBlade.Radius / analyticBlade.Radius(end), analyticBlade.Twist);
plot(scaledBlade.Radius / scaledBlade.Radius(end), scaledBlade.Twist);
% plot(r_n, twist_weighted);
plot(r_n, twist, 'linewidth', 2);
grid;
ylabel('Twist (deg)')
xlabel('r/R (-)')

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Chord);
plot(analyticBlade.Radius / analyticBlade.Radius(end), analyticBlade.Chord);
plot(scaledBlade.Radius / scaledBlade.Radius(end), scaledBlade.Chord);
% plot(r_n, chord_weighted)
plot(r_n, chord, 'linewidth', 2)
grid;
ylabel('Chord (m)')
xlabel('r/R (-)')


legend('NREL 5 MW', 'analytic 3.5 MW', 'scaled 3.5 MW', 'fitted 3.5 MW');


%% And finally, let's save the results.
BulgAir.Blade.Twist = twist;
BulgAir.Blade.Chord = chord;

% Fix airfoil naming
BulgAir.Airfoil.Name = cellfun(@(x) char(x), BulgAir.Airfoil.Name, 'UniformOutput', false);

save("BulgAirChordTwist.mat", "-struct", "BulgAir")


%% Do the rotor analysis with these parameters.
clearvars;

%% FASTTool analysis.
% Now load BulgAirChordTwist in the FASTTool and run 'Steady operating
% curves' for different pitch angles but not for different wind speeds
% (because we're only interested in the below-rated response, so the
% Cp-lambda curve).
% Let's save those results
save("BulgAir_CpLambdaPitch.mat")


%% Now make the plot for the rotor analysis.
% We just need to add the NREL 5 MW design to it.
% I hacked this a bit in the command window to produce 'cP-lamda-pitch'.

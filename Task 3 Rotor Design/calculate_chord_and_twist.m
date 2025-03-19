%% Reset and load the Airfoil data
close all; clearvars; clc;

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
R_old = 126;

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
% We want our scaled solution near the root because structural effects
% dominate. Near the tip, we want our twist and chord to be closer to the
% analytic solution for optimal aerodynamic performance. Finally, at the
% tip we want to follow the scaled solution because of tip unloading (good
% for noise I think).
r_n = BulgAir.Blade.Radius / BulgAir.Blade.Radius(end);  % r/R
analytic_weighting = zeros(size(r_n));
analytic_weighting(r_n > 0.6) = 1;
analytic_weighting(r_n > 0.95) = 0;

% Before we do the weighting, we need to replace the nans and infs in the 
% analytic calculation.
analyticBladeTwist = analyticBlade.Twist;
analyticBladeChord = analyticBlade.Chord;
analyticBladeTwist(isnan(analyticBladeTwist)) = 0;
analyticBladeChord(isinf(analyticBladeChord)) = 0;

twist_weighted = analytic_weighting .* analyticBladeTwist + (1 - analytic_weighting) .* scaledBlade.Twist;
chord_weighted = analytic_weighting .* analyticBladeChord + (1 - analytic_weighting) .* scaledBlade.Chord;

% Now fit to that data.
twist_p = polyfit(r_n, twist_weighted, 5);
chord_p = polyfit(r_n, chord_weighted, 5);

% And apply those coefficients.
twist = polyval(twist_p, r_n);
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

save("BulgAirChordTwist.mat", "-struct", "BulgAir")
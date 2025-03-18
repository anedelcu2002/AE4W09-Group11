%% Reset and load the NREL 5 MW turbine.
close all; clearvars; clc;

NREL5MW = load("FASTTool\NREL5MW.mat");


%% Run the functions.
lambda = 7;  % TODO: get value for NREL5MW turbine.
AnalyticBlade5MW = analytic_chord_and_twist_for_blade(NREL5MW.Blade, NREL5MW.Airfoil, lambda);
a=5;

%% Plot
figure;
tiledlayout(2, 1)

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Twist);
plot(AnalyticBlade5MW.Radius / AnalyticBlade5MW.Radius(end), AnalyticBlade5MW.Twist);
grid;
ylabel('Twist (deg)')
xlabel('r/R (-)')

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Chord);
plot(AnalyticBlade5MW.Radius / AnalyticBlade5MW.Radius(end), AnalyticBlade5MW.Chord);
grid;
ylabel('Chord (m)')
xlabel('r/R (-)')


legend('NREL 5 MW', 'analytic 5 MW');


%% Now apply scaling laws.
lambda_new = lambda * 1.1;  % TODO
R_new = NREL5MW.Blade.Radius(end) * 1.4;    % TODO

ScaledBlade = scale_chord_and_twist_for_blade(NREL5MW.Blade, lambda, lambda_new, R_new);


%% Plot.
figure;
tiledlayout(2, 1)

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Twist);
plot(AnalyticBlade5MW.Radius / AnalyticBlade5MW.Radius(end), AnalyticBlade5MW.Twist);
plot(ScaledBlade.Radius / ScaledBlade.Radius(end), ScaledBlade.Twist);
grid;
ylabel('Twist (deg)')
xlabel('r/R (-)')

nexttile; hold on;
plot(NREL5MW.Blade.Radius / NREL5MW.Blade.Radius(end), NREL5MW.Blade.Chord);
plot(AnalyticBlade5MW.Radius / AnalyticBlade5MW.Radius(end), AnalyticBlade5MW.Chord);
plot(ScaledBlade.Radius / ScaledBlade.Radius(end), ScaledBlade.Chord);
grid;
ylabel('Chord (m)')
xlabel('r/R (-)')


legend('NREL 5 MW', 'analytic 5 MW', 'scaled using lambda and radius');


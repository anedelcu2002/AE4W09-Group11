%% to-do:
% Feedback from the session
% - Look at flap and edgewise moment stesses separately for fatigue.


% determine wind speeds for each simulation - alex


% calculate damage equivalent load for fatigue - jesse
% check simulation duration for fatigue - jesse -> 60 minutes without
% transients per wind speed (page 50)


% plot some nice spectra - victor


% calculate thickness factor and recheck design with task 3, task 5 - victor
clc
clear
close all

function cumulativeDamageStemPlot(ni,Nfi, blade_nr) % stolen from https://nl.mathworks.com/help/signal/ug/practical-introduction-to-fatigue-analysis-using-rainflow-counting.html
    figure
    L = length(ni);
    damage = sum(ni./Nfi);
    stem(0,NaN,"Color",[0 1 0])
    title(["Cumulative Damage from Palmgren-Miner Rule, blade ", num2str(blade_nr)])
    xlabel("Cycle $i$","Interpreter","latex")
    ylabel("Cum. damage $D_{i} = \sum_{j=1}^{i}n_{j}/N_{f,j}$","Interpreter","latex")
    set(gca,"XLim",[0 L],"YLim",[0 damage])
    grid("on")
    iter = unique([1:round(L/100):L,L]);
    hold(gca,"on")
    for i = iter
        cdi = sum(ni(1:i)./Nfi(1:i)); % cumulative damge upto cycle i
        plt = stem(i,cdi,"filled");
        setStemColor(plt,cdi,0.95)
    end
end

function setStemColor(hplt,cumulativeDamage,gamma)
c = lines(5);
c = c([2,3,5],:);
if (cumulativeDamage > 1)
    color = c(1,:);
else
    if (cumulativeDamage > gamma)
        c1 = c(1,:);
        c2 = c(2,:);
    else
        c1 = c(3,:);
        c2 = c(2,:);
    end
    color = zeros(1,3);    
    for i = 1:3
        color(i) = c1(i)+(c2(i)-c1(i))*cumulativeDamage;
    end
end
hplt.Color = color;
end

function stress=calculate_stress(moment_1, moment_2, thickness, EI_1, EI_2, E)
    moment_sum=(moment_1.^2.+moment_2.^2).^0.5;
    distance_1=thickness.*(moment_1./moment_sum);
    distance_2=thickness.*(moment_2./moment_sum);
    stress=moment_1.*distance_1./EI_1.*E+moment_2.*distance_2./EI_2.*E;
end

function [D, DEL] = rainflow_counting(stress_timeseries, blade_nr, plot_b) % stolen from brightspace
    gam_m=1; % safety factor for material properties
    gam_f=1.2; % safety factor for loads
    gam_n=1.15; % safety factor for severity of effect

    S_s= gam_m.*gam_f.*gam_n.*stress_timeseries;

    % Perform rainflow count, blade 1
    [c,hist,edges,rmm,idx] = rainflow(S_s);

    % Plot histogram
    if plot_b
        figure;
        histogram('BinEdges',edges','BinCounts',sum(hist,2))
        title(['Rainflow counting, blade ', num2str(blade_nr)])
        xlabel('Stress Range')
        ylabel('Cycle Counts')
    end

    % Set S-N curve and use Minor's rule to determine damage
    %   by summing damages from each individual cycle or half cycle
    m=9; % slope S-N curve N*S^m=K (with S stress range)
    UCS=600*10^6; % ultimate compression strength (Nm2);
    K=(2*UCS)^m; % S-N in stress ranges, so factor 2 for amplitudes

    range= c(:,2);
    count= c(:,1);
    D=1/K*sum(count.*(range.^m));
    
    if plot_b
        cumulativeDamageStemPlot(count,K./(range.^m), blade_nr);
    end

    % Also calculate the damage equivalent load.
    Dt = 0.008;
    neq = length(S_s) / Dt;  % assume 1h of data.
    DEL = (sum(count .* (range.^m)) / neq).^(1/m);
end


%% initialize data
response_data=load('nrel_cert.mat');
response_data.Legend;

turbine_data=load('NREL5MW.mat');

E=14.7*10^9;
EI_edge=turbine_data.Blade.EIedge(1);
EI_flap=turbine_data.Blade.EIflap(1);
thickness=turbine_data.Blade.Thickness(1); %base root is a load-bearing cylinder

%% plot simulation outputs
plot_b=true;

if plot_b
    figure;
    plot(response_data.Time, response_data.OoPDefl1);
    hold on;
    plot(response_data.Time, response_data.OoPDefl2);
    hold on;
    plot(response_data.Time, response_data.OoPDefl3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('Out of plane tip deflection')
    ylabel('deflection (m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMEdg1);
    hold on;
    plot(response_data.Time, response_data.RootMEdg2);
    hold on;
    plot(response_data.Time, response_data.RootMEdg3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('Blade root moment, edge direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMFlp1);
    hold on;
    plot(response_data.Time, response_data.RootMFlp2);
    hold on;
    plot(response_data.Time, response_data.RootMFlp3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('Blade root moment, flap direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')
end

%% perform post-processing
analysis=0; % 1 for extreme load, 0 for fatigue

if analysis==1
    maximum_deflection=max([max(response_data.OoPDefl1), max(response_data.OoPDefl2), max(response_data.OoPDefl3)]); %maximum out of plane tip deflection in meters
    sf_m_defl=1.1; % safety factor for material properties
    sf_f_defl=1.35; % safety factor for loads
    sf_n_defl=1; % safety factor for severity of effect
    maximum_deflection=maximum_deflection*sf_m_defl*sf_f_defl*sf_n_defl; % maximum deflection in meters
    disp(['Maximum tip deflection in out-of-plane direction is ', num2str(maximum_deflection), ' meters'])

    root_stress1=calculate_stress(response_data.RootMFlp1*1000, response_data.RootMEdg1*1000, thickness/2, EI_flap, EI_edge, E);%blade 1 root stress timeseries
    root_stress2=calculate_stress(response_data.RootMFlp2*1000, response_data.RootMEdg2*1000, thickness/2, EI_flap, EI_edge, E);
    root_stress3=calculate_stress(response_data.RootMFlp3*1000, response_data.RootMEdg3*1000, thickness/2, EI_flap, EI_edge, E);

    if plot_b
        figure;
        plot(response_data.Time, root_stress1);
        hold on;
        plot(response_data.Time, root_stress2);
        hold on;
        plot(response_data.Time, root_stress3);
        legend('blade 1', 'blade 2', 'blade 3');
        title('Blade root stress amplitude')
        ylabel('stress (N/m2)')
        xlabel('time (s)')
    end

    gam_m=1.35; % safety factor for material properties
    gam_f=1.3; % safety factor for loads
    gam_n=1; % safety factor for severity of effect

    maximum_stress=gam_m*gam_f*gam_n*max([max(root_stress1), max(root_stress2), max(root_stress3)]); %maximum blade root stress after safety factors
    disp(['Maximum blade root stress is ', num2str(maximum_stress/10^6), ' MPa'])

elseif analysis==0
    % TODO: remove transient part of the simulation
    root_stress1=calculate_stress(response_data.RootMFlp1*1000, response_data.RootMEdg1*1000, thickness/2, EI_flap, EI_edge, E);%blade 1 root stress timeseries
    root_stress2=calculate_stress(response_data.RootMFlp2*1000, response_data.RootMEdg2*1000, thickness/2, EI_flap, EI_edge, E);
    root_stress3=calculate_stress(response_data.RootMFlp3*1000, response_data.RootMEdg3*1000, thickness/2, EI_flap, EI_edge, E);

    if plot_b
        figure;
        plot(response_data.Time, root_stress1);
        hold on;
        plot(response_data.Time, root_stress2);
        hold on;
        plot(response_data.Time, root_stress3);
        legend('blade 1', 'blade 2', 'blade 3');
        title('Blade root stress amplitude')
        ylabel('stress (N/m2)')
        xlabel('time (s)')
    end

    [D1, DEL1] = rainflow_counting(root_stress1, 1, plot_b);
    [D2, DEL2] = rainflow_counting(root_stress2, 2, plot_b);
    [D3, DEL3] = rainflow_counting(root_stress3, 3, plot_b);
    disp(['Cumulative fatigue damage on blade 1 is ', num2str(D1), ', with DEL ', num2str(DEL1)])
    disp(['Cumulative fatigue damage on blade 2 is ', num2str(D2), ', with DEL ', num2str(DEL2)])
    disp(['Cumulative fatigue damage on blade 3 is ', num2str(D3), ', with DEL ', num2str(DEL3)])

    fprintf("Average damage for the three blades is %.4g.\n", mean([D1, D2, D3]))

end


%% Lifetime fatigue calculation.
% Using equation H.3 from the IEC 61400-1.
function [lifetime_damage, accumulated_damage] = calc_lifetime_damage(lifetime_s, period_s, windspeeds, damages, wbl_scale, wbl_shape)
arguments
    lifetime_s double  % Lifetime of the turbine in seconds.
    period_s double    % Period over which the damages were calculated
    windspeeds double  % Wind speeds for which the damages were calculated
    damages double     % Damage for the wind speed
    wbl_scale double   % Weibull scale parameter 
    wbl_shape double   % Weibull shape parameter 
end

assert(length(damages) == length(windspeeds), "damages (%i) must be calculated for each windspeed (%i).", length(damages), length(windspeeds))

dw = diff(windspeeds);
assert(std(dw) == 0, "windspeeds must be evenly spaced.")
dw = dw(1);
windspeeds_lower_bounds = windspeeds - 0.5 * dw;
windspeeds_upper_bounds = windspeeds + 0.5 * dw;

windspeeds_probability = wblcdf(windspeeds_upper_bounds, wbl_scale, wbl_shape) - wblcdf(windspeeds_lower_bounds, wbl_scale, wbl_shape);

lifetime_damage = lifetime_s / period_s * sum(damages .* windspeeds_probability);
accumulated_damage = lifetime_s / period_s * cumsum(damages .* windspeeds_probability);

end

% Inputs to the lifetime damage calculation.
lifetime_s = seconds(years(20));
period_s = seconds(minutes(10));

windspeeds = 3:25;
% Hmm this is not so super nice, it'd be nicer to get the damages of all
% the simulations automatically.
damages = [1e-8, 1e-7, 1e-6, ...
    1e-6, 1e-6, 1e-6, 1e-6, 1e-6, ...
    1e-6, 1e-6, 1e-6, 1e-6, 1e-6, ...
    1e-6, 1e-6, 1e-7, 1e-6, 1e-6, ...
    1e-6, 1e-6, 1e-6, 1e-6, 1e-6];
    
% Weibull distribution at hub height (see bottom of MAIN.m in Task 2).
scale = 8.45387;
shape = 1.50374;

[lifetime_damage, accumulated_damage] = calc_lifetime_damage(lifetime_s, period_s, windspeeds, damages, scale, shape);
fprintf("Lifetime damage is %.3f (should be below 1).\n", lifetime_damage)

%% Make some nice plots.
% TODO: Look up what plots they want us to put in the report.
figure;
yyaxis left
plot(windspeeds, damages)
ylim([0, max(damages)*1.2])
set(gca, 'YScale', 'log')
ylabel('Damage (-)')

yyaxis right
all_windspeeds = linspace(0, 30, 1000);
plot(all_windspeeds, wblpdf(all_windspeeds, scale, shape))
ylabel('Wind speed probability (-)')

xlabel('Wind speed (m/s)')

% legend("Damage per wind speed", "Wind speed distribution")


%M=response_data.RootMFlp1; % blade root flap moment (kNm)
% Remove 60 s start-up transient
%Ind=find(response_data.Time<=60);
%M(Ind)=[];
% Set geometrical properties of the cross section and determine stress 
%I=1.37; % blade root inertia (m^4)
%y=1.5; % distance to neutral line (m)
%S=M*y/I;
% Set safety factors and apply these to stress
%gam_m=1; % safety factor for material properties
%gam_f=1; % safety factor for loads
%gam_n=1; % safety factor for severity of effect
%S_s = gam_m*gam_f*gam_n*S;

% Perform rainflow count
%[c,hist,edges,rmm,idx] = rainflow(S_s);

% Plot histogram
%histogram('BinEdges',edges','BinCounts',sum(hist,2))
%xlabel('Stress Range')
%ylabel('Cycle Counts')

% Set S-N curve and use Minor's rule to determine damage
%   by summing damages from each individual cycle or half cycle
%m=9; % slope S-N curve N*S^m=K (with S stress range)
%UCS=6e5; % ultimate compression strength (kNm);
%K=(2*UCS)^m; % S-N in stress ranges, so factor 2 for amplitudes
%range= c(:,2);
%count= c(:,1);
%D=1/K*sum(count.*(range.^m));

%% Plotting spectra for clean representation of results
% In-plane and out-of-plane blade root moment spectra and removing average
% offset (The average displays no dynamic behaviour, just the constant weight of the
% blades on the moment)
M_in_plane_1 = response_data.RootMEdg1 - mean(response_data.RootMEdg1);
M_out_plane_1 = response_data.RootMFlp1 - mean(response_data.RootMFlp1);

M_in_plane_2 = response_data.RootMEdg2 - mean(response_data.RootMEdg2);
M_out_plane_2 = response_data.RootMFlp2 - mean(response_data.RootMFlp2);

M_in_plane_3 = response_data.RootMEdg3 - mean(response_data.RootMEdg3);
M_out_plane_3 = response_data.RootMFlp3 - mean(response_data.RootMFlp3);

Fs = (length(M_in_plane_1)-1)/seconds(minutes(11)); %Sampling Frequency (I think its this)

% Append first minute of simulation
M_in_plane_1 = M_in_plane_1(60*Fs+1:end);
M_out_plane_1 = M_out_plane_1(60*Fs+1:end);

M_in_plane_2 = M_in_plane_2(60*Fs+1:end);
M_out_plane_2 = M_out_plane_2(60*Fs+1:end);

M_in_plane_3 = M_in_plane_3(60*Fs+1:end);
M_out_plane_3 = M_out_plane_3(60*Fs+1:end);

T = seconds(minutes(10)); %Total duration of sampling
N = Fs*T; %Number of samples
t = (0:N-1)/Fs; % Time that each sample is taken in seconds

% Fast Fourier Transform
Y_in_1 = fft(M_in_plane_1);
Y_out_1 = fft(M_out_plane_1);
Y_in_2 = fft(M_in_plane_2);
Y_out_2 = fft(M_out_plane_2);
Y_in_3 = fft(M_in_plane_3);
Y_out_3 = fft(M_out_plane_3);

% Compute the magnitude (only) of the Fourier values and normalize with N
P2_in_1 = abs(Y_in_1/N);
P2_out_1 = abs(Y_out_1/N);
P2_in_2 = abs(Y_in_2/N);
P2_out_2 = abs(Y_out_2/N);
P2_in_3 = abs(Y_in_3/N);
P2_out_3 = abs(Y_out_3/N);

% Convert to one-sided spectrum due to the complex-conjugate redundancy
P1_in_1 = P2_in_1(1:N/2+1);
P1_out_1 = P2_out_1(1:N/2+1);
P1_in_2 = P2_in_2(1:N/2+1);
P1_out_2 = P2_out_2(1:N/2+1);
P1_in_3 = P2_in_3(1:N/2+1);
P1_out_3 = P2_out_3(1:N/2+1);

% Multiply by 2 (except for DC and Nyquist)
P1_in_1(2:end-1) = 2*P1_in_1(2:end-1);
P1_out_1(2:end-1) = 2*P1_out_1(2:end-1);
P1_in_2(2:end-1) = 2*P1_in_2(2:end-1);
P1_out_2(2:end-1) = 2*P1_out_2(2:end-1);
P1_in_3(2:end-1) = 2*P1_in_3(2:end-1);
P1_out_3(2:end-1) = 2*P1_out_3(2:end-1);

% Frequency vector up to Nyquist Frequency
f = Fs*(0:(N/2))/N;

% Plot in-plane bending moment spectrum
figure;
hold on
plot(f, P1_in_1,"DisplayName","Blade 1")
plot(f, P1_in_2,"DisplayName","Blade 2")
plot(f, P1_in_3,"DisplayName","Blade 3")
title('Spectrum of In-Plane Bending Moment')
xlabel('Frequency (Hz)')
ylabel('Root Bending Moment (kNm)')
xlim([0 2])
legend

% Plot out-of-plane bending moment spectrum
figure;
hold on
plot(f, P1_out_1,"DisplayName","Blade 1")
plot(f, P1_out_2,"DisplayName","Blade 2")
plot(f, P1_out_3,"DisplayName","Blade 3")
title('Spectrum of Out-of-Plane Bending Moment')
xlabel('Frequency (Hz)')
ylabel('Root Bending Moment (kNm)')
xlim([0 2])
legend
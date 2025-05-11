%% to-do:
% figure out e-modulus
% calculate stress at each moment rather than maximum superimposed stresses
% add safety factors
% check whether to use ip/oop or edge/flap moments
% calculate fatigue damage and damage equivalent load for fatigue

response_data=load('nrel_cert.mat');
response_data.Legend;

turbine_data=load('NREL5MW.mat');



analysis='extreme_load';

if analysis=='extreme_load'


    figure;
    plot(response_data.Time, response_data.OoPDefl1);
    hold on;
    plot(response_data.Time, response_data.OoPDefl2);
    hold on;
    plot(response_data.Time, response_data.OoPDefl3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('out of plane tip deflection')
    ylabel('deflection (m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMEdg1);
    hold on;
    plot(response_data.Time, response_data.RootMEdg2);
    hold on;
    plot(response_data.Time, response_data.RootMEdg3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root moment edge direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    figure;
    plot(response_data.Time, response_data.RootMFlp1);
    hold on;
    plot(response_data.Time, response_data.RootMFlp2);
    hold on;
    plot(response_data.Time, response_data.RootMFlp3);
    legend('blade 1', 'blade 2', 'blade 3');
    title('blade root moment flap direction')
    ylabel('moment (kN m)')
    xlabel('time (s)')

    maximum_deflection=max([max(response_data.OoPDefl1), max(response_data.OoPDefl2), max(response_data.OoPDefl3)]); %maximum out of plane tip deflection in meters
    maximum_edge_moment=max([max(response_data.RootMEdg1), max(response_data.RootMEdg2), max(response_data.RootMEdg3)]); %maximum edge blade root moment in kNm
    maximum_flap_moment=max([max(response_data.RootMFlp1), max(response_data.RootMFlp2), max(response_data.RootMFlp3)]); %maximum edge blade root moment in kNm

    EI_edge=turbine_data.Blade.EIedge(1);
    EI_flap=turbine_data.Blade.EIflap(1);
    thickness=turbine_data.Blade.Thickness(1)/sqrt(2); %base root is a cylinder, assuming spar structure is square enclosed in circle

    maximum_stress=maximum_flap_moment*1000*thickness/2/EI_flap+maximum_edge_moment*1000*thickness/2/EI_edge; %maximum blade root stress in N/m2

    max_stress_mpa=maximum_stress/10^6*4.7*10^9;

elseif analysis=='fatigue_load'


    figure;
    plot(response_data.Time, response_data.OoPDefl1);
    hold on;
    plot(response_data.Time, response_data.OoPDefl2);
    hold on;
    plot(response_data.Time, response_data.OoPDefl3);
    legend('blade 1', 'blade 2', 'blade 3');
end
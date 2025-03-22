%% Structural scaling code, limit state analysis

% Assumptions:


%% Set up functions for stress and deflection calculation at every position

% make sure to unit test these functions, they are coded by alex so inherently wrong af!

function stress = stress_calculation(r, t, F_in_plane, F_out_of_plane, EI_in_plane, EI_out_of_plane)
    stations=length(r);
    stress=zeros(1,stations);
    for i=1:stations
        stress(i)=sum(transpose(F_in_plane(i:stations)).*(r(i:stations)-r(i)))*(t(i)/2)/EI_in_plane(i)+sum(transpose(F_out_of_plane(i:stations)).*...
            (r(i:stations)-r(i)))*(t(i)/2)/EI_out_of_plane(i);
    end
end

function [deflection_in_plane, deflection_out_of_plane] = deflection_calculation(r, F_in_plane, F_out_of_plane, EI_in_plane, EI_out_of_plane)
    EI_in_plane=transpose(EI_in_plane);
    EI_out_of_plane=transpose(EI_out_of_plane);
    shear_force_in_plane=cumtrapz(r, F_in_plane); % moment equivalent to EI*y'', only integrating once because shear force V is already calculated
    shear_force_out_of_plane=cumtrapz(r, F_out_of_plane);
    deflection_in_plane=cumtrapz(r, cumtrapz(r, shear_force_in_plane./EI_in_plane));
    deflection_out_of_plane=cumtrapz(r, cumtrapz(r, shear_force_out_of_plane./EI_out_of_plane));
end

%% Initialize parameters for BulgAir

calculate_chord_and_twist;
BEM;

station_array=BulgAir.Blade.Radius; % stations in meters
mass_array=BulgAir.Blade.Mass*(143/126)^3; % distributed mass in kilograms, scaled by R^3 from the original NREL 5MW turbine
in_plane_force_array=abs(FTang); % in-plane force in newtons
out_of_plane_force_array=abs(FAxial); % out-of-plane force in newtons
twist_array=Twist*2*pi/360; % twist array in radians
chord_array=Chord; % chord length in meters
flap_stiffness_array=BulgAir.Blade.EIflap; % flap stiffness in Nm2, oriented by structural twist angle, assumed equal to aerodynamic twist angle
edge_stiffness_array=BulgAir.Blade.EIedge;

%% Calculate airfoil thickness at each station

airfoil_thickness_array=BulgAir.Blade.NFoil;
for i=1:48
    if airfoil_thickness_array(i)==1 % Cylinder
        airfoil_thickness_array(i)=sqrt(2)/2; % for a cylinder, the spar cap is a square inscribed within the circle
    elseif airfoil_thickness_array(i)==2 % DU 97-W-300
        airfoil_thickness_array(i)=0.3;
    elseif airfoil_thickness_array(i)==3 % Blend 1
        airfoil_thickness_array(i)=0.2545; % for a blend, thickness ratio is average of the blended airfoils
    elseif airfoil_thickness_array(i)==4 % NACA 644
        airfoil_thickness_array(i)=0.209;
    elseif airfoil_thickness_array(i)==5 % Blend 2
        airfoil_thickness_array(i)=0.1945;
    elseif airfoil_thickness_array(i)==6 % NACA 643
        airfoil_thickness_array(i)=0.18;
    end
end
airfoil_thickness_array=airfoil_thickness_array.*Chord;

%% Calculate in-plane/out-of-plane stiffnesses and forces

in_plane_stiffness_array=flap_stiffness_array.*cos(abs(twist_array))+edge_stiffness_array.*sin(abs(twist_array));
out_of_plane_stiffness_array=edge_stiffness_array.*cos(abs(twist_array))+flap_stiffness_array.*sin(abs(twist_array));

in_plane_force_array=in_plane_force_array+transpose(mass_array.*9.81); % maximum load case scenario has the blade perpendicular to tower, distributed mass acts as bending load

%% Do the same for NREL 5MW scaled to 143m diameter

BEM_NREL;

station_array_NREL=NREL5MW.Blade.Radius*143/126; % stations in meters
mass_array_NREL=NREL5MW.Blade.Mass*(143/126)^3; % distributed mass in kilograms
in_plane_force_array_NREL=abs(FTang_NREL); % in-plane force in newtons
out_of_plane_force_array_NREL=abs(FAxial_NREL); % out-of-plane force in newtons
twist_array_NREL=NREL5MW.Blade.Twist*2*pi/360; % twist array in radians
chord_array_NREL=NREL5MW.Blade.Chord; % chord length in meters
flap_stiffness_array_NREL=NREL5MW.Blade.EIflap; 
edge_stiffness_array_NREL=NREL5MW.Blade.EIedge;
airfoil_thickness_array_NREL=NREL5MW.Blade.Thickness;

in_plane_stiffness_array_NREL=flap_stiffness_array_NREL.*cos(abs(twist_array))+edge_stiffness_array_NREL.*sin(abs(twist_array));
out_of_plane_stiffness_array_NREL=edge_stiffness_array_NREL.*cos(abs(twist_array))+flap_stiffness_array_NREL.*sin(abs(twist_array));

in_plane_force_array_NREL=in_plane_force_array_NREL+transpose(mass_array_NREL.*9.81);

% calculate transverse stress from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

maximum_stress_NREL=stress_calculation(station_array_NREL, airfoil_thickness_array_NREL, in_plane_force_array_NREL,... 
out_of_plane_force_array_NREL,in_plane_stiffness_array_NREL, out_of_plane_stiffness_array_NREL);

% calculate transverse deflection from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

[max_in_plane_deflection_NREL, max_out_of_plane_deflection_NREL]=deflection_calculation(station_array_NREL, in_plane_force_array_NREL,... 
out_of_plane_force_array_NREL, in_plane_stiffness_array_NREL, out_of_plane_stiffness_array_NREL);

%% Scaled turbine response calculation

% calculate transverse stress from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

maximum_stress=stress_calculation(station_array, airfoil_thickness_array, in_plane_force_array, out_of_plane_force_array,... 
in_plane_stiffness_array, out_of_plane_stiffness_array);

% calculate transverse deflection from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

[max_in_plane_deflection, max_out_of_plane_deflection]=deflection_calculation(station_array, in_plane_force_array, out_of_plane_force_array,... 
in_plane_stiffness_array, out_of_plane_stiffness_array);

%% Plot results
figure;
tiledlayout(3, 1)

nexttile; hold on;
plot(station_array, maximum_stress);
plot(station_array_NREL, maximum_stress_NREL);
grid;
ylabel('Maximum Stress/Elastic Modulus')
xlabel('Station')
title('Stress and deflection at every station')

nexttile; hold on;
plot(station_array, max_in_plane_deflection);
plot(station_array_NREL, max_in_plane_deflection_NREL);
grid;
ylabel('In-plane Deflection')
xlabel('Station')

nexttile; hold on;
plot(station_array, max_out_of_plane_deflection);
plot(station_array_NREL, max_out_of_plane_deflection_NREL);
grid;
ylabel('Out-of-plane Deflection')
xlabel('Station')

legend('BulgAir', 'Scaled NREL 5MW');


%% Maximum required thickness factor calculation

% comparison between response values for scaled NREL 5MW and designed turbine, and thickness factor scaling to stress ratio 

stress_factor=max(maximum_stress./maximum_stress_NREL)

in_plane_deflection_factor=max(max_in_plane_deflection/max_in_plane_deflection_NREL)

out_of_plane_deflection_factor=max(max_out_of_plane_deflection/max_out_of_plane_deflection_NREL)
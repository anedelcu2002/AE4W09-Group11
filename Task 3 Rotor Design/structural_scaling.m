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
    shear_force_in_plane=cumtrapz(r, cumtrapz(r, F_in_plane));
    shear_force_out_of_plane=cumtrapz(r, cumtrapz(r, F_out_of_plane)); % shear q equivalent to EI*y''
    deflection_in_plane=cumtrapz(r, cumtrapz(r, shear_force_in_plane./EI_in_plane));
    deflection_out_of_plane=cumtrapz(r, cumtrapz(r, shear_force_out_of_plane./EI_out_of_plane));
end


%% Initialize parameters for BulgAir

calculate_chord_and_twist;
BEM;

station_array=BulgAir.Blade.Radius; % stations in meters
mass_array=BulgAir.Blade.Mass; % distributed mass in kilograms
in_plane_force_array=abs(FTang); % in-plane force in newtons
out_of_plane_force_array=abs(FAxial); % out-of-plane force in newtons
twist_array=Twist*2*pi/180; % twist array in radians
chord_array=Chord; % chord length in meters
flap_stiffness_array=BulgAir.Blade.EIflap; % flap stiffness in Nm2, oriented by structural twist angle, assumed equal to aerodynamic twist angle
edge_stiffness_array=BulgAir.Blade.EIedge;

%% Calculate airfoil thickness at each station

airfoil_thickness_array=BulgAir.Blade.NFoil;
for i=airfoil_thickness_array
    if i==1 % Cylinder
        airfoil_thickness_array(i)=sqrt(2)/2; % for a cylinder, the spar cap is a square inscribed within the circle
    elseif i==2 % DU 97-W-300
        airfoil_thickness_array(i)=0.3;
    elseif i==3 % Blend 1
        airfoil_thickness_array(i)=0.2545; % for a blend, thickness ratio is average of the blended airfoils
    elseif i==4 % NACA 644
        airfoil_thickness_array(i)=0.209;
    elseif i==5 % Blend 2
        airfoil_thickness_array(i)=0.1945;
    elseif i==6 % NACA 643
        airfoil_thickness_array(i)=0.18;
    end
end
airfoil_thickness_array=airfoil_thickness_array.*Chord;

%% Calculate in-plane/out-of-plane stiffnesses and forces

in_plane_stiffness_array=flap_stiffness_array.*cos(abs(twist_array))+edge_stiffness_array.*sin(abs(twist_array));
out_of_plane_stiffness_array=edge_stiffness_array.*cos(abs(twist_array))+flap_stiffness_array.*sin(abs(twist_array));

in_plane_force_array=in_plane_force_array+transpose(mass_array.*9.81); % maximum load case scenario has the blade perpendicular to tower, distributed mass acts as bending load

%% Original turbine response calculation

% calculate transverse stress from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

% calculate transverse deflection from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

% calculate transverse natural frequency from aerodynamic loading (thrust, torque) and weight loading, worst case scenario


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
grid;
ylabel('Maximum Stress/Elastic Modulus')
xlabel('Station')
title('Stress and deflection at every station')

nexttile; hold on;
plot(station_array, max_in_plane_deflection);
grid;
ylabel('In-plane Deflection')
xlabel('Station')

nexttile; hold on;
plot(station_array, max_out_of_plane_deflection);
grid;
ylabel('Out-of-plane Deflection')
xlabel('Station')


%% Maximum required thickness factor calculation

% comparison between response values for scaled NREL 5MW and designed turbine, and thickness factor scaling to stress ratio 
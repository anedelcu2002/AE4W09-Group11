%% Structural scaling code, limit state analysis

% Assumptions:


%% Stress, deflection, and natural stiffness calculation functions

% make sure to unit test these functions, they are coded by alex so inherently wrong af!

function stress = stress_calculation(r, m, c, t, F_in_plane, F_out_of_plane, EI_in_plane, EI_out_of_plane)
    stations=length(r);
    stress=zeros(stations);
    for i=1:stations
        stress(i)=sum(F_in_plane(i:stations).*(r(i:stations)-r(i)))*(t/2)/EI_in_plane+sum(F_out_of_plane(i:stations).*...
            (r(i:stations)-r(i)))*(t/2)/EI_out_of_plane;
    end
end

function [deflection_in_plane, deflection_out_of_plane] = deflection_calculation(r, F_in_plane, F_out_of_plane, EI_in_plane, EI_out_of_plane)
    shear_force_in_plane=cumtrapz(r, cumtrapz(r, F_in_plane));
    shear_force_out_of_plane=cumtrapz(r, cumtrapz(r, F_out_of_plane)); % shear q equivalent to EI*y''
    deflection_in_plane=cumtrapz(r, cumtrapz(r, shear_force_in_plane/EI_in_plane));
    deflection_out_of_plane=cumtrapz(r, cumtrapz(r, shear_force_out_of_plane/EI_out_of_plane));
end


%% Parameter initialization

BEM % not too fond of doing this
MAIN

% need to initialize this data for original turbine and for our turbine
station_array=[];
mass_array=[];
thrust_array=CT*0.5*rho*U_rated^2*(pi*D^2/4); % divide by three
torque_array=CQ*0.5*rho*U_rated^3*(pi*D^2/4);
twist_array=[];
chord_array=[];
airfoil_thickness_array=[];
flap_stiffness_array=[]; % oriented by structural twist angle, assumed equal to aerodynamic twist angle
edge_stiffness_array=[];

% airfoil distribution 0-20 cylinder, 20-40, 40-70, 70-100

%% Original turbine response calculation

% calculate transverse stress from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

% calculate transverse deflection from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

% calculate transverse natural frequency from aerodynamic loading (thrust, torque) and weight loading, worst case scenario


%% Scaled turbine response calculation

% calculate transverse stress from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

% calculate transverse deflection from aerodynamic loading (thrust, torque) and weight loading, worst case scenario

% calculate transverse natural frequency from aerodynamic loading (thrust, torque) and weight loading, worst case scenario


%% Maximum required thickness factor calculation

% comparison between response values for scaled NREL 5MW and designed turbine, and thickness factor scaling to stress ratio 
%% Structural scaling code, limit state analysis

% Assumptions:


%% Stress, deflection, and natural stiffness calculation functions

function stress = stress_calculation(r, m, c, F_in_plane, F_out_of_plane, EI_in_plane, EI_out_of_plane)
    stations=length(r);
    stress=zeros(stations);
    for i=1:stations
        stress(i)=sum(F_in_plane(i:stations).*(r(i:stations)-r(i)))/EI_in_plane+sum(F_out_of_plane(i:stations).*...
            (r(i:stations)-r(i)))/EI_out_of_plane; % need to multiply by distance from centroid, which is airfoil thickness (=t/c * c)
    end
end


%% Parameter initialization

BEM % not too fond of doing this
MAIN

% need to initialize this data for original turbine and for our turbine
station_array=[];
mass_array=[];
thrust_array=CT*0.5*rho*U_rated^2*(pi*D^2/4);
torque_array=CQ*0.5*rho*U_rated^3*(pi*D^2/4);
twist_array=[];
flap_stiffness_array=[]; % oriented by structural twist angle, assumed equal to aerodynamic twist angle
edge_stiffness_array=[];
chord_array=[];

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
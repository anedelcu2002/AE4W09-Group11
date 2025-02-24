%% INPUT DATA
U_0=[]; % measured wind speed at altitude h_0 in m/s
z_0=1; % surface roughness during measurement in m
h_0=1; % altitude of wind speed measurement in m
P_rated=1; % selected rated power in W
rho=1.225; % air density at hub level in kg/m3
P_original=1; % original turbine rated power in W
D_original=1; % original turbine diameter in W
E_y_original=1; % original turbine yearly energy generation in J (may not be needed)
U_ci=1; % assumed cut-in speed in m/s
U_co=1; % assumed cut-out speed in m/s
c_p=0.5; % assumed power coefficient
D_0=1; % assumed initial diameter for iterator
alpha=0.5; % constant for power law wind speed profile
eff=0.9; % assumed turbine efficiency
lambda_design=8; % assumed design tip speed ratio
max_tip_speed_limit=100; % upper limit for maximum tip speed in rad/s

%% HUB HEIGHT WIND PROFILE CALCULATOR
h_hub=h_0*P_rated/P_original; % scale hub height linearly with rated power based on original turbine
U_array=Speed_profile(U_0, z_0, alpha, h_0, h_hub);

%% WEIBULL REGRESSOR
f_curve=Weibull_regressor(U_array);

%% DIAMETER OPTIMIZER
LPC_0=1/E_y_original; % calculate the LPC for the original turbine
LPC=LPC_calculator(E_y,D_0,D_original,P_rated,P_original); % calculate the LPC for the assumed diameter

while (LPC-LPC_0)/LPC>0.1
    %% POWER CURVE CALCULATOR
    [P_curve,U_rated] = Power_curve(P_rated,rho,U_ci,U_co,c_p,D_0,eff,U_array);

    %% YEARLY ENERGY GENERATION CALCULATOR
    E_y=Yearly_energy(P_curve,U_ci,U_co,f_curve);

    %% LPC CALCULATOR
    LPC = LPC_calculator(E_y,D,D_original,P_rated,P_original);
end

%% MINIMUM AND MAXIMUM TIP SPEED CALCULATOR
min_tip_speed = omega_min(U_ci, lambda, D);
max_tip_speed = min(max_tip_speed_limit, tipspeed_max(U_rated, lambda_design, D));

%% TORQUE CALCULATOR
Q = torque(P_curve, U_array, D, lambda_design);

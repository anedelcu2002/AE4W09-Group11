% problems to check: 
% currently all terms are defined and calculated as arrays of integers, starting from 1

%% INPUT DATA
U_0=[1,1,2,2,3,3,3,4,4,4,4,5,5,5,6,6,7,7,8,9,10,11,12]; % measured wind speed at altitude h_0 in m/s
z_0=0.001; % surface roughness during measurement in m
h_0=10; % altitude of wind speed measurement in m
P_rated= 7.5*10^6; % selected rated power in W
rho=1.225; % air density at hub level in kg/m3
P_original= 5*10^6; % original turbine rated power in W
D_original=70; % original turbine diameter in W
h_hub_original= 90;

U_ci=5; % assumed cut-in speed in m/s
U_co=20; % assumed cut-out speed in m/s
c_p=0.5; % assumed power coefficient
alpha=0.4; % constant for power law wind speed profile
eff=0.9; % assumed turbine efficiency
lambda_design=8; % assumed design tip speed ratio
max_tip_speed_limit=100; % upper limit for maximum tip speed in rad/s

%% DIAMETER ARRAY CALCULATOR
D_array=[];
LPC_array=[];

for D=0.5*D_original*P_rated/P_original:2.0*D_original*P_rated/P_original
    %% HUB HEIGHT WIND PROFILE CALCULATOR
    h_hub=Hub_height(D, D_original, h_hub_original);
    U_array=Speed_profile(U_0, z_0, alpha, h_0, h_hub);

    %% WEIBULL REGRESSOR
    f_curve=Weibull_regressor(U_array, U_co);

    %% POWER CURVE CALCULATOR
    [P_curve,U_rated] = Power_curve(P_rated,rho,U_ci,U_co,c_p,D,eff,U_array);

    %% YEARLY ENERGY GENERATION CALCULATOR
    E_y=Yearly_energy(P_curve,U_ci,U_co,f_curve);

    %% LPC CALCULATOR
    LPC = LPC_calculator(E_y,D,D_original,P_rated,P_original);

    %% VALUE STORAGE
    D_array=[D_array D];
    LPC_array=[LPC_array LPC];
end

%% DIAMETER SELECTION
[LPC_min, i]=min(LPC_array);
D = D_array(i)

%%DESIGN VALUES BASED ON SELECTED DIAMETER
h_hub=Hub_height(D, D_original, h_hub_original)
U_array=Speed_profile(U_0, z_0, alpha, h_0, h_hub);
U_r = (P_rated/(0.5*eff*rho*c_p*pi*(D/2)^2))^(1/3)

f_curve=Weibull_regressor(U_array, U_co);
[P_curve,U_rated] = Power_curve(P_rated,rho,U_ci,U_co,c_p,D,eff,U_array);
E_y=Yearly_energy(P_curve,U_ci,U_co,f_curve)
LPC = LPC_calculator(E_y,D,D_original,P_rated,P_original)

%% MINIMUM AND MAXIMUM TIP SPEED CALCULATOR
min_tip_speed = omega_min(U_ci, lambda_design, D)*(D/2);
max_tip_speed = tipspeed_max(U_rated, lambda_design, D);
if max_tip_speed > max_tip_speed_limit
    disp('WARNING: tip speeds to high, ?lambda? needs adjustment')
    max_tip_speed = max_tip_speed_limit;
end

%% TORQUE CALCULATOR
Q = torque(P_curve, U_array, D, lambda_design, U_co);


%% PLOT RESULTS
figure;
plot(1:U_co, P_curve);
title('Power curve');
xlabel('Power [W]');
ylabel('Wind speed [m/s]');
axis tight

figure;
title('Weibull distribution');
xlabel('Wind speed (m/s)');
ylabel('Probability');
axis tight
histfit(U_array, ceil(max(U_array)), 'wbl');

figure;
plot(D_array, LPC_array);
title('LPC-diameter plot');
xlabel('Diameter [m]');
ylabel('LPC [1/J]');
axis tight

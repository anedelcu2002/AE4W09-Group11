% problems to check: 
% currently all terms are defined and calculated as arrays of integers, starting from 1

clear all;

%% INPUT DATA
U_0=[1,1,2,2,3,3,3,4,4,4,4,5,5,5,6,6,7,7,8,9,10,11,12]; % measured wind speed at altitude h_0 in m/s
z_0=0.1; % surface roughness during measurement in m
h_0=10; % altitude of wind speed measurement in m
P_rated=3.5*10^6; % selected rated power in W
rho=1.225; % air density at hub level in kg/m3
P_original=5*10^6; % original turbine rated power in W
D_original=126.5; % original turbine diameter in W
h_original=90; % original turbine hub height in m
C_original=10^6; % original turbine installation costs in EUR
L_original=20; % original turbine lifetime in years
U_ci=2; % assumed cut-in speed in m/s
U_co=25; % assumed cut-out speed in m/s
c_p=0.482; % assumed power coefficient
alpha=0.143; % constant for power law wind speed profile
eff=0.944; % assumed turbine efficiency
lambda_design=8; % assumed design tip speed ratio
rot_speed_limit=100; % upper limit for maximum tip speed in rad/s

a=4.1; % weibull scale parameter to skip weibull regressor, leave zero if not known
k=1.43; % weibull shape parameter to skip weibull regressor, leave zero if not known

%% DIAMETER ARRAY CALCULATOR
D_array=[];
LPC_array=[];
CF_array=[];

for D=(0.8*D_original*P_rated/P_original):(2*D_original*P_rated/P_original)
    if a==0 && k==0
        %% HUB HEIGHT WIND PROFILE CALCULATOR
        h_hub=h_original*D/D_original; % scale hub height linearly with rated power based on original turbine, could use different rule
        U_array=Speed_profile(U_0, z_0, alpha, h_0, h_hub);

        %% WEIBULL REGRESSOR
        f_curve=Weibull_regressor(U_array);
    else 
        %% DIRECT REGRESSOR AND SCALER
        h_hub=h_original*D/D_original; % scale hub height linearly with rated power based on original turbine, could use different rule
        wind_speed_factor=Speed_profile(1, z_0, alpha, h_0, h_hub); % shift speed to hub height
        f_curve = [];
        for i = 1:U_co*2
            f_curve(i) = (k/a)*((i/wind_speed_factor/a)^(k-1))*exp(-(i/wind_speed_factor/a)^k); % if shape and scale parameters are known, define Weibull curve by them
        end
        f_curve=f_curve/sum(f_curve);
    end

    %% POWER CURVE CALCULATOR
    [P_curve,U_rated] = Power_curve(P_rated,rho,U_ci,U_co,c_p,D,eff);

    %% YEARLY ENERGY GENERATION CALCULATOR
    E_y=Yearly_energy(P_curve,U_ci,U_co,f_curve);

    %% LPC CALCULATOR
    LPC = LPC_calculator(E_y,D,D_original,P_rated,P_original);

    %% CAPACITY FACTOR CALCULATOR
    CF=E_y/P_rated/(3600*24*365);

    %% VALUE STORAGE
    D_array=[D_array D];
    LPC_array=[LPC_array LPC];
    CF_array=[CF_array CF];
end

%% DIAMETER SELECTION
[LPC_min, i]=min(LPC_array);
D=D_array(i);

%% RECALCULATE VARIABLES WITH SELECTED DIAMETER
if a==0 && k==0
    h_hub=h_original*D/D_original; 
    U_array=Speed_profile(U_0, z_0, alpha, h_0, h_hub);
    f_curve=Weibull_regressor(U_array);
else 
    h_hub=h_original*D/D_original
    wind_speed_factor=Speed_profile(1, z_0, alpha, h_0, h_hub);
    f_curve = [];
    for i = 1:U_co*2
        f_curve(i) = (k/a)*((i*(wind_speed_factor-1)/a)^(k-1))*exp(-(i*(wind_speed_factor-1)/a)^k);
    end
    f_curve=f_curve/sum(f_curve);
end

[P_curve,U_rated] = Power_curve(P_rated,rho,U_ci,U_co,c_p,D,eff);

E_y=Yearly_energy(P_curve,U_ci,U_co,f_curve);

LPC = LPC_calculator(E_y,D,D_original,P_rated,P_original);

CF=E_y/P_rated/(3600*24*365);

%% ENERGY COST CALCULATOR
C=C_original*LPC*E_y; % installation cost for designed turbine, assuming linear scaling (?)
E_cost=C/E_y/L_original; % energy cost for designed turbine, assuming lifetime of x years (?)


%% MINIMUM AND MAXIMUM TIP SPEED CALCULATOR
min_rot_speed = omega_min(U_ci, lambda_design, D); % calculate minimum rotation speed
max_rot_speed = omega_max(U_rated, lambda_design, D); % calculate maximum rotation speed

if max_rot_speed>rot_speed_limit
    disp('Calculated maximum rotation speed is higher than the set rotation speed limit!')
end

min_tip_speed=min_rot_speed*D/2;
max_tip_speed=max_rot_speed*D/2;

%% TORQUE CALCULATOR
Q = torque(P_curve, U_co, D, lambda_design);


%% PLOT RESULTS
figure;
plot(1:U_co, P_curve);
title('Power curve');
xlabel('Power [W]');
ylabel('Wind speed [m/s]');
axis tight

figure;
if a==0 && k==0
    histfit(U_array, ceil(max(U_array)), 'wbl');
else
    plot(f_curve);
end
title('Weibull distribution');
xlabel('Wind speed (m/s)');
ylabel('Probability');
axis tight

figure;
plot(D_array, LPC_array);
title('LPC-diameter plot');
xlabel('Diameter [m]');
ylabel('LPC [1/J]');
axis tight

figure;
plot(D_array, CF_array);
title('CF-diameter plot');
xlabel('Diameter [m]');
ylabel('Capacity factor');
axis tight

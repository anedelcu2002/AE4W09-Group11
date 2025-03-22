%CONVERTER DATA
f_min = 4; %Hz
f_max = 50; %Hz
eff_con = 0.97;
M_con = 3500; %kg

%GENERATOR DATA
n = 4;
F_d = 15*10^3; %N/m^2
r_s = 0.4; %m 
eff_gen = 0.975;
rho = 7130;

%GEARBOX DATA
stage_types = [true,false,false];

%TASK 2 DATA
U_ci = 2;
P_r = 3.5*10^6; %W
P_r_gencon = P_r/eff_con;
P_r_HSS = P_r_gencon/eff_gen;
load('../Task 2 System Design/P_curve.mat')
P_min = P_curve(U_ci);
P_min_gencon = P_min/eff_con;
P_min_HSS = P_min_gencon/eff_gen;

Omega_LSS_max = 1.0327; %rad/s at rated
Omega_LSS_min = 0.224; %rad/s at cut-in


%GEARBOX
[Omega_HSS_min_allowed, Omega_HSS_max_allowed] = Omega_HSS_allowed_range(n, f_min, f_max); %rad/s
f_HSS_min_allowed = Omega_HSS_min_allowed/(2*pi);
f_HSS_max_allowed = Omega_HSS_max_allowed/(2*pi);

r_goal = gbx_ratio_total(Omega_LSS_max, Omega_LSS_min, Omega_HSS_max_allowed, Omega_HSS_min_allowed, 150)
r_list = gbx_design(r_goal, 3, stage_types);
r_total = prod(r_list)

[Omega_HSS_min, Omega_HSS_max] = Omega_HSS_used_range(r_total, Omega_LSS_min, Omega_LSS_max);
f_HSS_min = Omega_HSS_min/(2*pi);
f_HSS_max = Omega_HSS_max/(2*pi);

eff_gbx = gbx_eff(r_list,stage_types)

P_r_LSS = P_r_HSS/eff_gbx
P_min_LSS = P_min_HSS/eff_gbx;

Q_LSS_min = P_min_LSS/Omega_LSS_min %At cut-in
Q_LSS_max = P_r_LSS/Omega_LSS_max %At rated

Q_HSS_min = P_min_HSS/Omega_HSS_min %At cut-in
Q_HSS_max = P_r_HSS/Omega_HSS_max %At rated

M_gbx = gbx_mass(Q_LSS_max) %Q_RATED


%GENERATOR 
V_gen = gen_sizing(F_d, P_r, Omega_HSS_max)
l_s = gen_length(V_gen, r_s)

[M_gen, J_gen] = gen_mass_inertia(r_s, l_s, rho)

M_nacel = M_con + M_gen + M_gbx
eff_total = eff_con * eff_gen * eff_gbx


